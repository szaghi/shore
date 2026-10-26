"""OCCT backend of the ``Geometry`` protocol (ADR 0002).

Exact B-rep geometry from STEP, BREP or IGES via OpenCASCADE (OCP, installed
with build123d: ``pip install 'shore-mesh[cad]'``). Entities are the shape's
faces (``SURFACE``), edges (``CURVE``) and vertices (``POINT``), indexed in
``TopExp.MapShapes`` order; indices are stable for one load only (persistent
identity is ADR 0003).

Projection has two tiers:

- ``exact=False`` (tier 1, vectorised): closest point on a tessellation that
  is built once at load and keeps (u, v) and the oriented surface normal at
  every node; parameters and normals are interpolated barycentrically.
  Closest points come from trimesh on the tessellation (R-trees built once at
  load, ~2.5e4-1.3e5 points/s). BRepMesh enforces its deflection on edges;
  triangle interiors sag up to ~2x further. For a query at distance ``h``
  from the surface the error is about ``2 * deflection + h * (facet normal
  error)``: a couple of deflections for points near the surface (smoothing,
  dragging), more for distant points, which should use ``exact=True``.
  Where a parameter is singular (a sphere pole, a cone apex: the surface
  derivative along it vanishes) the arbitrary value at the singular node is
  left out of that parameter's interpolation.
- ``exact=True`` (tier 2, per point): Newton projection onto the exact face
  surface, seeded from tier 1, with a trimming check; points whose foot lies
  outside the trimmed face fall back to an exact shape-distance query, so the
  answer stays on the face (or its boundary). Curves use a clamped Newton
  iteration on the exact edge curve.

``ray_cast`` answers on the tessellation (tier 1), through the same
``Projector`` the STL backend uses.

These calls run in the caller's thread and hold the GIL; a GUI must run long
loads in a subprocess (ADR 0006).
"""

from __future__ import annotations

import functools
import logging
from pathlib import Path

import numpy as np
import trimesh
import trimesh.proximity
import trimesh.triangles
from OCP.Bnd import Bnd_Box
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex, BRepBuilderAPI_Sewing
from OCP.BRepClass import BRepClass_FaceClassifier
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepGProp import BRepGProp
from OCP.BRepLib import BRepLib_ToolTriangulatedShape
from OCP.BRepLProp import BRepLProp_SLProps
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.collections import IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher as _ShapeMap
from OCP.GCPnts import GCPnts_QuasiUniformDeflection
from OCP.gp import gp_Pnt, gp_Pnt2d, gp_Vec
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.IGESControl import IGESControl_Reader
from OCP.ShapeAnalysis import ShapeAnalysis_Surface
from OCP.ShapeFix import ShapeFix_Shape
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_OUT, TopAbs_REVERSED, TopAbs_VERTEX
from OCP.TopExp import TopExp
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS, TopoDS_Shape

from shore.exceptions import MeshInputError, ParameterError
from shore.geometry.base import BODY, EntityKind, EntityRef, Projection
from shore.geometry.keys import EntityKey, make_key
from shore.geometry.stl import StlGeometry
from shore.surface.projector import ProjectionResult

log = logging.getLogger(__name__)

_NEWTON_ITERS = 20


def _xyz(p: gp_Pnt) -> tuple[float, float, float]:
    return p.X(), p.Y(), p.Z()


class OcctGeometry:
    """B-rep body as a ``Geometry``.

    Parameters
    ----------
    shape :
        OCCT shape (``TopoDS_Shape``), e.g. ``build123d_object.wrapped``.
    deflection :
        Absolute linear deflection of the tier-1 tessellation. Default:
        ``rel_deflection`` times the bounding-box diagonal.
    rel_deflection :
        Linear deflection relative to the bounding-box diagonal.
    angular_deflection :
        Angular deflection of the tessellation, radians.
    heal :
        Run ``ShapeFix_Shape`` on the input first.
    """

    def __init__(
        self,
        shape: TopoDS_Shape,
        *,
        deflection: float | None = None,
        rel_deflection: float = 1e-3,
        angular_deflection: float = 0.2,
        heal: bool = True,
    ) -> None:
        if shape is None or shape.IsNull():
            raise MeshInputError("OcctGeometry: empty shape.")
        if heal:
            fixer = ShapeFix_Shape(shape)
            fixer.Perform()
            shape = fixer.Shape()
        self.shape = shape

        box = Bnd_Box()
        BRepBndLib.AddOptimal_s(shape, box, False, False)
        lo, hi = box.CornerMin(), box.CornerMax()
        self._bounds = np.array([_xyz(lo), _xyz(hi)], dtype=np.float64)
        self.diag = float(np.linalg.norm(self._bounds[1] - self._bounds[0]))
        self.deflection = float(deflection) if deflection else rel_deflection * self.diag
        self.angular_deflection = angular_deflection
        self.tol = 1e-9 * self.diag  # 3-D precision for exact queries

        self.faces = self._map(TopAbs_FACE, TopoDS.Face)
        self.edges = self._map(TopAbs_EDGE, TopoDS.Edge)
        self.vertices = self._map(TopAbs_VERTEX, TopoDS.Vertex)
        if not self.faces:
            raise MeshInputError("OcctGeometry: shape has no faces.")

        self._tessellate()
        self._edge_polys = [self._sample_edge(e) for e in self.edges]
        self._face_proxies: dict[int, tuple[trimesh.Trimesh, np.ndarray]] = {}
        self._surf_tools: dict[int, tuple[ShapeAnalysis_Surface, BRepAdaptor_Surface]] = {}
        self._keys: dict[EntityRef, EntityKey] = {}

    # -- construction ----------------------------------------------------------

    @classmethod
    def from_file(cls, path: str | Path, **kwargs: object) -> OcctGeometry:
        """Load ``.step``/``.stp`` (build123d), ``.brep``/``.brp`` (build123d)
        or ``.iges``/``.igs`` (OCP ``IGESControl_Reader``)."""
        path = Path(path)
        if not path.exists():
            raise MeshInputError(f"Geometry file not found: {path}")
        suffix = path.suffix.lower()
        if suffix in (".step", ".stp", ".brep", ".brp"):
            import build123d as bd

            reader = bd.import_step if suffix in (".step", ".stp") else bd.import_brep
            try:
                shape = reader(str(path)).wrapped
            except Exception as err:  # build123d raises bare ValueError/RuntimeError
                raise MeshInputError(f"{path.name}: cannot read CAD file ({err}).") from err
        elif suffix in (".iges", ".igs"):
            reader = IGESControl_Reader()
            if reader.ReadFile(str(path)) != IFSelect_RetDone:
                raise MeshInputError(f"{path.name}: cannot read IGES file.")
            reader.TransferRoots()
            # IGES stores faces as independent trimmed surfaces: sew them so
            # that adjacent faces share edges (needed for association)
            sewing = BRepBuilderAPI_Sewing(1e-6)
            sewing.Add(reader.OneShape())
            sewing.Perform()
            shape = sewing.SewedShape()
        else:
            raise MeshInputError(f"{path.name}: not a CAD file (STEP, BREP or IGES).")
        return cls(shape, **kwargs)  # type: ignore[arg-type]

    def _map(self, kind: object, cast: object) -> list:
        shapes = _ShapeMap()
        TopExp.MapShapes_s(self.shape, kind, shapes)
        return [cast(shapes.FindKey(i)) for i in range(1, shapes.Extent() + 1)]

    def _tessellate(self) -> None:
        """Tier-1 proxy: per-face nodes with (u, v) and oriented normals."""
        BRepMesh_IncrementalMesh(self.shape, self.deflection, False, self.angular_deflection, True)
        verts, uvs, normals, tris, owner, singular = [], [], [], [], [], []
        base = 0
        for k, face in enumerate(self.faces):
            loc = TopLoc_Location()
            tri = BRep_Tool.Triangulation_s(face, loc)
            if tri is None:
                log.warning("face %d has no triangulation; skipped in the proxy", k)
                continue
            BRepLib_ToolTriangulatedShape.ComputeNormals_s(face, tri)
            trsf = loc.Transformation()
            sign = -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0
            adaptor = BRepAdaptor_Surface(face)
            pnt, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
            n = tri.NbNodes()
            v = np.empty((n, 3))
            uv = np.empty((n, 2))
            nn = np.empty((n, 3))
            sing = np.zeros((n, 2), dtype=bool)
            for i in range(n):
                v[i] = _xyz(tri.Node(i + 1).Transformed(trsf))
                p2 = tri.UVNode(i + 1)
                uv[i] = (p2.X(), p2.Y())
                d = gp_Vec(tri.Normal(i + 1)).Transformed(trsf)
                nn[i] = (sign * d.X(), sign * d.Y(), sign * d.Z())
                adaptor.D1(uv[i, 0], uv[i, 1], pnt, du, dv)
                sing[i] = (du.Magnitude() <= self.tol, dv.Magnitude() <= self.tol)
            t = np.array([tri.Triangle(i + 1).Get() for i in range(tri.NbTriangles())]) - 1
            if sign < 0.0:
                t = t[:, ::-1]
            verts.append(v)
            uvs.append(uv)
            normals.append(nn)
            singular.append(sing)
            tris.append(t + base)
            owner.append(np.full(len(t), k))
            base += n
        self._verts = np.concatenate(verts)
        self._uv = np.concatenate(uvs)
        self._uv_singular = np.concatenate(singular)  # (n, 2): u / v singular at node
        self._normals = np.concatenate(normals)
        self._tris = np.concatenate(tris)
        self._tri_face = np.concatenate(owner)
        self._proxy = trimesh.Trimesh(self._verts, self._tris, process=False)
        self._proxy_geometry = StlGeometry(self._proxy)
        _ = self._proxy.triangles_tree  # build the R-tree once, at load

    def _sample_edge(self, edge: object) -> tuple[np.ndarray, np.ndarray]:
        """Tier-1 polyline of an edge: points (n, 3) and parameters (n,)."""
        curve = BRepAdaptor_Curve(edge)
        t0, t1 = curve.FirstParameter(), curve.LastParameter()
        if BRep_Tool.Degenerated_s(edge):
            p = _xyz(BRep_Tool.Pnt_s(TopoDS.Vertex(_first_vertex(edge))))
            return np.array([p, p]), np.array([t0, t1])
        sampler = GCPnts_QuasiUniformDeflection(curve, self.deflection)
        if not sampler.IsDone() or sampler.NbPoints() < 2:
            ts = np.linspace(t0, t1, 33)
        else:
            ts = np.array([sampler.Parameter(i) for i in range(1, sampler.NbPoints() + 1)])
        pts = np.array([_xyz(curve.Value(t)) for t in ts])
        return pts, ts

    # -- Geometry protocol -----------------------------------------------------

    @property
    def bounds(self) -> np.ndarray:
        return self._bounds.copy()

    def entities(self, kind: EntityKind) -> tuple[EntityRef, ...]:
        n = {
            EntityKind.BODY: 1,
            EntityKind.SURFACE: len(self.faces),
            EntityKind.CURVE: len(self.edges),
            EntityKind.POINT: len(self.vertices),
        }[kind]
        return tuple(EntityRef(kind, i) for i in range(n))

    def ray_cast(
        self, anchor: np.ndarray, dirs: np.ndarray, on_failure: str = "raise"
    ) -> ProjectionResult:
        return self._proxy_geometry.ray_cast(anchor, dirs, on_failure)

    def project(
        self,
        points: np.ndarray,
        entity: EntityRef | None = None,
        *,
        exact: bool = False,
    ) -> Projection:
        pts = np.asarray(points, dtype=np.float64)
        shape = pts.shape[:-1]
        flat = pts.reshape(-1, 3)
        if entity is None or entity == BODY:
            proj = self._project_body(flat, exact)
        elif entity.kind is EntityKind.SURFACE:
            proj = self._project_face(flat, self._check(entity), exact)
        elif entity.kind is EntityKind.CURVE:
            proj = self._project_edge(flat, self._check(entity), exact)
        elif entity.kind is EntityKind.POINT:
            p = np.array(_xyz(BRep_Tool.Pnt_s(self.vertices[self._check(entity)])))
            pp = np.broadcast_to(p, flat.shape).copy()
            proj = Projection(pp, np.linalg.norm(flat - pp, axis=-1), np.full_like(pp, np.nan))
        else:
            raise ParameterError(f"Unsupported entity {entity}.")
        return _reshape(proj, shape)

    def evaluate(self, entity: EntityRef, params: np.ndarray) -> np.ndarray:
        k = self._check(entity)
        prm = np.asarray(params, dtype=np.float64)
        if entity.kind is EntityKind.SURFACE:
            _, adaptor = self._surface_tools(k)
            flat = prm.reshape(-1, 2)
            out = np.array([_xyz(adaptor.Value(u, v)) for u, v in flat])
            return out.reshape(*prm.shape[:-1], 3)
        if entity.kind is EntityKind.CURVE:
            curve = BRepAdaptor_Curve(self.edges[k])
            flat = prm.reshape(-1)
            out = np.array([_xyz(curve.Value(t)) for t in flat])
            return out.reshape(*prm.shape, 3)
        if entity.kind is EntityKind.POINT:
            p = np.array(_xyz(BRep_Tool.Pnt_s(self.vertices[k])))
            return np.broadcast_to(p, (*prm.shape, 3)).copy()
        raise ParameterError("The body has no parametrisation; evaluate a SURFACE or CURVE.")

    def fingerprint(self, entity: EntityRef) -> EntityKey:
        """Fingerprint: OCCT geometric type, centroid, area/length, bbox (cached)."""
        if entity == BODY:
            shape, kind, gtype = self.shape, EntityKind.BODY, "brep"
        else:
            k = self._check(entity)
            kind = entity.kind
            if kind is EntityKind.SURFACE:
                shape = self.faces[k]
                gtype = BRepAdaptor_Surface(shape).GetType().name.removeprefix("GeomAbs_")
            elif kind is EntityKind.CURVE:
                shape = self.edges[k]
                gtype = (
                    "Degenerate"
                    if BRep_Tool.Degenerated_s(shape)
                    else BRepAdaptor_Curve(shape).GetType().name.removeprefix("GeomAbs_")
                )
            else:
                shape, gtype = self.vertices[k], "Point"
        if entity in self._keys:
            return self._keys[entity]
        props = GProp_GProps()
        if kind in (EntityKind.BODY, EntityKind.SURFACE):
            BRepGProp.SurfaceProperties_s(shape, props)
        elif kind is EntityKind.CURVE and gtype != "Degenerate":
            BRepGProp.LinearProperties_s(shape, props)
        if kind is EntityKind.POINT or (kind is EntityKind.CURVE and gtype == "Degenerate"):
            vertex = shape if kind is EntityKind.POINT else TopoDS.Vertex(_first_vertex(shape))
            centroid, measure = np.array(_xyz(BRep_Tool.Pnt_s(vertex))), 0.0
        else:
            centroid, measure = np.array(_xyz(props.CentreOfMass())), props.Mass()
        box = Bnd_Box()
        BRepBndLib.AddOptimal_s(shape, box, False, False)
        bbox = np.array([_xyz(box.CornerMin()), _xyz(box.CornerMax())])
        key = make_key(kind, gtype, centroid, measure, bbox, self._bounds)
        self._keys[entity] = key
        return key

    # -- OCCT specifics ----------------------------------------------------------

    @property
    def proxy(self) -> trimesh.Trimesh:
        """The tier-1 tessellation (faces not merged across edges)."""
        return self._proxy

    @functools.cached_property
    def merged_proxy(self) -> trimesh.Trimesh:
        """The tessellation with coincident vertices merged and degenerate
        triangles dropped: watertight for a closed solid.

        For mesh-based consumers that need a closed triangulation
        (inside/outside tests, plane sections, dihedral-angle TE detection).
        Its vertex normals are smoothed across face boundaries, so do not
        read normals from it: use ``project``.
        """
        mesh = trimesh.Trimesh(self._verts, self._tris, process=True)
        # singular edges (sphere poles, fillet corner patches) are meshed with
        # several coincident nodes; merged, their fan triangles collapse and
        # leave non-manifold edges until removed
        mesh.update_faces(mesh.nondegenerate_faces())
        mesh.remove_unreferenced_vertices()
        return mesh

    def _check(self, entity: EntityRef) -> int:
        n = len(self.entities(entity.kind))
        if not 0 <= entity.index < n:
            raise ParameterError(
                f"{entity} out of range: this geometry has {n} {entity.kind.value}s."
            )
        return entity.index

    def _surface_tools(self, k: int) -> tuple[ShapeAnalysis_Surface, BRepAdaptor_Surface]:
        if k not in self._surf_tools:
            face = self.faces[k]
            self._surf_tools[k] = (
                ShapeAnalysis_Surface(BRep_Tool.Surface_s(face)),
                BRepAdaptor_Surface(face),
            )
        return self._surf_tools[k]

    def _face_proxy(self, k: int) -> tuple[trimesh.Trimesh, np.ndarray]:
        """Tessellation of face *k* alone, with its global triangle ids (lazy)."""
        if k not in self._face_proxies:
            ids = np.flatnonzero(self._tri_face == k)
            mesh = trimesh.Trimesh(self._verts, self._tris[ids], process=False)
            _ = mesh.triangles_tree
            self._face_proxies[k] = (mesh, ids)
        return self._face_proxies[k]

    def _interpolate(
        self, mesh: trimesh.Trimesh, tri_ids: np.ndarray | None, flat: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Tier 1: points, distances, normals, uv, global triangle ids."""
        closest, dist, local = trimesh.proximity.closest_point(mesh, flat)
        tid = local if tri_ids is None else tri_ids[local]
        corners = self._tris[tid]
        tv = self._verts[corners]
        with np.errstate(divide="ignore", invalid="ignore"):
            bary = trimesh.triangles.points_to_barycentric(tv, closest)
        bad = ~np.isfinite(bary).all(axis=1)
        if bad.any():  # zero-area triangle (e.g. a pole fan): weight the nearest corner
            near = np.linalg.norm(tv[bad] - closest[bad, None], axis=-1).argmin(axis=1)
            bary[bad] = np.eye(3)[near]
        normals = (bary[..., None] * self._normals[corners]).sum(axis=1)
        mag = np.linalg.norm(normals, axis=-1, keepdims=True)
        normals = normals / np.where(mag < 1e-15, 1.0, mag)
        # per parameter, drop corners where it is singular (pole, apex) if any
        # regular corner remains; their value there is arbitrary
        w = bary[..., None] * ~self._uv_singular[corners]  # (n, 3, 2)
        wsum = w.sum(axis=1, keepdims=True)
        w = np.where(wsum > 0.0, w / np.where(wsum > 0.0, wsum, 1.0), bary[..., None])
        uv = (w * self._uv[corners]).sum(axis=1)
        return closest, np.asarray(dist), normals, uv, tid

    def _project_body(self, flat: np.ndarray, exact: bool) -> Projection:
        pts, dist, normals, uv, tri_ids = self._interpolate(self._proxy, None, flat)
        owner = self._tri_face[tri_ids]
        if exact:
            for i in range(len(flat)):
                pts[i], uv[i], normals[i] = self._exact_on_face(
                    int(owner[i]), flat[i], uv[i], normals[i]
                )
            dist = np.linalg.norm(flat - pts, axis=-1)
        return Projection(pts, dist, normals, params=None, entity_index=owner)

    def _project_face(self, flat: np.ndarray, k: int, exact: bool) -> Projection:
        mesh, ids = self._face_proxy(k)
        pts, dist, normals, uv, _ = self._interpolate(mesh, ids, flat)
        if exact:
            for i in range(len(flat)):
                pts[i], uv[i], normals[i] = self._exact_on_face(k, flat[i], uv[i], normals[i])
            dist = np.linalg.norm(flat - pts, axis=-1)
        return Projection(pts, dist, normals, params=uv)

    def _exact_on_face(
        self, k: int, p: np.ndarray, uv0: np.ndarray, n0: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Exact closest point on trimmed face *k*: (point, uv, normal)."""
        sas, adaptor = self._surface_tools(k)
        face = self.faces[k]
        gp = gp_Pnt(*p)
        uv = sas.NextValueOfUV(gp_Pnt2d(float(uv0[0]), float(uv0[1])), gp, self.tol)
        tol = max(BRep_Tool.Tolerance_s(face), self.tol)
        if BRepClass_FaceClassifier(face, uv, tol).State() == TopAbs_OUT:
            # foot outside the trimmed face: exact distance to the face (lands on its boundary)
            dss = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp).Vertex(), face)
            q = dss.PointOnShape2(1)
            uv = sas.NextValueOfUV(uv, q, self.tol)
            point = np.array(_xyz(q))
        else:
            point = np.array(_xyz(adaptor.Value(uv.X(), uv.Y())))
        props = BRepLProp_SLProps(adaptor, uv.X(), uv.Y(), 1, self.tol)
        if props.IsNormalDefined():
            d = props.Normal()
            sign = -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0
            normal = sign * np.array([d.X(), d.Y(), d.Z()])
        else:  # singular point (e.g. a pole): keep the tier-1 normal
            normal = n0
        return point, np.array([uv.X(), uv.Y()]), normal

    def _project_edge(self, flat: np.ndarray, k: int, exact: bool) -> Projection:
        poly, ts = self._edge_polys[k]
        pts, ts_proj = _closest_on_polyline(flat, poly, ts)
        if exact and not BRep_Tool.Degenerated_s(self.edges[k]):
            curve = BRepAdaptor_Curve(self.edges[k])
            t0, t1 = curve.FirstParameter(), curve.LastParameter()
            for i in range(len(flat)):
                ts_proj[i] = _newton_on_curve(curve, flat[i], ts_proj[i], t0, t1, self.tol)
                pts[i] = _xyz(curve.Value(ts_proj[i]))
        dist = np.linalg.norm(flat - pts, axis=-1)
        return Projection(pts, dist, np.full_like(pts, np.nan), params=ts_proj)


def _first_vertex(edge: object) -> object:
    return TopExp.FirstVertex_s(edge)


def _reshape(proj: Projection, shape: tuple[int, ...]) -> Projection:
    params = proj.params
    if params is not None:
        params = params.reshape(*shape, *params.shape[1:])
    owner = proj.entity_index
    return Projection(
        points=proj.points.reshape(*shape, 3),
        distances=proj.distances.reshape(shape),
        normals=proj.normals.reshape(*shape, 3),
        params=params,
        entity_index=None if owner is None else owner.reshape(shape),
    )


def _closest_on_polyline(
    pts: np.ndarray, poly: np.ndarray, ts: np.ndarray, chunk: int = 4096
) -> tuple[np.ndarray, np.ndarray]:
    """Closest point on a polyline and its interpolated parameter (vectorised)."""
    a, b = poly[:-1], poly[1:]
    ab = b - a
    ab2 = np.einsum("ij,ij->i", ab, ab)
    out_p = np.empty_like(pts)
    out_t = np.empty(len(pts))
    for s in range(0, len(pts), chunk):
        p = pts[s : s + chunk, None, :]
        w = np.einsum("nsj,sj->ns", p - a[None], ab) / np.where(ab2 > 0.0, ab2, 1.0)
        w = np.clip(w, 0.0, 1.0)
        foot = a[None] + w[..., None] * ab[None]
        d2 = np.einsum("nsj,nsj->ns", p - foot, p - foot)
        j = np.argmin(d2, axis=1)
        r = np.arange(len(j))
        out_p[s : s + chunk] = foot[r, j]
        out_t[s : s + chunk] = ts[j] + w[r, j] * (ts[j + 1] - ts[j])
    return out_p, out_t


def _newton_on_curve(
    curve: BRepAdaptor_Curve, p: np.ndarray, t: float, t0: float, t1: float, tol: float
) -> float:
    """Minimise |C(t) - p|^2 by Newton, clamped to [t0, t1]."""
    c, d1, d2 = gp_Pnt(), gp_Vec(), gp_Vec()
    for _ in range(_NEWTON_ITERS):
        curve.D2(t, c, d1, d2)
        r = np.array(_xyz(c)) - p
        v1 = np.array([d1.X(), d1.Y(), d1.Z()])
        v2 = np.array([d2.X(), d2.Y(), d2.Z()])
        g = r @ v1
        h = v1 @ v1 + r @ v2
        step = g / h if h > 0.0 else g / max(v1 @ v1, 1e-300)
        t_new = min(max(t - step, t0), t1)
        if abs(t_new - t) * np.sqrt(v1 @ v1) < tol:
            return t_new
        t = t_new
    return t
