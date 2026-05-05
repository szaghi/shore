"""OOP mesh layer for SHORE."""

from shore.mesh.block import Block, HexBlock
from shore.mesh.body import Body
from shore.mesh.edge import Edge
from shore.mesh.face import AxisMap, Face, FaceBC
from shore.mesh.inspect import (
    check_geo,
    export_geo,
    export_surface,
    info_geo,
    info_surface,
    view_geo_file,
    view_surface_file,
)
from shore.mesh.mesh import Mesh
from shore.mesh.spacing import Spacing1D
from shore.mesh.topology import (
    CGridTipCap,
    CGridTopology,
    CHTopology,
    CubedSphereTopology,
    OGridTopology,
    OHTopology,
    Topology,
)
from shore.mesh.vertex import Vertex

__all__ = [
    "AxisMap",
    "Block",
    "Body",
    "CGridTipCap",
    "CGridTopology",
    "CHTopology",
    "CubedSphereTopology",
    "Edge",
    "Face",
    "FaceBC",
    "HexBlock",
    "Mesh",
    "OGridTopology",
    "OHTopology",
    "Spacing1D",
    "Topology",
    "Vertex",
    "check_geo",
    "export_geo",
    "export_surface",
    "info_geo",
    "info_surface",
    "view_geo_file",
    "view_surface_file",
]
