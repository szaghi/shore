"""Generate the Python API reference of the docs from the package's docstrings.

Writes one Markdown page per public module under ``docs/reference/api/``,
an index page and ``_sidebar.json`` (read by ``docs/.vitepress/config.mts``).
The source is parsed with :mod:`ast` only: nothing is imported, so optional
extras (CAD, GUI) need not be installed.

Run from the project root::

    python docs/scripts/api_docs.py          # (re)write the pages
    python docs/scripts/api_docs.py --check  # exit 1 if the pages are stale
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "shore"
OUT = ROOT / "docs" / "reference" / "api"
REPO = "https://github.com/szaghi/shore/blob/main/src/shore"

SKIP = {"shore.cli", "shore.__main__"}
SECTIONS = (
    "Parameters",
    "Other Parameters",
    "Attributes",
    "Returns",
    "Yields",
    "Raises",
    "Warns",
    "See Also",
    "Notes",
    "References",
    "Examples",
)
LIST_SECTIONS = set(SECTIONS[:7])
GROUPS = (
    ("Package", ("shore", "shore.exceptions", "shore.quality", "shore.pipeline", "shore.config")),
    ("Geometry", ("shore.geometry", "shore.shapes", "shore.profiles")),
    ("Block model", ("shore.mesh",)),
    ("Surface", ("shore.surface",)),
    ("Volume and topologies", ("shore.volume", "shore.topology", "shore.primitive")),
    ("Collars", ("shore.collar",)),
    ("Projects", ("shore.project",)),
    ("Chimera and Xall", ("shore.io", "shore.split", "shore.balance", "shore.boxes")),
    ("GUI and visualisation", ("shore.gui", "shore.vis")),
)


@dataclass
class Page:
    name: str
    path: Path
    tree: ast.Module
    objects: list[str] = field(default_factory=list)


def module_name(path: Path) -> str:
    parts = path.relative_to(SRC.parent).with_suffix("").parts
    return ".".join(parts[:-1] if parts[-1] == "__init__" else parts)


def is_private(name: str) -> bool:
    return any(p.startswith("_") for p in name.split("."))


def exported(tree: ast.Module) -> list[str] | None:
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "__all__" for t in node.targets
        ):
            return [e.value for e in node.value.elts if isinstance(e, ast.Constant)]
    return None


def definitions(tree: ast.Module) -> dict[str, ast.AST]:
    out = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            out[node.name] = node
    return out


def imports(tree: ast.Module, package: str) -> dict[str, str]:
    """Name -> dotted source module of ``from x import name`` statements."""
    out = {}
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                parent = package.split(".")
                parent = parent[: len(parent) - node.level + 1]
                base = ".".join([*parent, base]) if base else ".".join(parent)
            for alias in node.names:
                out[alias.asname or alias.name] = base
    return out


# ---------------------------------------------------------------- docstrings


def inline(text: str) -> str:
    """RST inline markup -> Markdown; escape what Vue or HTML would eat."""
    parts = re.split(r"(``.+?``|`[^`]+`)", text)
    out = []
    for p in parts:
        if p.startswith("``") and p.endswith("``") and len(p) > 4:
            out.append("`" + p[2:-2] + "`")
        elif p.startswith("`") and p.endswith("`") and len(p) > 2:
            out.append("`" + p[1:-1].lstrip("~") + "`")
        else:
            p = re.sub(r":(?:func|class|meth|mod|attr|data|exc|obj):`~?([^`]+)`", r"`\1`", p)
            p = re.sub(r":math:`([^`]+)`", r"$\1$", p)
            p = p.replace("\\ ", "").replace("<", "&lt;").replace(">", "&gt;")
            p = p.replace("{{", "&#123;&#123;").replace("}}", "&#125;&#125;")
            p = p.replace("*", "\\*") if re.search(r"\w\*\w|\s\*\s", p) else p
            out.append(p)
    return "".join(out)


def dedent(lines: list[str]) -> list[str]:
    ind = [len(x) - len(x.lstrip()) for x in lines if x.strip()]
    cut = min(ind) if ind else 0
    return [x[cut:] for x in lines]


def split_sections(doc: str) -> tuple[list[str], list[tuple[str, list[str]]]]:
    lines = doc.expandtabs().splitlines()
    head, sections, cur, name = [], [], None, None
    i = 0
    while i < len(lines):
        line = lines[i]
        nxt = lines[i + 1] if i + 1 < len(lines) else ""
        title = line.strip()
        if title in SECTIONS and set(nxt.strip()) == {"-"} and len(nxt.strip()) >= 3:
            if name is not None:
                sections.append((name, cur))
            name, cur = title, []
            i += 2
            continue
        (head if name is None else cur).append(line)
        i += 1
    if name is not None:
        sections.append((name, cur))
    return head, sections


def prose(lines: list[str]) -> str:
    """Paragraphs, bullet lists, literal blocks (``::``) and doctests."""
    out, code, doctest = [], False, False
    lines = dedent(lines)
    i = 0
    while i < len(lines):
        line = lines[i]
        s = line.strip()
        if s.startswith(">>>") and not code:
            if not doctest:
                out.append("```python")
                doctest = True
            out.append(line.strip())
        elif doctest and s and not s.startswith(">>>"):
            out.append(s if s.startswith("...") else "# " + s)
        elif doctest and not s:
            out.append("```")
            doctest = False
            out.append("")
        elif code:
            if s and not line.startswith(" "):
                out.append("```")
                code = False
                continue
            out.append(line[4:] if line.startswith("    ") else line.lstrip())
        else:
            if s.endswith("::"):
                text = s[:-2].rstrip()
                if text:
                    out.append(inline(text + (":" if not text.endswith(":") else "")))
                out.append("")
                out.append("```")
                code = True
                i += 1
                while i < len(lines) and not lines[i].strip():
                    i += 1
                continue
            out.append(inline(line))
        i += 1
    if code or doctest:
        while out and not out[-1].strip():
            out.pop()
        out.append("```")
    text = "\n".join(out).strip()
    return re.sub(r"\n{3,}", "\n\n", text)


def items(lines: list[str]) -> str:
    """``name : type`` entries with indented descriptions -> a list."""
    lines = dedent(lines)
    entries: list[tuple[str, list[str]]] = []
    for line in lines:
        if line and not line.startswith(" "):
            entries.append((line.strip(), []))
        elif entries:
            entries[-1][1].append(line)
        elif line.strip():
            entries.append(("", [line]))
    out = []
    for term, desc in entries:
        body = " ".join(x.strip() for x in desc if x.strip())
        if " : " in term or term.endswith(" :"):
            name, _, typ = term.partition(" :")
            head = f"**`{name.strip()}`**"
            if typ.strip():
                head += f" ({inline(typ.strip())})"
        elif term:
            head = f"`{term}`"
        else:
            head = ""
        out.append(f"- {head}{' — ' if head and body else ''}{inline(body)}")
    return "\n".join(out)


def render_doc(doc: str | None) -> str:
    if not doc:
        return ""
    head, sections = split_sections(doc)
    out = [prose(head)]
    for name, lines in sections:
        if not any(x.strip() for x in lines):
            continue
        out.append(f"**{name}**")
        out.append(items(lines) if name in LIST_SECTIONS else prose(lines))
    return "\n\n".join(x for x in out if x)


# ---------------------------------------------------------------- signatures


def signature(fn: ast.FunctionDef, *, drop_self: bool) -> str:
    args = fn.args
    pos = [*args.posonlyargs, *args.args]
    if drop_self and pos and pos[0].arg in ("self", "cls"):
        pos = pos[1:]
        n_pos_only = max(len(args.posonlyargs) - 1, 0)
    else:
        n_pos_only = len(args.posonlyargs)
    defaults = [None] * (len(pos) - len(args.defaults)) + list(args.defaults)
    parts = []
    for k, (a, d) in enumerate(zip(pos, defaults[-len(pos) :] if pos else [], strict=False)):
        parts.append(_arg(a, d))
        if n_pos_only and k == n_pos_only - 1:
            parts.append("/")
    if args.vararg:
        parts.append("*" + _arg(args.vararg, None))
    elif args.kwonlyargs:
        parts.append("*")
    for a, d in zip(args.kwonlyargs, args.kw_defaults, strict=True):
        parts.append(_arg(a, d))
    if args.kwarg:
        parts.append("**" + _arg(args.kwarg, None))
    ret = f" -> {ast.unparse(fn.returns)}" if fn.returns is not None else ""
    one = f"{fn.name}({', '.join(parts)}){ret}"
    if len(one) <= 88:
        return one
    return f"{fn.name}(\n" + "".join(f"    {p},\n" for p in parts) + f"){ret}"


def _arg(a: ast.arg, default: ast.AST | None) -> str:
    s = a.arg
    if a.annotation is not None:
        s += f": {ast.unparse(a.annotation)}"
        if default is not None:
            s += f" = {ast.unparse(default)}"
    elif default is not None:
        s += f"={ast.unparse(default)}"
    return s


def is_dataclass(node: ast.ClassDef) -> bool:
    for d in node.decorator_list:
        target = d.func if isinstance(d, ast.Call) else d
        if ast.unparse(target).split(".")[-1] == "dataclass":
            return True
    return False


def decorators(node: ast.AST) -> set[str]:
    targets = (d.func if isinstance(d, ast.Call) else d for d in node.decorator_list)
    return {ast.unparse(t).split(".")[-1] for t in targets}


# ---------------------------------------------------------------- pages


def anchor(name: str) -> str:
    return name.replace(".", "-")


def source_link(rel: str, node: ast.AST) -> str:
    return f'<a class="api-source" href="{REPO}/{rel}#L{node.lineno}">source</a>'


def render_function(fn: ast.FunctionDef, rel: str, level: int, qual: str, kind: str) -> list[str]:
    tags = decorators(fn)
    label = "property" if "property" in tags else kind
    if "staticmethod" in tags:
        label = "staticmethod"
    elif "classmethod" in tags:
        label = "classmethod"
    badge = f'<Badge type="info" text="{label}" />'
    out = [f"{'#' * level} `{qual}` {badge} {{#{anchor(qual)}}}", ""]
    if "property" in tags:
        ret = f": {ast.unparse(fn.returns)}" if fn.returns is not None else ""
        out += ["```python", f"{fn.name}{ret}", "```"]
    else:
        out += ["```python", signature(fn, drop_self=kind == "method"), "```"]
    out.append(source_link(rel, fn))
    doc = render_doc(ast.get_docstring(fn))
    if doc:
        out += ["", doc]
    return [*out, ""]


def render_class(node: ast.ClassDef, rel: str, qual: str) -> list[str]:
    bases = ", ".join(ast.unparse(b) for b in node.bases)
    kind = "dataclass" if is_dataclass(node) else "class"
    out = [f'## `{qual}` <Badge type="tip" text="{kind}" /> {{#{anchor(qual)}}}', ""]
    methods = [n for n in node.body if isinstance(n, ast.FunctionDef)]
    init = next((n for n in methods if n.name == "__init__"), None)
    if init is not None:
        sig = signature(init, drop_self=True).replace("__init__", node.name, 1)
        sig = re.sub(r" -> None$", "", sig)
    elif is_dataclass(node):
        fields = [
            n for n in node.body if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name)
        ]
        parts = []
        for f in fields:
            if f.target.id.startswith("_") or "ClassVar" in ast.unparse(f.annotation):
                continue
            default = f" = {ast.unparse(f.value)}" if f.value is not None else ""
            parts.append(f"{f.target.id}: {ast.unparse(f.annotation)}{default}")
        one = f"{node.name}({', '.join(parts)})"
        many = f"{node.name}(\n" + "".join(f"    {p},\n" for p in parts) + ")"
        sig = one if len(one) <= 88 else many
    else:
        sig = f"class {node.name}({bases})" if bases else f"class {node.name}"
    out += ["```python", sig, "```"]
    if bases and (init is not None or is_dataclass(node)):
        out.append(f"Bases: `{bases}`<br>")
    out.append(source_link(rel, node))
    doc = render_doc(ast.get_docstring(node))
    if doc:
        out += ["", doc]
    out.append("")
    if init is not None and ast.get_docstring(init):
        out += [render_doc(ast.get_docstring(init)), ""]
    for member in node.body:
        if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if member.name.startswith("_"):
                continue
            tags = decorators(member)
            if "setter" in tags or "deleter" in tags or "overload" in tags:
                continue
            out += render_function(member, rel, 3, f"{qual}.{member.name}", "method")
    return out


def public_objects(page: Page, pages: dict[str, Page]) -> list[tuple[str, ast.AST, str]]:
    """(name, node, source-relative path) for the objects documented on ``page``."""
    defs = definitions(page.tree)
    names = exported(page.tree)
    rel = page.path.relative_to(SRC).as_posix()
    out = []
    if names is None:
        for name, node in defs.items():
            if not name.startswith("_"):
                out.append((name, node, rel))
        return out
    imp = imports(page.tree, page.name)
    for name in names:
        if name in defs:
            out.append((name, defs[name], rel))
            continue
        src = imp.get(name)
        if src and src in pages and is_private(src.removeprefix("shore.")):
            other = pages[src]
            node = definitions(other.tree).get(name)
            if node is not None:
                out.append((name, node, other.path.relative_to(SRC).as_posix()))
    return out


def render_page(page: Page, pages: dict[str, Page]) -> str:
    out = ["---", f"title: {page.name}", "outline: [2, 3]", "---", "", f"# `{page.name}`", ""]
    doc = render_doc(ast.get_docstring(page.tree))
    if doc:
        out += [doc, ""]
    children = sorted(
        n
        for n in pages
        if n.rsplit(".", 1)[0] == page.name
        and n != page.name
        and not is_private(n.removeprefix("shore"))
    )
    if children:
        out += ["**Modules**", ""]
        out += [f"- [`{c}`](./{c}) — {summary(pages[c])}" for c in children]
        out.append("")
    objs = public_objects(page, pages)
    names = exported(page.tree)
    if names is not None:
        imp = imports(page.tree, page.name)
        reexports = [n for n in names if n not in {o[0] for o in objs}]
        if reexports:
            out += ["**Re-exported**", ""]
            for n in reexports:
                src = imp.get(n, "")
                link = f"./{src}#{anchor(n)}" if src in pages else None
                out.append(f"- [`{n}`]({link})" if link else f"- `{n}`")
            out.append("")
    if objs:
        for name, node, rel in objs:
            if isinstance(node, ast.ClassDef):
                out += render_class(node, rel, name)
            else:
                out += render_function(node, rel, 2, name, "function")
    page.objects = [o[0] for o in objs]
    text = "\n".join(line.rstrip() for line in "\n".join(out).splitlines())
    return re.sub(r"\n{3,}", "\n\n", text).rstrip() + "\n"


def summary(page: Page) -> str:
    doc = ast.get_docstring(page.tree) or ""
    first = doc.strip().split("\n\n")[0].replace("\n", " ")
    first = re.sub(r"\s+", " ", first)
    first = re.sub(r"\(ADR [^)]*\)", "", first).strip()
    first = re.split(r"(?<=[.:;])\s", first)[0].rstrip(".:;")
    return inline(first) if first else ""


def collect() -> dict[str, Page]:
    pages = {}
    for path in sorted(SRC.rglob("*.py")):
        name = module_name(path)
        if name in SKIP:
            continue
        pages[name] = Page(name, path, ast.parse(path.read_text(), filename=str(path)))
    return pages


def group_of(name: str) -> str:
    best, length = "Package", 0
    for title, prefixes in GROUPS:
        for p in prefixes:
            if (name == p or name.startswith(p + ".")) and len(p) > length:
                best, length = title, len(p)
    return best


def generate(out: Path) -> list[Path]:
    pages = collect()
    public = {n: p for n, p in pages.items() if not is_private(n.removeprefix("shore"))}
    out.mkdir(parents=True, exist_ok=True)
    written = []
    for name, page in public.items():
        target = out / f"{name}.md"
        target.write_text(render_page(page, pages))
        written.append(target)
    groups: dict[str, list[str]] = {t: [] for t, _ in GROUPS}
    for name in public:
        groups[group_of(name)].append(name)
    index = [
        "---",
        "title: Python API",
        "---",
        "",
        "# Python API",
        "",
        "Generated from the docstrings by `docs/scripts/api_docs.py`; do not edit",
        "these pages by hand. The command-line interface is documented under",
        "[CLI commands](/reference/), the project journal under",
        "[Project scripts](/reference/run).",
        "",
    ]
    sidebar = []
    for title, names in groups.items():
        if not names:
            continue
        index += [f"## {title}", "", "| Module | Summary |", "|---|---|"]
        index += [f"| [`{n}`](./{n}) | {summary(public[n])} |" for n in names]
        index.append("")
        sidebar.append(
            {
                "text": title,
                "collapsed": True,
                "items": [
                    {"text": n.removeprefix("shore.") or "shore", "link": f"/reference/api/{n}"}
                    for n in names
                ],
            }
        )
    (out / "index.md").write_text("\n".join(x.rstrip() for x in index).rstrip() + "\n")
    (out / "_sidebar.json").write_text(json.dumps(sidebar, indent=2) + "\n")
    return [*written, out / "index.md", out / "_sidebar.json"]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="exit 1 if the pages are stale")
    args = ap.parse_args(argv)
    if not args.check:
        for old in OUT.glob("*"):
            old.unlink()
        files = generate(OUT)
        print(f"wrote {len(files)} files under {OUT.relative_to(ROOT)}")
        return 0
    with tempfile.TemporaryDirectory() as tmp:
        files = generate(Path(tmp))
        fresh = {f.name: f.read_text() for f in files}
    current = {f.name: f.read_text() for f in OUT.glob("*")} if OUT.exists() else {}
    stale = sorted(n for n in fresh.keys() | current.keys() if fresh.get(n) != current.get(n))
    if stale:
        print("API reference is stale; run python docs/scripts/api_docs.py:", *stale, sep="\n  ")
        return 1
    print(f"API reference up to date ({len(fresh)} files)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
