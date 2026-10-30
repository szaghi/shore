# `shore grd-merge`

Concatenate multiple `.grd` files into one for an Xall multi-mesh run.

```bash
shore grd-merge INPUT1.grd INPUT2.grd ... -o OUTPUT.grd
```

| Argument / option | Description |
|---|---|
| `INPUTS...` | The `.grd` files, in component order. |
| `-o` / `--output PATH` | Required. The merged `.grd`. |

## Motivation

When an Xall job combines several independent Chimera components — for
example a near-body wall mesh plus a far-field background — each
component is generated as its own `.grd` + `cc.par` pair. Running them
through Xall as one job needs both pairs concatenated into a single
`assembly.grd` and a single `assembly.cc.par` whose patch indices are
shifted to refer to the merged block list.

`shore grd-merge` does the file half:

```bash
shore grd-merge background.grd wall.grd -o assembly.grd
```

The merged file lists every block from `inputs[0]` first, then
`inputs[1]`, etc. The header block count is rewritten to the sum;
shape and coordinate records are passed through byte-for-byte (no
node values are touched). Order is the caller's choice — it matters
because Xall's overset solver uses block ordering to break level /
priority ties — but it **must** match the order used for the
companion [`shore cc-par-merge`](./cc-par-merge).

## Coordinated invocation

`grd-merge` and [`cc-par-merge`](./cc-par-merge) are a pair. The
`cc.par` patch `block` and box `associated_block` fields are 1-based
indices into the `.grd` block list, so the merged `cc.par` only makes
sense if its block ordering matches the merged `.grd`. Always invoke
them with the same input ordering:

```bash
shore grd-merge    background.grd    wall.grd    -o assembly.grd
shore cc-par-merge background.cc.par wall.cc.par -o assembly.cc.par \
                   --grd assembly.grd
```

The example at `examples/05-chimera-assembly/chimera_assembly.py`
runs both halves in sequence; the bash sibling
`examples/05-chimera-assembly/chimera_assembly.sh` does the same
through `shore` CLI calls only.  See the
[cc.par + boxes pipeline](/tutorials/xall-bundle) walkthrough for the
end-to-end flow.

## Python API

```python
from shore.io.grd import merge_grd

merge_grd(
    ["background.grd", "wall.grd"],
    "assembly.grd",
)
```

Raises `ValueError` for fewer than 2 inputs or a malformed `.grd` file.

## See also

- [`shore grd`](./grd) — generate one `.grd` per Chimera component.
- [`shore cc-par-merge`](./cc-par-merge) — companion `cc.par` merger;
  invoke as a pair.
- [`shore cc-par`](./cc-par) — generate one `cc.par` per component.
