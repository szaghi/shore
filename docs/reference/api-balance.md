# shore.balance

Python API for the **multi-way equal-partition** split planner. The CLI wrapper
[`shore balance`](./balance) is a thin adapter around `plan_balance`.

The algorithm, the multigrid floor, and the design rationale are documented at
[`shore balance`](./balance) and
[Design decisions](/guide/design-decisions#multi-way-equal-partition-not-bisection-issue-6).

## Public exports

| Name | Role |
|---|---|
| `plan_balance` | The planner. Reads a `.grd` (optionally `.adjacency.json` + wired blocks) and returns a `BalancePlan`. |
| `BalancePlan` | Dataclass holding the planned splits + before/after summaries. `to_toml()` renders a `shore split` config; `to_bsplit_par()` renders a `bsplit.par`. |

## `plan_balance`

```python
def plan_balance(
    grd_path: str | Path,
    *,
    np: int,
    blocks_per_rank: float = 1.0,
    ngr: int = 1,
    min_coarse_cells: int = 8,
    tolerance: float = 0.05,
    adjacency_path: str | Path | None = None,
    blocks: list | None = None,
) -> BalancePlan
```

Partition every block into roughly equal, roughly cubic chunks toward
`np * blocks_per_rank`, keeping every chunk axis above the multigrid floor.

| Parameter | Description |
|---|---|
| `grd_path` | Path to the input `.grd` (block ordering + cell counts). |
| `np` | Total MPI rank count (>= 1). |
| `blocks_per_rank` | Target average blocks per rank (default `1.0`). Higher → finer granularity / better balance, more blocks. |
| `ngr` | Multigrid levels the mesh must support (default `1`). Chunk axes are kept divisible by `2^(ngr-1)`. |
| `min_coarse_cells` | Minimum cells per axis on the coarsest (level-`ngr`) grid (default `8`). Chunk axes are kept `>= min_coarse_cells * 2^(ngr-1)` fine cells. |
| `tolerance` | Target fractional imbalance `(max - min) / mean` for the `converged` flag (default `0.05`). Never causes a failure. |
| `adjacency_path` | Optional `.adjacency.json` sidecar — used for labels / block-count validation. |
| `blocks` | Optional list of wired `HexBlock`s. When given, the planner floor-checks against **seam-propagated** cuts (via [`shore.split.propagate_cuts`](./api-split)) so the plan matches the split output. Without it, the floor is checked per block in isolation, and the plan uses synthetic `block_<n>` labels (required for `to_bsplit_par`). |

**Returns** a `BalancePlan`. Best-effort: when the floor prevents reaching
`tolerance`, it warns (`UserWarning`) and returns the best plan with
`converged=False`.

**Raises**

- `ParameterError` — `np < 1`, `ngr < 1`, `min_coarse_cells < 1`,
  `tolerance <= 0`, or `blocks_per_rank <= 0`.
- `FileNotFoundError`, `ValueError` — forwarded from `read_grd_metadata` and the adjacency reader.

## `BalancePlan`

```python
@dataclass
class BalancePlan:
    splits:           list[tuple[str, str, int]]  # (label, axis, vertex_index)
    iterations:       int
    converged:        bool
    initial_summary:  BalanceSummary
    final_summary:    BalanceSummary
    tolerance:        float
```

| Attribute | Description |
|---|---|
| `splits` | Cuts in decision order: `(original_block_label, axis, original_block_vertex_index)`. Vertex indices are relative to the original block, not any intermediate chunk — this is what `shore split`'s TOML schema expects. |
| `iterations` | Number of split iterations the planner actually ran. |
| `converged` | `True` iff the final imbalance is within `tolerance`. |
| `initial_summary` | [`BalanceSummary`](./api-proc-input#balancesummary) before any splits. |
| `final_summary` | `BalanceSummary` after the last split. |
| `tolerance` | The tolerance the planner targeted. |

### `BalancePlan.to_toml`

```python
def to_toml(self) -> str
```

Render the plan as a `shore split` config TOML (`version = 1`,
`kind = "split"`, one `[[splits]]` per `(label, axis)` pair with a
sorted `at = [...]` list).  Same schema as the example in the
[`shore split`](./split) reference.

### `BalancePlan.to_bsplit_par`

```python
def to_bsplit_par(
    self, *,
    ngr: int = 1,
    sol_root: str = "sol",
    grd_root: str = "cc",
    icc_root: str = "cc",
    varn: int = 0,
    varo: int = 0,
    debug: int = 0,
) -> str
```

Render the plan as a `bsplit.par` for the external `bsplit` tool, which splits
the **framed overset outputs** *after* an overset run (the post-overset
analogue of `to_toml` + `shore split`). See
[After-overset splitting](/guide/after-overset-bsplit).

Requires a **grd-only** plan (no `adjacency_path`/`blocks`): its synthetic
`block_<n>` labels carry the 1-based grd block order, which is exactly bsplit's
block number. Raises `ParameterError` if any split label is not `block_<n>`.

## Example

```python
from shore.balance import plan_balance

plan = plan_balance(
    "wall.grd",
    np=16,
    tolerance=0.05,
    adjacency_path="wall.adjacency.json",
)

print(f"converged: {plan.converged}")
print(f"initial: {plan.initial_summary.imbalance * 100:.1f}%")
print(f"final:   {plan.final_summary.imbalance * 100:.1f}%")
print(f"splits:  {len(plan.splits)} cuts across "
      f"{len({s[0] for s in plan.splits})} blocks")

with open("wall.splits.toml", "w") as f:
    f.write(plan.to_toml())
```

## See also

- [`shore balance` CLI](./balance) — wrapper for `plan_balance`.
- [`shore.split`](./api-split) — applies the plan.
- [`shore.io.proc_input`](./api-proc-input) — final rank assignment.
- [`shore.io.grd.read_grd_metadata`](./api-grd) — the header-only
  reader the planner uses.
