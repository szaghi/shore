# shore.balance

Python API for the iterative split planner.  The CLI wrapper
[`shore balance`](./balance) is a thin adapter around
`plan_balance`.

The algorithm and limitations are documented at
[`shore balance`](./balance).

## Public exports

| Name | Role |
|---|---|
| `plan_balance` | The planner.  Reads a `.grd` (and optionally an `.adjacency.json`) and returns a `BalancePlan`. |
| `BalancePlan` | Dataclass holding the planned splits + before/after summaries.  `to_toml()` renders a `shore split`-compatible config. |

## `plan_balance`

```python
def plan_balance(
    grd_path: str | Path,
    *,
    np: int,
    tolerance: float = 0.05,
    max_iterations: int = 100,
    adjacency_path: str | Path | None = None,
) -> BalancePlan
```

Plan the splits needed to bring the load imbalance below
`tolerance`.

| Parameter | Description |
|---|---|
| `grd_path` | Path to the input `.grd`. |
| `np` | Total MPI rank count (>= 1). |
| `tolerance` | Target fractional imbalance `(max - min) / mean`. Default `0.05` (5%). |
| `max_iterations` | Cap on planning iterations. Default `100` for the Python API; the CLI defaults to `200`. |
| `adjacency_path` | Optional path to a `.adjacency.json` sidecar. When supplied, the planner refuses to break SHARED / DIRICHLET seams without paired co-splits and prefers axes free of seam coupling. |

**Returns** a `BalancePlan`.

**Raises**

- `ParameterError` — `np < 1`, `tolerance <= 0`, or `max_iterations < 1`.
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
