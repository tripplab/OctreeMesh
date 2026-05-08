# Octet operation (`octet`)

`octet` filters mesh content to a single Cartesian octant around the global origin `(0,0,0)`.

## Syntax

- `--op octet:+x+y+z`
- `--op octet:+x+y-z`
- `--op octet:+x-y+z`
- `--op octet:+x-y-z`
- `--op octet:-x+y+z`
- `--op octet:-x+y-z`
- `--op octet:-x-y+z`
- `--op octet:-x-y-z`

Accepted grammar: `^[+-]x[+-]y[+-]z$`.

## Selection rules (v1)

- `+x` means `x >= 0`, `-x` means `x < 0`.
- `+y` means `y >= 0`, `-y` means `y < 0`.
- `+z` means `z >= 0`, `-z` means `z < 0`.

Boundary nodes on `x=0`, `y=0`, `z=0` are always assigned to the `+` side.

## Element policy (v1)

Elements are retained only when **all** referenced nodes are inside the octet.
Cross-boundary elements are dropped (no geometric clipping in v1).

## Output behavior

- Node IDs and element IDs are preserved for retained entities.
- Unreferenced nodes are removed from output.
- Empty selections succeed and print `mesh.octet.warning=selection produced no elements`.

## Ordering in operation pipelines

Operations are evaluated in command order. For example:

- `--op octet:+x+y+z --op scale:2` filters first, then scales.
- `--op scale:2 --op octet:+x+y+z` scales first, then filters.
