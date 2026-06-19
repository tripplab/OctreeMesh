# meshpp

`meshpp` is a standalone `.post.msh` processing tool.

## Build

```bash
cd MESHPP
make
```

## CLI tools

### 1) Roundtrip validation

```bash
./meshpp_roundtrip <input.post.msh> <output.post.msh> [--validate] [--perf_stats]
```

- `--validate`: prints node/element counts
- `--perf_stats`: prints stage timings (`read`, `validate`, `write`)

### 2) Apply operation pipeline

```bash
./meshpp_apply --in <input.post.msh> --out <output.post.msh> --op <spec> [--op <spec> ...] [--mesh_stats] [--perf_stats]
```

Supported operations:
- `scale:<factor>`
- `translate:<dx>,<dy>,<dz>`
- `rotate:<ux>,<uy>,<uz>,<degrees>` rotates all node coordinates around the world origin using an arbitrary axis vector; the axis is normalized internally and positive angles use the right-hand rule
- `align_axes:pca` computes PCA principal directions from all nodes, rotates the mesh around the world origin so the second/third/first principal directions align with X/Y/Z respectively, and prints the PCA eigenvalues, PCA eigenvectors, and rotation matrix used
- `align_axes:simple` exports the first ten elements found in the input mesh, uses the first selected element first node and its three nearest element neighbors to build a right-handed local cube frame, rotates the exported nodes so those local axes align with world X/Y/Z, and reports the frame vectors, checks, and rotation matrix
- `align_axes:kabsch` exports the first ten elements found in the input mesh, uses all eight corners of the first element to fit a canonical axis-aligned cube with the Kabsch / orthogonal-Procrustes algorithm, applies the resulting det=+1 orthogonal rotation about the world origin to the exported nodes, and reports fit diagnostics including RMSD and orthogonality error
- `align_axes:kabsch:global` computes the seed Kabsch rotation first, uses that seed only to bucket every edge of every hexahedron to its nearest signed X/Y/Z axis, then refits one global orientation-only rotation from all cell edge directions. It applies the final global rotation to the full mesh, reports `mesh.align_axes.global.*` residual diagnostics, and leaves snapping off.
- `align_axes:kabsch:global:snap` runs the global rotation refit, then snaps the rotated full mesh onto the rectilinear lattice using per-axis median edge lengths from the global edge buckets. Modifier order is independent, so `align_axes:kabsch:snap:global` is equivalent.
- `align_axes:kabsch:snap` runs the seed-only Kabsch rotation, then snaps the rotated exported nodes onto the uniform rectilinear lattice defined by the Kabsch edge length `mesh.align_axes.L`; it reports `mesh.align_axes.snap.*` diagnostics and rejects non-separable or collapsing inputs instead of silently changing non-grid geometry
- `octet:<spec>` where `<spec>` is one of `+x+y+z`, `+x+y-z`, `+x-y+z`, `+x-y-z`, `-x+y+z`, `-x+y-z`, `-x-y+z`, `-x-y-z`
- `cylinder:<radius_ang>` keeps only content in a cylinder around the +Z axis through origin; node rule is `x^2 + y^2 <= radius^2` and `z > 0`, and elements are kept only when all element nodes are kept (strict containment)
- `--mesh_stats` appends a reporting operation that prints deterministic mesh summary keys (`mesh.stats.*`) without mutating geometry (fixed 6-decimal formatting for floating-point fields)
- optional operation form remains available for advanced pipelines: `--op mesh_stats` and `--op mesh_stats:format=text`
- `--perf_stats` prints pipeline stage timings and node/element counts

Operations are applied in the exact order provided. This matters when combining `cylinder` with `align_axes:kabsch:*` because each operation mutates the mesh seen by the next operation. Let `N` be the number of elements in the input mesh at the point `align_axes:kabsch:*` runs, `C_before` be the number of elements kept by `cylinder:35` on the original/pre-alignment coordinates, and `C_after` be the number of elements kept by `cylinder:35` after global Kabsch rotation/snap. Then:

| Pipeline | Final exported mesh elements | Notes |
| --- | ---: | --- |
| `--op cylinder:35 --op align_axes:kabsch:global:snap` | `C_before` | Cylinder filters first; global Kabsch/snap runs only on the cylinder-kept mesh. If the cylinder leaves no elements, or removes the seed element/nodes required by Kabsch, the alignment fails instead of exporting a mesh. |
| `--op align_axes:kabsch:global:snap --op cylinder:35` | `C_after` | Global Kabsch/snap runs on the full mesh first; cylinder then filters the rotated/snapped mesh. The align stage may report `N` exported before the later cylinder operation reduces the final output. |
| `--op align_axes:kabsch:global:snap` | `N` | Global mode does not apply the legacy first-ten-element export filter. |
| `--op align_axes:kabsch:snap` | `min(N, 10)` | Seed-only Kabsch/snap keeps the legacy first-ten-element export behavior. |

For `rotate`, the syntax is `rotate:<ux>,<uy>,<uz>,<degrees>`. The axis vector must be finite and non-zero, is normalized internally, and positive angles follow the right-hand rule.

For `align_axes:pca`, PCA axes must be non-degenerate. For `align_axes:simple`, element order follows the input file; the first selected element must provide three non-zero orthonormal edge neighbors that form a right-handed local frame. For `align_axes:kabsch`, the first element must provide the eight distinct cube corners and node ids 1, 5, 2, and 4 seed the provisional classification frame. The dominant principal direction aligns with global Z, the second with global X, and the third with global Y. The optional `:snap` and `:global` modifiers are supported only with `kabsch`. Pipeline order is fixed regardless of modifier order: seed Kabsch (`R0`, `L0`), optional global edge-direction refit, rotation about the world origin, optional snap, then export. Without `:global`, `align_axes:kabsch` keeps the legacy first-ten-element export behavior and reports `mesh.align_axes.global=off`. With `:global`, the full hex connectivity is required and the final matrix diagnostics are overwritten with the global rotation; the command reports iteration count, edge count, final bucket changes, dominance, residual mean/p95/max in degrees, and per-axis median edge lengths. The global refit fails instead of emitting a bad rotation when all edges are degenerate or edge bucketing falls below the 0.80 dominance guard, and emits a quality warning if the maximum residual exceeds 1 degree. The optional `:snap` stage is post-rotation: seed-only snap uses the single Kabsch `L` on every axis, while global snap feeds the per-axis median `L.x`, `L.y`, and `L.z` into the unchanged axis-independent snap logic. When snap is not requested, the operation prints `mesh.align_axes.snap=off` and leaves the rotation-only behavior unchanged.

For `octet`, half-space boundaries are deterministic: `+` means `>= 0`, `-` means `< 0`.
Elements are kept only when all element nodes satisfy the selected octet.

For `cylinder`, boundary inclusion is deterministic (`<= radius`) and the cylinder is infinite along Z with the additional `z > 0` constraint for node inclusion.

## Exit codes

- `0`: success
- `2`: usage/configuration error
- `3`: parse/validation error
- `4`: unsupported input
- `5`: topology/consistency error
- `6`: I/O error

## Test and sanitizer checks

```bash
./tests/ci.sh
make sanitize
./meshpp_roundtrip_asan ../tests/fixtures/meshpp/valid/single_hex.post.msh /tmp/out.post.msh --validate
./meshpp_apply_asan --in ../tests/fixtures/meshpp/valid/single_hex.post.msh --out /tmp/out2.post.msh --op scale:1.1 --op translate:0,0,0
```
