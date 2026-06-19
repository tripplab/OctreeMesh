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
- `align_axes:kabsch` exports the first ten elements found in the current mesh, selects the first usable hexahedral seed from the current mesh, derives its seed frame from local element connectivity, fits a canonical axis-aligned cube with the Kabsch / orthogonal-Procrustes algorithm, applies the resulting det=+1 orthogonal rotation about the world origin to the exported nodes, and reports fit diagnostics including RMSD, orthogonality error, and seed-selection diagnostics
- `align_axes:kabsch:global` computes the seed Kabsch rotation first, uses that seed only to bucket every edge of every hexahedron to its nearest signed X/Y/Z axis, then refits one global orientation-only rotation from all cell edge directions. It applies the final global rotation to the full mesh, reports `mesh.align_axes.global.*` residual diagnostics including `mesh.align_axes.global.model=single_orientation`, and leaves snapping off.
- `align_axes:kabsch:global:snap` runs the global rotation refit, then snaps the rotated full mesh onto the rectilinear lattice using per-axis median edge lengths from the global edge buckets, routing axes with non-uniform cell sizes to per-plane clustering snap. Modifier order is independent, so `align_axes:kabsch:snap:global` is equivalent.
- `align_axes:kabsch:global:snap:cube` shares the snap plane bucketing, then writes every axis with one center-anchored cubic edge length `L*`. The default `L*` criterion is volume preservation: `L* = cbrt(step.x * step.y * step.z)`, so the bounding-box center and integer grid dimensions stay fixed while each axis span changes symmetrically by `L*/step.axis`. It refuses non-uniform global inputs and reports `mesh.align_axes.snap.cube.*` diagnostics including pre-cube steps, `L_star`, anisotropy, unchanged grid dimensions, held centers, span ratios, volume before/after, and maximum displacement.
- `align_axes:kabsch:snap` runs the seed-only Kabsch rotation, then snaps the rotated exported nodes onto the uniform rectilinear lattice defined by the Kabsch edge length `mesh.align_axes.L`; it reports `mesh.align_axes.snap.*` diagnostics and rejects non-separable or collapsing inputs instead of silently changing non-grid geometry
- `octet:<spec>` where `<spec>` is one of `+x+y+z`, `+x+y-z`, `+x-y+z`, `+x-y-z`, `-x+y+z`, `-x+y-z`, `-x-y+z`, `-x-y-z`
- `cylinder:<radius_ang>` keeps only content in a cylinder around the +Z axis through origin; node rule is `x^2 + y^2 <= radius^2` and `z > 0`, and elements are kept only when all element nodes are kept (strict containment)
- `--mesh_stats` appends a reporting operation that prints deterministic mesh summary keys (`mesh.stats.*`) without mutating geometry (fixed 6-decimal formatting for floating-point fields)
- optional operation form remains available for advanced pipelines: `--op mesh_stats` and `--op mesh_stats:format=text`
- `--perf_stats` prints pipeline stage timings and node/element counts
- Long-running CLI phases print human-oriented start messages to stderr; reading/parsing also prints coarse `meshpp read: progress=<percent>%` updates to stderr so machine-readable stdout reports remain stable.

Operations are applied in the exact order provided. This matters when combining `cylinder` with `align_axes:kabsch:*` because each operation mutates the mesh seen by the next operation. Let `N` be the number of elements in the input mesh at the point `align_axes:kabsch:*` runs, `C_before` be the number of elements kept by `cylinder:35` on the original/pre-alignment coordinates, and `C_after` be the number of elements kept by `cylinder:35` after global Kabsch rotation/snap. Then:

| Pipeline | Final exported mesh elements | Notes |
| --- | ---: | --- |
| `--op cylinder:35 --op align_axes:kabsch:global:snap` | `C_before` | Cylinder filters first; global Kabsch/snap runs only on the cylinder-kept mesh. Alignment seeds from a usable hexahedron in that current filtered mesh and fails only if no current element can provide a valid hexahedral seed. |
| `--op align_axes:kabsch:global:snap --op cylinder:35` | `C_after` | Global Kabsch/snap runs on the full mesh first; cylinder then filters the rotated/snapped mesh. The align stage may report `N` exported before the later cylinder operation reduces the final output. |
| `--op align_axes:kabsch:global:snap` | `N` | Global mode does not apply the legacy first-ten-element export filter. |
| `--op align_axes:kabsch:snap` | `min(N, 10)` | Seed-only Kabsch/snap keeps the legacy first-ten-element export behavior. |

For `rotate`, the syntax is `rotate:<ux>,<uy>,<uz>,<degrees>`. The axis vector must be finite and non-zero, is normalized internally, and positive angles follow the right-hand rule.

For `align_axes:pca`, PCA axes must be non-degenerate. For `align_axes:simple`, element order follows the input file; the first selected element must provide three non-zero orthonormal edge neighbors that form a right-handed local frame. For `align_axes:kabsch`, the operation scans the current mesh for the first usable hexahedral seed, resolves that element's eight local corners, and derives the provisional classification frame from local element connectivity rather than fixed global node ids. The dominant principal direction aligns with global Z, the second with global X, and the third with global Y. The optional `:snap` and `:global` modifiers are supported only with `kabsch`. Pipeline order is fixed regardless of modifier order: seed Kabsch (`R0`, `L0`) from the current mesh, optional global edge-direction refit over the current mesh, optional seed-only first-ten export filtering, rotation about the world origin, optional snap, then export. Without `:global`, `align_axes:kabsch` keeps the legacy first-ten-element export behavior and reports `mesh.align_axes.global=off`; seed discovery still happens before that export filter so filtered meshes can align using any usable current element. With `:global`, the full current hex connectivity is required and the final matrix diagnostics are overwritten with the global rotation; the command reports `mesh.align_axes.global.model=single_orientation`, iteration count, edge count, final bucket changes, converged-frame dominance, residual mean/p95/max in degrees, per-axis median edge lengths, and per-axis nonuniform flags. The global refit fails instead of emitting a bad rotation when all edges are degenerate or the converged frame falls below the 0.80 dominance guard, and emits a quality warning if the maximum residual exceeds 1 degree. The global refit uses all cell edges so R is the noise-averaged orientation over all cells rather than a single noisy cell's, yielding a lower, seed-independent residual. The optional `:snap` stage is post-rotation: seed-only snap uses the single Kabsch `L` on every axis, while global snap feeds each axis's median `L.x`, `L.y`, and `L.z` into the snap logic. If a global axis is non-uniform, that axis emits a `W_QUALITY` warning and uses per-plane clustering snap instead of uniform phase+k*L snap. Adding `:cube` makes the final snap write isotropic: snap still assigns the same per-node integer plane index per axis, but coordinates are rewritten as `center.axis + (k - k_mid) * L*` with one `L*` for all axes. This keeps the object stationary by holding the bounding-box center per axis exactly and keeps grid cell counts unchanged; because those requirements are over-constrained with cubic cells, per-axis spans become `n.axis * L*` and the log reports the resulting span ratios and maximum displacement. `:cube` refuses non-uniform global inputs with `E_USAGE: align_axes:snap:cube requires uniform-scale input`, warns if forcing isotropy distorts axes by more than 2%, and asserts that every hex still spans exactly two planes on each axis. When snap is not requested, the operation prints `mesh.align_axes.snap=off` and leaves the rotation-only behavior unchanged.

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
