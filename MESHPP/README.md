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
- `octet:<spec>` where `<spec>` is one of `+x+y+z`, `+x+y-z`, `+x-y+z`, `+x-y-z`, `-x+y+z`, `-x+y-z`, `-x-y+z`, `-x-y-z`
- `cylinder:<radius_ang>` keeps only content in a cylinder around the +Z axis through origin; node rule is `x^2 + y^2 <= radius^2` and `z > 0`, and elements are kept only when all element nodes are kept (strict containment)
- `--mesh_stats` appends a reporting operation that prints deterministic mesh summary keys (`mesh.stats.*`) without mutating geometry (fixed 6-decimal formatting for floating-point fields)
- optional operation form remains available for advanced pipelines: `--op mesh_stats` and `--op mesh_stats:format=text`
- `--perf_stats` prints pipeline stage timings and node/element counts

Operations are applied in the exact order provided.

For `rotate`, the syntax is `rotate:<ux>,<uy>,<uz>,<degrees>`. The axis vector must be finite and non-zero, is normalized internally, and positive angles follow the right-hand rule.

For `align_axes:pca`, PCA axes must be non-degenerate. For `align_axes:simple`, element order follows the input file; the first selected element must provide three non-zero orthonormal edge neighbors that form a right-handed local frame. For `align_axes:kabsch`, the first element must provide the eight distinct cube corners and node ids 1, 5, 2, and 4 seed the provisional classification frame. The dominant principal direction aligns with global Z, the second with global X, and the third with global Y.

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
