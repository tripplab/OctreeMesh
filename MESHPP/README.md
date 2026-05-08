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
- `--mesh_stats` appends a reporting operation that prints deterministic mesh summary keys (`mesh.stats.*`) without mutating geometry (fixed 6-decimal formatting for floating-point fields)
- optional operation form remains available for advanced pipelines: `--op mesh_stats` and `--op mesh_stats:format=text`
- `--perf_stats` prints pipeline stage timings and node/element counts

Operations are applied in the exact order provided.

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
