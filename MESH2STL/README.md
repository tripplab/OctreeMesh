# mesh2stl

Standalone v1 CLI tool to convert GiD hexahedral meshes to ASCII STL.

## Build

```bash
cd MESH2STL
make
```

## Usage

```bash
./mesh2stl <input_mesh.gid> <output_mesh.stl> [options]
```

### Options
- `--surface-only` (default): export boundary faces only.
- `--all-faces`: export all element faces.
- `--validate`: print parser/topology/write summary.
- `--units <name>`: report units in logs only.
- `--help`: print usage.

## Examples

```bash
# Boundary-only export
./mesh2stl ../tests/fixtures/mesh2stl/valid/single_hex_cube.gid /tmp/single.stl

# Validation report with all faces
./mesh2stl ../tests/fixtures/mesh2stl/valid/two_hex_shared_face.gid /tmp/shared_all.stl --all-faces --validate
```

## Current limitations (v1)
- ASCII STL only.
- Hexahedra only (`ElemType Hexahedra`).
- No outward-normal correction guarantees for malformed meshes.

## Troubleshooting
- Exit `2`: CLI usage/option error.
- Exit `3`: parse or input validation error.
- Exit `4`: unsupported input (e.g., non-hexa elements).
- Exit `5`: topology validation error (e.g., non-manifold face in `--surface-only`).
- Exit `6`: I/O error opening input/output files.
