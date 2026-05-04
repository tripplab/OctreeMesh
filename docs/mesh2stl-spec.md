# mesh2stl v1 Specification (Ticket 1)

## 1) Purpose and Scope

### Goal
`mesh2stl` is a standalone command-line tool that reads a GiD hexahedral mesh file and writes an STL surface mesh for downstream visualization/geometry workflows.

### Non-goals (v1)
- No FEM results transfer (that remains in `mesh2pdb`).
- No surface smoothing, remeshing, or simplification.
- No volumetric export formats.
- No support for non-hexahedral element types in v1.

## 2) Input and Output

### Input
- **Required:** GiD mesh file containing:
  - node coordinates
  - hexahedral element connectivity

### Output
- **Required:** STL file path (`.stl`)
- **Default format:** ASCII STL
- **Optional format:** Binary STL is out of scope for v1.

## 3) CLI Contract

## Command
```bash
mesh2stl <input_mesh.gid> <output_mesh.stl> [options]
```

### Options (v1)
- `--surface-only` (default): export only boundary faces.
- `--all-faces`: export all element faces (debug mode).
- `--validate`: run mesh validation and print report before conversion.
- `--units <name>`: annotate units in logs only (STL remains unitless).
- `--help`: print usage and exit.

### Mutually exclusive behavior
- `--surface-only` and `--all-faces` are mutually exclusive.
- If neither is provided, `--surface-only` is assumed.

### CLI examples
```bash
# Default boundary extraction
mesh2stl capside.mesh capside.stl

# Explicit boundary extraction with validation
mesh2stl capside.mesh capside.stl --surface-only --validate

# Debug: export all element faces
mesh2stl capside.mesh capside_allfaces.stl --all-faces
```

## 4) Geometry and Topology Semantics

### Supported elements
- Only 8-node hexahedra are supported in v1.

### Coordinate policy
- Preserve the coordinates from the input mesh (no normalization or rotation in v1).

### Surface extraction mode
- In `--surface-only`, a face is boundary if it is owned by exactly one element.
- Faces owned by exactly two elements are interior and excluded.
- Faces with ownership >2 are non-manifold and treated as topology errors (see errors section).

### Face triangulation
- Each quadrilateral boundary face is split into two triangles using a fixed local rule to ensure deterministic output.

### Normals and winding
- Triangle normals are computed from vertex cross products.
- Output winding must be deterministic and consistent across runs for identical input.
- Perfect outward orientation is best-effort in v1; determinism is required.

### Degenerate geometry
- Zero-area triangles are skipped and counted.
- If the skipped count exceeds 0, tool still succeeds unless `--validate` is set to strict mode in future versions (strict mode not part of v1).

## 5) Error Handling and Exit Codes

- `0`: success.
- `2`: command-line usage error (missing args, invalid option combos).
- `3`: input parse/format error.
- `4`: unsupported feature in v1 (e.g., non-hexa element types).
- `5`: topology validation error (e.g., non-manifold face in `--surface-only`).
- `6`: filesystem I/O error (cannot read input or write output).

Error messages must:
- include a stable error code label,
- include file path when relevant,
- include line/record context when available.

## 6) Logging and Progress

v1 runtime logs should include:
- input file path, output file path
- node count, element count
- mode (`surface-only` or `all-faces`)
- generated quad-face count and final triangle count
- skipped degenerate triangle count

## 7) Determinism Requirements

For identical input and options:
- face extraction and triangulation order must be deterministic,
- output STL must be deterministic (same facet order and numeric formatting policy).

## 8) Acceptance Checklist for Ticket 1

- [ ] Product scope and non-goals documented.
- [ ] CLI contract documented with defaults and mutually exclusive options.
- [ ] At least 3 concrete CLI examples included.
- [ ] Surface extraction semantics specified.
- [ ] Coordinate, normal, and degeneracy policies specified.
- [ ] Exit code table specified.
- [ ] Determinism requirement specified.

## 9) Explicit Out-of-Scope Items for Later Tickets

- Binary STL output (`--binary`).
- JSON validation report.
- Non-hexa element support (tet/prism/pyramid).
- Guaranteed outward normal correction for arbitrary malformed meshes.
- Mesh optimization/simplification.
