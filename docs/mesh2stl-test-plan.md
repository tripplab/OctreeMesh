# mesh2stl Test Plan (Ticket 2)

This plan defines fixture-driven validation for `mesh2stl` per `docs/mesh2stl-spec.md`.

## 1) Goals

- Freeze expected behavior before implementation.
- Provide deterministic, reviewable fixture meshes.
- Define acceptance checks for parsing, topology extraction, and output accounting.

## 2) Fixture Matrix

| Fixture | Class | Purpose | Expected Result |
|---|---|---|---|
| `valid/single_hex_cube.gid` | valid | Minimal 1-hexahedron mesh | Parse success; 6 boundary quads; 12 triangles in `--surface-only` |
| `valid/two_hex_shared_face.gid` | valid | Interior face cancellation | Parse success; 10 boundary quads; 20 triangles in `--surface-only`; 12 quads/24 triangles in `--all-faces` |
| `valid/two_components.gid` | valid | Disconnected components | Parse success; 12 boundary quads; 24 triangles in `--surface-only` |
| `invalid/missing_node_reference.gid` | invalid | Element references missing node id | Parse/topology failure with input validation error |
| `invalid/non_hexa_element.gid` | invalid | Unsupported element type | Fail with unsupported feature error (v1) |
| `invalid/non_manifold_face.gid` | invalid | >2 elements share one face | Fail with topology validation error in `--surface-only` |
| `invalid/malformed_mesh_section.gid` | invalid | Broken GiD section syntax | Parse/format error |

## 3) Expected Behavioral Checks

### 3.1 Parse and validation
- Valid fixtures parse into node + hexa element collections.
- Invalid fixtures fail deterministically with the corresponding exit-class from the spec.

### 3.2 Topology accounting (`--surface-only`)
- Boundary faces are faces with exactly one owner.
- Interior faces (2 owners) are excluded.
- Non-manifold faces (>2 owners) are rejected as topology errors.

### 3.3 Determinism
For each valid fixture and fixed CLI options:
- face extraction order is deterministic,
- triangle count is deterministic,
- STL output ordering is deterministic (byte-for-byte or canonical equivalent per implementation policy).

## 4) Golden Assertions by Fixture

- `single_hex_cube.gid`: 1 element => 6 boundary quads => 12 triangles.
- `two_hex_shared_face.gid`: 2 elements sharing 1 face => 12 total element faces, 2 interior-face owners removed => 10 boundary quads => 20 triangles.
- `two_components.gid`: 2 disconnected elements => 12 boundary quads => 24 triangles.

## 5) Error-Class Assertions

- Missing node reference -> input parse/format or validation error (spec exit code class for invalid input).
- Non-hexa element -> unsupported-feature error class.
- Non-manifold face in `--surface-only` -> topology validation error class.
- Malformed section -> input parse/format error class.

## 6) Review Checklist (Ticket 2)

- [ ] Every fixture has a documented expected outcome.
- [ ] At least one malformed input exists for each major failure class:
  - parse/format
  - unsupported element type
  - topology validation
- [ ] Determinism assertions are stated for valid fixtures.
- [ ] Counts are manually derivable from fixture topology.

## 7) CI / Automation Hook (Ticket 6)

For repeatable end-to-end checks, run:

```bash
cd MESH2STL
./tests/ci.sh
```

This performs:
- clean rebuild,
- valid fixture conversions,
- invalid fixture exit-code assertions,
- non-manifold mode split checks (`--surface-only` fails, `--all-faces` passes).
