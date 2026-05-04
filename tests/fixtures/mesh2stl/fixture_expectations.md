# Fixture Expectations

## Valid Fixtures

### valid/single_hex_cube.gid
- Nodes: 8
- Elements: 1 hex
- `--surface-only`: 6 quads, 12 triangles
- `--all-faces`: 6 quads, 12 triangles

### valid/two_hex_shared_face.gid
- Nodes: 12
- Elements: 2 hex
- Shared interior face: 1
- `--surface-only`: 10 quads, 20 triangles
- `--all-faces`: 12 quads, 24 triangles

### valid/two_components.gid
- Nodes: 16
- Elements: 2 hex (disconnected)
- `--surface-only`: 12 quads, 24 triangles
- `--all-faces`: 12 quads, 24 triangles

## Invalid Fixtures

### invalid/missing_node_reference.gid
- Intent: element references node id not present in coordinate list.
- Expected class: invalid input/validation error.

### invalid/non_hexa_element.gid
- Intent: include tetra element declaration.
- Expected class: unsupported feature (v1 hex-only).

### invalid/non_manifold_face.gid
- Intent: three elements share one face.
- Expected class: topology validation error in `--surface-only`.

### invalid/malformed_mesh_section.gid
- Intent: malformed `Coordinates` section terminator.
- Expected class: parse/format error.
