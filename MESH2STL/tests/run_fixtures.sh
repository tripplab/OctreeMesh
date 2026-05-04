#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
BIN="$ROOT_DIR/MESH2STL/mesh2stl"
FIX="$ROOT_DIR/tests/fixtures/mesh2stl"
TMP_DIR="${TMPDIR:-/tmp}/mesh2stl_fixture_runs"
mkdir -p "$TMP_DIR"

expect_ok() {
  local input="$1"
  shift
  "$BIN" "$input" "$TMP_DIR/out.stl" "$@" >/tmp/mesh2stl_stdout.txt 2>/tmp/mesh2stl_stderr.txt
}

expect_code() {
  local expected="$1"
  shift
  set +e
  "$BIN" "$@" >/tmp/mesh2stl_stdout.txt 2>/tmp/mesh2stl_stderr.txt
  local got=$?
  set -e
  if [[ "$got" -ne "$expected" ]]; then
    echo "Expected exit $expected but got $got for: $*" >&2
    cat /tmp/mesh2stl_stderr.txt >&2 || true
    exit 1
  fi
}

# Valid fixtures
expect_ok "$FIX/valid/single_hex_cube.gid" --validate
expect_ok "$FIX/valid/two_hex_shared_face.gid" --validate
expect_ok "$FIX/valid/two_hex_shared_face.gid" --all-faces --validate
expect_ok "$FIX/valid/two_components.gid" --validate

# Invalid fixtures
expect_code 3 "$FIX/invalid/malformed_mesh_section.gid" "$TMP_DIR/out.stl" --validate
expect_code 3 "$FIX/invalid/missing_node_reference.gid" "$TMP_DIR/out.stl" --validate
expect_code 4 "$FIX/invalid/non_hexa_element.gid" "$TMP_DIR/out.stl" --validate
expect_code 5 "$FIX/invalid/non_manifold_face.gid" "$TMP_DIR/out.stl" --validate

# Non-manifold accepted in all-faces mode
expect_ok "$FIX/invalid/non_manifold_face.gid" --all-faces --validate

echo "mesh2stl fixture checks: PASS"
