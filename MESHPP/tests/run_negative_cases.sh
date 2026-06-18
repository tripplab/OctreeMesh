#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
APPLY_BIN="$ROOT_DIR/MESHPP/meshpp_apply"
ROUNDTRIP_BIN="$ROOT_DIR/MESHPP/meshpp_roundtrip"
FIX="$ROOT_DIR/tests/fixtures/meshpp/invalid"
TMP_DIR="${TMPDIR:-/tmp}/meshpp_negative"
mkdir -p "$TMP_DIR"

expect_code() {
  local expected="$1"
  shift
  set +e
  "$@" >/tmp/meshpp_neg_stdout.txt 2>/tmp/meshpp_neg_stderr.txt
  local got=$?
  set -e
  if [[ "$got" -ne "$expected" ]]; then
    echo "Expected exit $expected, got $got for: $*" >&2
    cat /tmp/meshpp_neg_stderr.txt >&2 || true
    exit 1
  fi
}

expect_code 3 "$ROUNDTRIP_BIN" "$FIX/missing_mesh_header.post.msh" "$TMP_DIR/out.post.msh"
expect_code 3 "$ROUNDTRIP_BIN" "$FIX/duplicate_node.post.msh" "$TMP_DIR/out.post.msh"
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op nope
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op translate:1,2

expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op rotate:0,0,1
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op rotate:z,90
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op rotate:0,0,1,90,extra
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op rotate:0,0,0,90
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op rotate:nan,0,1,90
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op rotate:0,0,1,inf

expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op align_axes
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op align_axes:foo
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op align_axes:pca,extra
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op align_axes:pca
KABSCH_BAD="$TMP_DIR/kabsch_bad_topology.post.msh"
python3 - <<'PY' >"$KABSCH_BAD"
print('MESH "bad" dimension 3 ElemType Hexahedra Nnode 8')
print()
print('Coordinates')
coords = [(0,0,0), (1,0,0), (1,0,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1)]
for i, xyz in enumerate(coords, 1):
    print(i, *xyz)
print('End Coordinates')
print('Elements')
print(1, *range(1, 9))
print('End Elements')
PY
expect_code 2 "$APPLY_BIN" --in "$KABSCH_BAD" --out "$TMP_DIR/out_kabsch_bad.post.msh" --op align_axes:kabsch


expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op mesh_stats:format=xml
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op mesh_stats:foo=bar
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op mesh_stats:format=json



expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op octet:x+y+z
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op octet:+x+y
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op octet:+x+q+z
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op cylinder:0
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op cylinder:-3
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --op cylinder:abc
expect_code 2 "$APPLY_BIN" --in "$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh" --out "$TMP_DIR/out2.post.msh" --cylinder foo

echo "meshpp negative checks: PASS"
