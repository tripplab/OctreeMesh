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

echo "meshpp negative checks: PASS"
