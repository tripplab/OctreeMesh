#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
BIN="$ROOT_DIR/MESHPP/meshpp_roundtrip"
FIX="$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh"
TMP_DIR="${TMPDIR:-/tmp}/meshpp_fixture_runs"
mkdir -p "$TMP_DIR"

"$BIN" "$FIX" "$TMP_DIR/out.post.msh" --validate >/tmp/meshpp_stdout.txt
"$BIN" "$TMP_DIR/out.post.msh" "$TMP_DIR/out2.post.msh" --validate >/tmp/meshpp_stdout2.txt

rg -q "nodes: 8" /tmp/meshpp_stdout.txt
rg -q "elements: 1" /tmp/meshpp_stdout.txt

echo "meshpp roundtrip checks: PASS"
