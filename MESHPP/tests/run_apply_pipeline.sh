#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
BIN="$ROOT_DIR/MESHPP/meshpp_apply"
ROUNDTRIP="$ROOT_DIR/MESHPP/meshpp_roundtrip"
FIX="$ROOT_DIR/tests/fixtures/meshpp/valid/single_hex.post.msh"
TMP_DIR="${TMPDIR:-/tmp}/meshpp_pipeline_runs"
mkdir -p "$TMP_DIR"

OUT="$TMP_DIR/out_scaled_translated.post.msh"
"$BIN" --in "$FIX" --out "$OUT" --op scale:2 --op translate:1,0,-1 >/tmp/meshpp_apply_stdout.txt

# Parse output and make sure it stays valid
"$ROUNDTRIP" "$OUT" "$TMP_DIR/out2.post.msh" --validate >/tmp/meshpp_apply_validate.txt
rg -q "nodes: 8" /tmp/meshpp_apply_validate.txt
rg -q "elements: 1" /tmp/meshpp_apply_validate.txt

echo "meshpp operation pipeline checks: PASS"
