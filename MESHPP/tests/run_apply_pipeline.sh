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


STATS_OUT="$TMP_DIR/stats_stdout.txt"
"$BIN" --in "$FIX" --out "$TMP_DIR/out_stats.post.msh" --mesh_stats >"$STATS_OUT"
rg -n "^mesh\.stats\.nodes=8$" "$STATS_OUT"
rg -n "^mesh\.stats\.elements=1$" "$STATS_OUT"
rg -n "^mesh\.stats\.min\.x=0\.000000$" "$STATS_OUT"
rg -n "^mesh\.stats\.max\.z=1\.000000$" "$STATS_OUT"

TRANSFORM_STATS_OUT="$TMP_DIR/transform_stats_stdout.txt"
"$BIN" --in "$FIX" --out "$TMP_DIR/out_transform_stats.post.msh" --op scale:2 --op translate:1,0,-1 --mesh_stats >"$TRANSFORM_STATS_OUT"
rg -n "^mesh\.stats\.min\.x=1\.000000$" "$TRANSFORM_STATS_OUT"
rg -n "^mesh\.stats\.min\.z=-1\.000000$" "$TRANSFORM_STATS_OUT"
rg -n "^mesh\.stats\.max\.x=3\.000000$" "$TRANSFORM_STATS_OUT"
rg -n "^mesh\.stats\.center\.x=2\.000000$" "$TRANSFORM_STATS_OUT"
rg -n "^mesh\.stats\.bbox\.diag=3\.464102$" "$TRANSFORM_STATS_OUT"

OCTET_FIX="$ROOT_DIR/tests/fixtures/meshpp/valid/two_hex_x_split.post.msh"
OCTET_OUT="$TMP_DIR/out_octet_plusx.post.msh"
OCTET_STDOUT="$TMP_DIR/out_octet_plusx.stdout"
"$BIN" --in "$OCTET_FIX" --out "$OCTET_OUT" --op octet:+x+y+z >"$OCTET_STDOUT"
"$ROUNDTRIP" "$OCTET_OUT" "$TMP_DIR/out_octet_plusx_roundtrip.post.msh" --validate >"$TMP_DIR/out_octet_plusx.validate"
rg -q "nodes: 8" "$TMP_DIR/out_octet_plusx.validate"
rg -q "elements: 1" "$TMP_DIR/out_octet_plusx.validate"
rg -n "^mesh\.octet\.spec=\+x\+y\+z$" "$OCTET_STDOUT"
rg -n "^mesh\.octet\.nodes\.kept=8$" "$OCTET_STDOUT"
rg -n "^mesh\.octet\.elements\.kept=1$" "$OCTET_STDOUT"

OCTET_EMPTY_OUT="$TMP_DIR/out_octet_minusx_minusy_minusz.post.msh"
OCTET_EMPTY_STDOUT="$TMP_DIR/out_octet_minusx_minusy_minusz.stdout"
"$BIN" --in "$FIX" --out "$OCTET_EMPTY_OUT" --op octet:-x-y-z >"$OCTET_EMPTY_STDOUT"
rg -n "^mesh\.octet\.warning=selection produced no elements$" "$OCTET_EMPTY_STDOUT"

FIRST_LINE=$(sed -n "1p" "$STATS_OUT")
SECOND_LINE=$(sed -n "2p" "$STATS_OUT")
THIRD_LINE=$(sed -n "3p" "$STATS_OUT")
if [[ "$FIRST_LINE" != "mesh.stats.nodes=8" || "$SECOND_LINE" != "mesh.stats.elements=1" || "$THIRD_LINE" != "mesh.stats.min.x=0.000000" ]]; then
  echo "stats output order changed" >&2
  cat "$STATS_OUT" >&2
  exit 1
fi


echo "meshpp operation pipeline checks: PASS"
