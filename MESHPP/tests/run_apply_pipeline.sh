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

CYL_BOUNDARY_OUT="$TMP_DIR/out_cylinder_boundary.post.msh"
CYL_BOUNDARY_STDOUT="$TMP_DIR/out_cylinder_boundary.stdout"
"$BIN" --in "$FIX" --out "$CYL_BOUNDARY_OUT" --op translate:0,0,1 --op cylinder:1.4142135623730951 >"$CYL_BOUNDARY_STDOUT"
rg -n "^mesh\.cylinder\.radius=1\.41421" "$CYL_BOUNDARY_STDOUT"
rg -n "^mesh\.cylinder\.nodes_kept=8$" "$CYL_BOUNDARY_STDOUT"
rg -n "^mesh\.cylinder\.elements_kept=1$" "$CYL_BOUNDARY_STDOUT"

CYL_STRICT_OUT="$TMP_DIR/out_cylinder_strict.post.msh"
CYL_STRICT_STDOUT="$TMP_DIR/out_cylinder_strict.stdout"
"$BIN" --in "$FIX" --out "$CYL_STRICT_OUT" --op translate:0,0,1 --op cylinder:1 >"$CYL_STRICT_STDOUT"
rg -n "^mesh\.cylinder\.nodes_kept=0$" "$CYL_STRICT_STDOUT"
rg -n "^mesh\.cylinder\.elements_kept=0$" "$CYL_STRICT_STDOUT"

CYL_ORDER_A_STDOUT="$TMP_DIR/out_cylinder_order_a.stdout"
"$BIN" --in "$FIX" --out "$TMP_DIR/out_cylinder_order_a.post.msh" --op translate:0,0,1 --op cylinder:1.4142135623730951 >"$CYL_ORDER_A_STDOUT"
rg -n "^mesh\.cylinder\.elements_kept=1$" "$CYL_ORDER_A_STDOUT"

assert_node_xyz() {
  local mesh_file="$1"
  local node_id="$2"
  local expected_x="$3"
  local expected_y="$4"
  local expected_z="$5"
  python3 - "$mesh_file" "$node_id" "$expected_x" "$expected_y" "$expected_z" <<'PY'
import sys
mesh_file, node_id, ex, ey, ez = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
tol = 1e-6
in_coordinates = False
with open(mesh_file, encoding="utf-8") as fh:
    for line in fh:
        t = line.strip()
        if t == "Coordinates":
            in_coordinates = True
            continue
        if t == "End Coordinates":
            break
        if not in_coordinates or not t:
            continue
        parts = t.split()
        if int(parts[0]) == node_id:
            got = tuple(float(v) for v in parts[1:4])
            exp = (ex, ey, ez)
            if any(abs(a - b) > tol for a, b in zip(got, exp)):
                raise SystemExit(f"node {node_id}: expected {exp}, got {got}")
            break
    else:
        raise SystemExit(f"node {node_id} not found in {mesh_file}")
PY
}

CYL_ORDER_B_STDOUT="$TMP_DIR/out_cylinder_order_b.stdout"
"$BIN" --in "$FIX" --out "$TMP_DIR/out_cylinder_order_b.post.msh" --op cylinder:1 --op translate:0,0,1 >"$CYL_ORDER_B_STDOUT"
rg -n "^mesh\.cylinder\.elements_kept=0$" "$CYL_ORDER_B_STDOUT"

ROT_Z90_OUT="$TMP_DIR/out_rotate_z90.post.msh"
"$BIN" --in "$FIX" --out "$ROT_Z90_OUT" --op rotate:0,0,1,90 >/tmp/meshpp_rotate_z90_stdout.txt
assert_node_xyz "$ROT_Z90_OUT" 2 0 1 0
assert_node_xyz "$ROT_Z90_OUT" 3 -1 1 0
assert_node_xyz "$ROT_Z90_OUT" 7 -1 1 1
"$ROUNDTRIP" "$ROT_Z90_OUT" "$TMP_DIR/out_rotate_z90_roundtrip.post.msh" --validate >"$TMP_DIR/out_rotate_z90.validate"
rg -q "nodes: 8" "$TMP_DIR/out_rotate_z90.validate"
rg -q "elements: 1" "$TMP_DIR/out_rotate_z90.validate"

ROT_X180_OUT="$TMP_DIR/out_rotate_x180.post.msh"
"$BIN" --in "$FIX" --out "$ROT_X180_OUT" --op rotate:1,0,0,180 >/tmp/meshpp_rotate_x180_stdout.txt
assert_node_xyz "$ROT_X180_OUT" 3 1 -1 0
assert_node_xyz "$ROT_X180_OUT" 7 1 -1 -1

ROT_Z90_NONUNIT_OUT="$TMP_DIR/out_rotate_z90_nonunit.post.msh"
"$BIN" --in "$FIX" --out "$ROT_Z90_NONUNIT_OUT" --op rotate:0,0,10,90 >/tmp/meshpp_rotate_z90_nonunit_stdout.txt
assert_node_xyz "$ROT_Z90_NONUNIT_OUT" 2 0 1 0
assert_node_xyz "$ROT_Z90_NONUNIT_OUT" 3 -1 1 0

ROT_ZERO_OUT="$TMP_DIR/out_rotate_zero.post.msh"
"$BIN" --in "$FIX" --out "$ROT_ZERO_OUT" --op rotate:0,0,1,0 >/tmp/meshpp_rotate_zero_stdout.txt
assert_node_xyz "$ROT_ZERO_OUT" 2 1 0 0
assert_node_xyz "$ROT_ZERO_OUT" 7 1 1 1

ROT_ORDER_STATS_OUT="$TMP_DIR/rotate_order_stats_stdout.txt"
"$BIN" --in "$FIX" --out "$TMP_DIR/out_rotate_order_stats.post.msh" --op translate:1,0,0 --op rotate:0,0,1,90 --mesh_stats >"$ROT_ORDER_STATS_OUT"
rg -n "^mesh\.stats\.min\.x=-1\.000000$" "$ROT_ORDER_STATS_OUT"
rg -n "^mesh\.stats\.max\.x=0\.000000$" "$ROT_ORDER_STATS_OUT"
rg -n "^mesh\.stats\.min\.y=1\.000000$" "$ROT_ORDER_STATS_OUT"
rg -n "^mesh\.stats\.max\.y=2\.000000$" "$ROT_ORDER_STATS_OUT"


RECT_FIX="$TMP_DIR/rectangular_prism.post.msh"
cat >"$RECT_FIX" <<'EOF'
MESH "rect" dimension 3 ElemType Hexahedra Nnode 8

Coordinates
1 0 0 0
2 1 0 0
3 1 2 0
4 0 2 0
5 0 0 4
6 1 0 4
7 1 2 4
8 0 2 4
End Coordinates
Elements
1 1 2 3 4 5 6 7 8
End Elements
EOF

ALIGN_ROTATED_OUT="$TMP_DIR/out_align_axes_rotated.post.msh"
ALIGN_STDOUT="$TMP_DIR/out_align_axes.stdout"
"$BIN" --in "$RECT_FIX" --out "$ALIGN_ROTATED_OUT" --op rotate:0,0,1,30 --op align_axes:pca --mesh_stats >"$ALIGN_STDOUT"
rg -n "^mesh\.align_axes\.method=pca$" "$ALIGN_STDOUT"
rg -n "^mesh\.align_axes\.eigenvalue\.0=" "$ALIGN_STDOUT"
rg -n "^mesh\.align_axes\.eigenvector\.0\.x=" "$ALIGN_STDOUT"
rg -n "^mesh\.align_axes\.eigenvector\.2\.z=" "$ALIGN_STDOUT"
rg -n "^mesh\.align_axes\.matrix\.r00=" "$ALIGN_STDOUT"
rg -n "^mesh\.align_axes\.matrix\.r22=" "$ALIGN_STDOUT"
METHOD_LINE=$(rg -n "^mesh\.align_axes\.method=pca$" "$ALIGN_STDOUT" | cut -d: -f1)
EIGEN_LINE=$(rg -n "^mesh\.align_axes\.eigenvalue\.0=" "$ALIGN_STDOUT" | cut -d: -f1)
MATRIX_LINE=$(rg -n "^mesh\.align_axes\.matrix\.r00=" "$ALIGN_STDOUT" | cut -d: -f1)
if (( METHOD_LINE >= EIGEN_LINE || EIGEN_LINE >= MATRIX_LINE )); then
  echo "align_axes PCA report order changed" >&2
  exit 1
fi
rg -n "^mesh\.stats\.nodes=8$" "$ALIGN_STDOUT"
rg -n "^mesh\.stats\.elements=1$" "$ALIGN_STDOUT"
rg -n "^mesh\.stats\.bbox\.dx=2\.000000$" "$ALIGN_STDOUT"
rg -n "^mesh\.stats\.bbox\.dy=1\.000000$" "$ALIGN_STDOUT"
rg -n "^mesh\.stats\.bbox\.dz=4\.000000$" "$ALIGN_STDOUT"
"$ROUNDTRIP" "$ALIGN_ROTATED_OUT" "$TMP_DIR/out_align_axes_roundtrip.post.msh" --validate >"$TMP_DIR/out_align_axes.validate"
rg -q "nodes: 8" "$TMP_DIR/out_align_axes.validate"
rg -q "elements: 1" "$TMP_DIR/out_align_axes.validate"

ALIGN_SIMPLE_FIX="$TMP_DIR/align_axes_simple_many.post.msh"
python3 - <<'PY' >"$ALIGN_SIMPLE_FIX"
print('MESH "many" dimension 3 ElemType Hexahedra Nnode 8')
print()
print('Coordinates')
node_id = 1
elements = []
for element_id in range(1, 13):
    x0 = element_id * 2
    ids = []
    for dx, dy, dz in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)):
        ids.append(node_id)
        print(node_id, x0 + dx, dy, dz)
        node_id += 1
    elements.append((element_id, ids))
print('End Coordinates')
print('Elements')
for element_id, ids in elements:
    print(element_id, *ids)
print('End Elements')
PY

ALIGN_SIMPLE_OUT="$TMP_DIR/out_align_axes_simple.post.msh"
ALIGN_SIMPLE_STDOUT="$TMP_DIR/out_align_axes_simple.stdout"
"$BIN" --in "$ALIGN_SIMPLE_FIX" --out "$ALIGN_SIMPLE_OUT" --op align_axes:simple --mesh_stats >"$ALIGN_SIMPLE_STDOUT"
rg -n "^mesh\.align_axes\.method=simple$" "$ALIGN_SIMPLE_STDOUT"
rg -n "^mesh\.align_axes\.elements\.exported=10$" "$ALIGN_SIMPLE_STDOUT"
rg -n "^mesh\.align_axes\.elements\.dropped=2$" "$ALIGN_SIMPLE_STDOUT"
rg -n "^mesh\.stats\.nodes=80$" "$ALIGN_SIMPLE_STDOUT"
rg -n "^mesh\.stats\.elements=10$" "$ALIGN_SIMPLE_STDOUT"
"$ROUNDTRIP" "$ALIGN_SIMPLE_OUT" "$TMP_DIR/out_align_axes_simple_roundtrip.post.msh" --validate >"$TMP_DIR/out_align_axes_simple.validate"
rg -q "nodes: 80" "$TMP_DIR/out_align_axes_simple.validate"
rg -q "elements: 10" "$TMP_DIR/out_align_axes_simple.validate"
if awk '/^Elements$/{in_elements=1; next} /^End Elements$/{in_elements=0} in_elements && ($1 == 11 || $1 == 12){found=1} END{exit found ? 0 : 1}' "$ALIGN_SIMPLE_OUT"; then
  echo "align_axes:simple exported more than ten elements" >&2
  exit 1
fi

FIRST_LINE=$(sed -n "1p" "$STATS_OUT")
SECOND_LINE=$(sed -n "2p" "$STATS_OUT")
THIRD_LINE=$(sed -n "3p" "$STATS_OUT")
if [[ "$FIRST_LINE" != "mesh.stats.nodes=8" || "$SECOND_LINE" != "mesh.stats.elements=1" || "$THIRD_LINE" != "mesh.stats.min.x=0.000000" ]]; then
  echo "stats output order changed" >&2
  cat "$STATS_OUT" >&2
  exit 1
fi


echo "meshpp operation pipeline checks: PASS"
