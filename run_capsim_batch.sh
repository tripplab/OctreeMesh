#!/bin/bash
set -u

usage() {
  cat <<USAGE
Usage: $0 --pdb LIST --folds LIST --res SPEC --threads N --bin PATH --vdb-dir DIR [--strict-skips] [--smoke]

Required:
  --pdb         Comma list of pdb IDs (allowed: 1cwp,3j4u,3izg,4g93)
  --folds       Comma list from: 2_0,2_1,3_0,3_1,5_0
  --res         Resolution spec (e.g. 1-16 or 1-4,8,10-12)
  --threads     Solver threads [1-30]
  --bin         Path to OctreeMesh bin directory
  --vdb-dir     Directory where <vdb>.vdb is expected first (fallback: cwd)
Optional:
  --strict-skips  Missing VDB skips trigger non-zero exit
  --smoke         Build/validate/summarize only, do not run simulations
USAGE
  exit 1
}

for req in awk sed date tr nohup; do
  command -v "$req" >/dev/null 2>&1 || { echo "Missing required tool: $req"; exit 1; }
done

PDB_LIST=""; FOLD_LIST=""; RES_SPEC=""; THREADS=""; BIN_PATH=""; VDB_DIR=""; STRICT_SKIPS=0; SMOKE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pdb) PDB_LIST="$2"; shift 2 ;;
    --folds) FOLD_LIST="$2"; shift 2 ;;
    --res) RES_SPEC="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --bin) BIN_PATH="$2"; shift 2 ;;
    --vdb-dir) VDB_DIR="$2"; shift 2 ;;
    --strict-skips) STRICT_SKIPS=1; shift ;;
    --smoke) SMOKE=1; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

[[ -n "$PDB_LIST" && -n "$FOLD_LIST" && -n "$RES_SPEC" && -n "$THREADS" && -n "$BIN_PATH" && -n "$VDB_DIR" ]] || usage
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "--threads must be integer"; exit 1; }
(( THREADS >= 1 && THREADS <= 30 )) || { echo "--threads must be in [1,30]"; exit 1; }
[[ -d "$BIN_PATH" ]] || { echo "--bin directory not found: $BIN_PATH"; exit 1; }
[[ -d "$VDB_DIR" ]] || { echo "--vdb-dir not found: $VDB_DIR"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_CONFIG="$SCRIPT_DIR/capsim_config.sh"
RUN_SCRIPT="$SCRIPT_DIR/run_capsim.sh"
[[ -f "$BASE_CONFIG" ]] || { echo "Missing base config: $BASE_CONFIG"; exit 1; }
[[ -x "$RUN_SCRIPT" || -f "$RUN_SCRIPT" ]] || { echo "Missing run script: $RUN_SCRIPT"; exit 1; }

# validate required bin entries inferred from run_capsim.sh
for b in extract_ATOM octree_mesh meshsolver mesh2pdb apply-matrix.awk rotate_back; do
  [[ -e "$BIN_PATH/$b" ]] || { echo "Missing required BIN entry: $BIN_PATH/$b"; exit 1; }
done

# detach once
if [[ "${CAPSIM_BATCH_CHILD:-0}" != "1" ]]; then
  ts="$(date -u +%Y%m%dT%H%M%SZ)"
  mkdir -p "$SCRIPT_DIR/runs"
  master_log="$SCRIPT_DIR/runs/batch_${ts}_$$.log"
  CAPSIM_BATCH_CHILD=1 nohup "$0" "$@" > "$master_log" 2>&1 &
  pid=$!
  echo "Spawned batch PID: $pid"
  echo "Master log: $master_log"
  echo "Tail command: tail -f $master_log"
  exit 0
fi

set -e
mkdir -p "$SCRIPT_DIR/runs"
ts="$(date -u +%Y%m%dT%H%M%SZ)"
out_prefix="$SCRIPT_DIR/runs/batch_${ts}_$$"
tsv_file="${out_prefix}.tsv"
csv_file="${out_prefix}.csv"
[[ ! -e "$tsv_file" && ! -e "$csv_file" ]] || { echo "Summary files already exist"; exit 1; }
work_dir="$SCRIPT_DIR/runs/batch_${ts}_$$"
mkdir -p "$work_dir/tmp_configs"

echo -e "status\texit_code\truntime_sec\tpdb\tvdb\tres\tyoung\tfold_type\tfold_index\tthreads\trun_dir" > "$tsv_file"
echo "job_name,total_proteins,total_atoms,nodes,elements,mesh_volume,volume_loaded,octree_mesh_sec,meshsolver_sec,mesh2pdb_sec" > "$csv_file"

young_for_pdb() {
  case "$1" in
    1cwp) echo "0.0193" ;;
    3j4u) echo "0.0096" ;;
    3izg) echo "0.0245" ;;
    4g93) echo "0.0147" ;;
    *) return 1 ;;
  esac
}

normalize_num() { echo "$1" | tr -d ','; }

parse_metric() {
  local file="$1"; shift
  local regex="$1"
  local v
  v=$(sed -nE "s/${regex}/\\1/p" "$file" | tail -n1 || true)
  if [[ -z "$v" ]]; then echo "NA"; else normalize_num "$v"; fi
}

step_timing() {
  local run_dir="$1"; local step="$2"
  local v
  v=$(awk -F'\t' -v rd="$run_dir" -v st="$step" 'NR>1 && $1==rd && $2==st && $3=="done" {val=$4} END{if(val!="") print val}' "$SCRIPT_DIR/.checkpoints/timings.ts" 2>/dev/null || true)
  [[ -n "$v" ]] && echo "$v" || echo "NA"
}

# parse and validate lists
IFS=',' read -r -a pdb_arr <<< "$PDB_LIST"
valid_pdb=(1cwp 3j4u 3izg 4g93)
for i in "${!pdb_arr[@]}"; do
  p=$(echo "${pdb_arr[$i]}" | tr '[:upper:]' '[:lower:]' | xargs)
  ok=0
  for vp in "${valid_pdb[@]}"; do [[ "$p" == "$vp" ]] && ok=1; done
  (( ok == 1 )) || { echo "Invalid PDB: ${pdb_arr[$i]}"; exit 1; }
  pdb_arr[$i]="$p"
done

IFS=',' read -r -a fold_arr <<< "$FOLD_LIST"
for f in "${fold_arr[@]}"; do
  [[ "$f" =~ ^(2_0|2_1|3_0|3_1|5_0)$ ]] || { echo "Invalid fold token: $f"; exit 1; }
done

res_values=()
IFS=',' read -r -a res_parts <<< "$RES_SPEC"
for part in "${res_parts[@]}"; do
  part="$(echo "$part" | xargs)"
  if [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then
    a=${BASH_REMATCH[1]}; b=${BASH_REMATCH[2]}
    (( a<=b )) || { echo "Invalid res range: $part"; exit 1; }
    for ((r=a; r<=b; r++)); do res_values+=("$r"); done
  elif [[ "$part" =~ ^[0-9]+$ ]]; then
    res_values+=("$part")
  else
    echo "Invalid res token: $part"; exit 1
  fi
done
# dedupe preserve order and validate range
uniq_res=()
seen=""
for r in "${res_values[@]}"; do
  (( r>=1 && r<=16 )) || { echo "Resolution out of range [1,16]: $r"; exit 1; }
  if [[ ",$seen," != *",$r,"* ]]; then uniq_res+=("$r"); seen+="$r,"; fi
done

fail_count=0
skip_count=0

for pdb_l in "${pdb_arr[@]}"; do
  pdb_u="$(echo "$pdb_l" | tr '[:lower:]' '[:upper:]')"
  vdb="${pdb_l}_full"
  young="$(young_for_pdb "$pdb_l")"
  for fold in "${fold_arr[@]}"; do
    fold_type="${fold%_*}"; fold_index="${fold#*_}"
    for res_i in "${uniq_res[@]}"; do
      res_f=$(printf "%.2f" "$res_i")
      job_start=$(date +%s)
      run_dir="NA"
      status="SMOKE_OK"
      exit_code=0
      runtime="NA"
      job_name="NA"
      tp="NA"; ta="NA"; nodes="NA"; elems="NA"; mv="NA"; vl="NA"; t2="NA"; t3="NA"; t4="NA"

      vdb_path="$VDB_DIR/${vdb}.vdb"
      if [[ ! -f "$vdb_path" ]]; then
        vdb_path="$SCRIPT_DIR/${vdb}.vdb"
      fi
      if [[ ! -f "$vdb_path" ]]; then
        status="SKIPPED_MISSING_VDB"
        exit_code=0
        runtime="NA"
        ((skip_count++))
      else
        cfg="$work_dir/tmp_configs/${pdb_u}_F${fold_type}_${fold_index}_R${res_f}.sh"
        cp "$BASE_CONFIG" "$cfg"
        sed -i -E "s|^PDB=.*$|PDB=${pdb_u}|" "$cfg"
        sed -i -E "s|^VDB=.*$|VDB=${vdb}|" "$cfg"
        sed -i -E "s|^Res=.*$|Res=${res_f}|" "$cfg"
        sed -i -E "s|^Young=.*$|Young=${young}|" "$cfg"
        sed -i -E "s|^FOLD_TYPE=.*$|FOLD_TYPE=${fold_type}|" "$cfg"
        sed -i -E "s|^FOLD_INDEX=.*$|FOLD_INDEX=${fold_index}|" "$cfg"
        sed -i -E "s|^SOLVER_THREADS=.*$|SOLVER_THREADS=${THREADS}|" "$cfg"
        sed -i -E "s|^BIN=.*$|BIN=${BIN_PATH}|" "$cfg"

        if (( SMOKE == 0 )); then
          pre_ckpt=$(find "$SCRIPT_DIR/runs" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null | sort)
          tmp_log="$work_dir/current_job.log"
          set +e
          printf 'n\n' | "$RUN_SCRIPT" -c "$cfg" -t "$THREADS" > "$tmp_log" 2>&1
          rc=$?
          set -e
          post_ckpt=$(find "$SCRIPT_DIR/runs" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null | sort)
          run_tag=$(comm -13 <(echo "$pre_ckpt") <(echo "$post_ckpt") | tail -n1)
          if [[ -n "$run_tag" && -d "$SCRIPT_DIR/runs/$run_tag" ]]; then
            run_dir="$SCRIPT_DIR/runs/$run_tag"
            cp "$tmp_log" "$run_dir/batch_run.log"
            if [[ -f "$run_dir/manifest.json" ]]; then
              job_name=$(awk -F'"' '/"run_tag"/ {print $4}' "$run_dir/manifest.json" | head -n1)
            fi
          else
            run_dir="NA"
          fi
          [[ "$job_name" == "NA" || -z "$job_name" ]] && job_name="${pdb_u}_F${fold_type}id${fold_index}_R${res_f}_NA"

          runtime=$(awk -v s="$job_start" 'BEGIN{print systime()-s}')
          exit_code=$rc
          if (( rc != 0 )); then
            fail_step=$(sed -nE 's/.*✗ Step ([0-9]+).*/\1/p' "$tmp_log" | tail -n1)
            [[ -z "$fail_step" ]] && fail_step="NA"
            status="FAILED_STEP_${fail_step}"
            ((fail_count++))
          else
            status="DONE"
          fi

          src_log="$tmp_log"
          tp=$(parse_metric "$src_log" 's/.*[Tt]otal[ _-]*proteins[[:space:]]*[:=][[:space:]]*([^[:space:]]+).*/\1/p')
          ta=$(parse_metric "$src_log" 's/.*[Tt]otal[ _-]*atoms[[:space:]]*[:=][[:space:]]*([^[:space:]]+).*/\1/p')
          nodes=$(parse_metric "$src_log" 's/.*[# ]*[Nn]odes[[:space:]]*[:=][[:space:]]*([^[:space:]]+).*/\1/p')
          elems=$(parse_metric "$src_log" 's/.*[# ]*[Ee]lements[[:space:]]*[:=][[:space:]]*([^[:space:]]+).*/\1/p')
          mv=$(parse_metric "$src_log" 's/.*[Mm]esh[ _-]*[Vv]olume[[:space:]]*[:=][[:space:]]*([^[:space:]]+).*/\1/p')
          vl=$(parse_metric "$src_log" 's/.*[Vv]olume[ _-]*[Ll]oaded[[:space:]]*[:=][[:space:]]*([^[:space:]]+).*/\1/p')
          if [[ "$run_dir" != "NA" ]]; then
            t2=$(step_timing "$run_dir" "octree_mesh")
            t3=$(step_timing "$run_dir" "meshsolver")
            t4=$(step_timing "$run_dir" "mesh2pdb")
          fi
        else
          runtime=0
          status="SMOKE_OK"
          job_name="${pdb_u}_F${fold_type}id${fold_index}_R${res_f}_SMOKE"
        fi
      fi

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$status" "$exit_code" "$runtime" "$pdb_u" "$vdb" "$res_f" "$young" "$fold_type" "$fold_index" "$THREADS" "$run_dir" >> "$tsv_file"
      echo "${job_name},${tp},${ta},${nodes},${elems},${mv},${vl},${t2},${t3},${t4}" >> "$csv_file"
    done
  done
done

if (( fail_count > 0 )); then
  exit 1
fi
if (( STRICT_SKIPS == 1 && skip_count > 0 )); then
  exit 1
fi
exit 0
