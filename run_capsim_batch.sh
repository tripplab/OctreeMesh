#!/bin/bash
set -u

usage() {
  cat <<USAGE
Usage: $0 --pdb LIST --folds LIST --res SPEC --threads N --bin PATH --vdb-dir DIR [--steps SPEC] [--cleanup-checkpoints MODE] [--mesh-orientation MODE] [--patch-radius R] [--strict-skips] [--smoke] [--run-parallel-jobs N]

Required:
  --pdb         Comma list of pdb IDs (allowed: 1cwp,3j4u,3izg,4g93)
  --folds       Comma list from: 2_0,2_1,3_0,3_1,5_0
  --res         Resolution spec (e.g. 1-16 or 1-4,8,10-12)
  --threads     Solver threads [1-30]
  --bin         Path to OctreeMesh bin directory
  --vdb-dir     Directory where <vdb>.vdb is expected first (fallback: cwd)
Optional:
  --steps SPEC      Steps passed to run_capsim.sh (default: 1-5; e.g. 1-5 or 2,3,4)
  --cleanup-checkpoints MODE
                    Checkpoint cleanup mode passed to run_capsim.sh: yes or no (default: no)
  --mesh-orientation MODE
                    octree_mesh orientation mode passed to run_capsim.sh:
                    --rotate_meshed_atoms/rotate_meshed_atoms or --mesh_rotated_atoms/mesh_rotated_atoms
                    (default: --rotate_meshed_atoms)
  --patch-radius R  Patch radius in Å; computes cone angle per PDB diameter
  --strict-skips    Missing VDB skips trigger non-zero exit
  --smoke           Build/validate/summarize only, do not run simulations
  --run-parallel-jobs N
                    Run generated jobs concurrently as detached background processes,
                    with at most N jobs running at once
USAGE
  exit 1
}

for req in awk sed date tr nohup; do
  command -v "$req" >/dev/null 2>&1 || { echo "Missing required tool: $req"; exit 1; }
done

ORIG_ARGS=("$@")

PDB_LIST=""; FOLD_LIST=""; RES_SPEC=""; THREADS=""; BIN_PATH=""; VDB_DIR=""; STEPS_SPEC="1-5"; CLEANUP_CHECKPOINTS="no"; MESH_ORIENTATION_MODE="--rotate_meshed_atoms"; PATCH_RADIUS=""; STRICT_SKIPS=0; SMOKE=0; PARALLEL_JOBS=0; PARALLEL_JOBS_SET=0; INTERNAL_WORKER_SPEC=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pdb) PDB_LIST="$2"; shift 2 ;;
    --pdb=*) PDB_LIST="${1#*=}"; shift ;;
    --folds) FOLD_LIST="$2"; shift 2 ;;
    --folds=*) FOLD_LIST="${1#*=}"; shift ;;
    --res) RES_SPEC="$2"; shift 2 ;;
    --res=*) RES_SPEC="${1#*=}"; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --threads=*) THREADS="${1#*=}"; shift ;;
    --bin) BIN_PATH="$2"; shift 2 ;;
    --bin=*) BIN_PATH="${1#*=}"; shift ;;
    --vdb-dir) VDB_DIR="$2"; shift 2 ;;
    --vdb-dir=*) VDB_DIR="${1#*=}"; shift ;;
    --steps) STEPS_SPEC="$2"; shift 2 ;;
    --steps=*) STEPS_SPEC="${1#*=}"; shift ;;
    --cleanup-checkpoints) CLEANUP_CHECKPOINTS="$2"; shift 2 ;;
    --cleanup-checkpoints=*) CLEANUP_CHECKPOINTS="${1#*=}"; shift ;;
    --mesh-orientation) MESH_ORIENTATION_MODE="$2"; shift 2 ;;
    --mesh-orientation=*) MESH_ORIENTATION_MODE="${1#*=}"; shift ;;
    --patch-radius) PATCH_RADIUS="$2"; shift 2 ;;
    --patch-radius=*) PATCH_RADIUS="${1#*=}"; shift ;;
    --strict-skips) STRICT_SKIPS=1; shift ;;
    --smoke) SMOKE=1; shift ;;
    --run-parallel-jobs) [[ $# -ge 2 ]] || { echo "--run-parallel-jobs requires N"; usage; }; PARALLEL_JOBS="$2"; PARALLEL_JOBS_SET=1; shift 2 ;;
    --run-parallel-jobs=*) PARALLEL_JOBS="${1#*=}"; PARALLEL_JOBS_SET=1; shift ;;
    --_capsim-batch-job-worker) [[ $# -ge 2 ]] || { echo "--_capsim-batch-job-worker requires SPEC"; exit 1; }; INTERNAL_WORKER_SPEC="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$INTERNAL_WORKER_SPEC" ]]; then
  [[ -n "$PDB_LIST" && -n "$FOLD_LIST" && -n "$RES_SPEC" && -n "$THREADS" && -n "$BIN_PATH" && -n "$VDB_DIR" ]] || usage
fi
if [[ -z "$INTERNAL_WORKER_SPEC" ]]; then
  [[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "--threads must be integer"; exit 1; }
  (( THREADS >= 1 && THREADS <= 30 )) || { echo "--threads must be in [1,30]"; exit 1; }
  if (( PARALLEL_JOBS_SET == 1 )); then
    [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || { echo "--run-parallel-jobs must be a positive integer"; exit 1; }
    (( PARALLEL_JOBS >= 1 )) || { echo "--run-parallel-jobs must be >= 1"; exit 1; }
  fi
fi
if [[ -z "$INTERNAL_WORKER_SPEC" ]]; then
[[ "$STEPS_SPEC" =~ ^[0-9]+-[0-9]+$ || "$STEPS_SPEC" =~ ^[0-9]+(,[0-9]+)*$ ]] || { echo "--steps must be a range (e.g. 1-5) or comma list (e.g. 2,3,4)"; exit 1; }
if [[ "$STEPS_SPEC" == *-* ]]; then
  step_start="${STEPS_SPEC%-*}"
  step_end="${STEPS_SPEC#*-}"
  (( step_start >= 1 && step_start <= step_end && step_end <= 5 )) || { echo "--steps range must be within [1,5]"; exit 1; }
else
  IFS=',' read -r -a step_arr <<< "$STEPS_SPEC"
  for step_i in "${step_arr[@]}"; do
    (( step_i >= 1 && step_i <= 5 )) || { echo "--steps entries must be within [1,5]"; exit 1; }
  done
fi
case "$CLEANUP_CHECKPOINTS" in
  yes|no) ;;
  *) echo "--cleanup-checkpoints must be yes or no"; exit 1 ;;
esac
case "$MESH_ORIENTATION_MODE" in
  --rotate_meshed_atoms|--mesh_rotated_atoms) ;;
  rotate_meshed_atoms) MESH_ORIENTATION_MODE="--rotate_meshed_atoms" ;;
  mesh_rotated_atoms) MESH_ORIENTATION_MODE="--mesh_rotated_atoms" ;;
  *) echo "--mesh-orientation must be --rotate_meshed_atoms or --mesh_rotated_atoms"; exit 1 ;;
esac
[[ -d "$BIN_PATH" ]] || { echo "--bin directory not found: $BIN_PATH"; exit 1; }
[[ -d "$VDB_DIR" ]] || { echo "--vdb-dir not found: $VDB_DIR"; exit 1; }
if [[ -n "$PATCH_RADIUS" ]]; then
  awk -v v="$PATCH_RADIUS" 'BEGIN { exit ((v+0 > 0 && v !~ /[^0-9.]/ && v !~ /\..*\./) ? 0 : 1) }' || { echo "--patch-radius must be a positive number in Å"; exit 1; }
fi

fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_CONFIG="$SCRIPT_DIR/capsim_config.sh"
RUN_SCRIPT="$SCRIPT_DIR/run_capsim.sh"
if [[ -z "$INTERNAL_WORKER_SPEC" ]]; then
  [[ -f "$BASE_CONFIG" ]] || { echo "Missing base config: $BASE_CONFIG"; exit 1; }
  [[ -x "$RUN_SCRIPT" || -f "$RUN_SCRIPT" ]] || { echo "Missing run script: $RUN_SCRIPT"; exit 1; }
fi

# validate required bin entries inferred from run_capsim.sh
# Step 5 is named "rotate_back" but implemented via apply-matrix.awk.
if [[ -z "$INTERNAL_WORKER_SPEC" ]]; then
  for b in extract_ATOM octree_mesh meshsolver mesh2pdb apply-matrix.awk; do
    [[ -e "$BIN_PATH/$b" ]] || { echo "Missing required BIN entry: $BIN_PATH/$b"; exit 1; }
  done
fi

# detach once
if [[ -z "$INTERNAL_WORKER_SPEC" && "${CAPSIM_BATCH_CHILD:-0}" != "1" ]]; then
  ts="$(date -u +%Y%m%dT%H%M%SZ)"
  mkdir -p "$SCRIPT_DIR/runs"
  master_log="$SCRIPT_DIR/runs/batch_${ts}_$$.log"
  CAPSIM_BATCH_CHILD=1 nohup "$0" "${ORIG_ARGS[@]}" > "$master_log" 2>&1 &
  pid=$!
  echo "Spawned batch PID: $pid"
  echo "Master log: $master_log"
  echo "Tail command: tail -f $master_log"
  exit 0
fi

set -e
if [[ -z "$INTERNAL_WORKER_SPEC" ]]; then
  mkdir -p "$SCRIPT_DIR/runs"
  ts="$(date -u +%Y%m%dT%H%M%SZ)"
  out_prefix="$SCRIPT_DIR/runs/batch_${ts}_$$"
  tsv_file="${out_prefix}.tsv"
  csv_file="${out_prefix}.csv"
  [[ ! -e "$tsv_file" && ! -e "$csv_file" ]] || { echo "Summary files already exist"; exit 1; }
  work_dir="$SCRIPT_DIR/runs/batch_${ts}_$$"
  mkdir -p "$work_dir/tmp_configs"

  echo -e "status\texit_code\truntime_sec\tpdb\tvdb\tres\tyoung\tfold_type\tfold_index\tthreads\tsteps\tmesh_orientation_mode\tpatch_radius\tcapsid_diameter\tcone_deg\trun_dir" > "$tsv_file"
  echo "job_name,total_proteins,total_atoms,nodes,elements,mesh_volume,volume_loaded,octree_mesh_sec,meshsolver_sec,mesh2pdb_sec" > "$csv_file"
fi

young_for_pdb() {
  case "$1" in
    1cwp) echo "0.0193" ;;
    3j4u) echo "0.0096" ;;
    3izg) echo "0.0245" ;;
    4g93) echo "0.0147" ;;
    *) return 1 ;;
  esac
}

diameter_for_pdb() {
  case "$1" in
    1cwp) echo "280" ;;
    4g93) echo "344" ;;
    3izg) echo "527" ;;
    3j4u) echo "631" ;;
    *) return 1 ;;
  esac
}

cone_for_patch_radius() {
  local patch_radius="$1"
  local capsid_diameter="$2"
  awk -v r="$patch_radius" -v d="$capsid_diameter" '
    BEGIN {
      x = 2*r/d
      if (x <= 0 || x > 1) exit 2
      pi = atan2(0, -1)
      theta = atan2(x, sqrt(1 - x*x)) * 180 / pi
      printf "%.6f\n", theta
    }
  '
}

normalize_num() { echo "$1" | tr -d ','; }

parse_metric() {
  local file="$1"; shift
  local sed_expr="$1"
  local v
  v=$(sed -nE "$sed_expr" "$file" | tail -n1 || true)
  if [[ -z "$v" ]]; then echo "NA"; else normalize_num "$v"; fi
}

step_timing() {
  local run_dir="$1"; local step="$2"
  local v
  v=$(awk -F'\t' -v rd="$run_dir" -v st="$step" '
    NR>1 && $2==st && $3=="done" {
      key=$1
      sub(/^\.\//, "", key)
      sub(/^\.\//, "", rd)
      short=rd
      sub(/^.*\/runs\//, "runs/", short)
      if (key==rd || key==short) val=$4
    }
    END{if(val!="") print val}
  ' "$SCRIPT_DIR/.checkpoints/timings.ts" 2>/dev/null || true)
  [[ -n "$v" ]] && echo "$v" || echo "NA"
}

shell_quote() { printf "%q" "$1"; }

write_spec_var() {
  local spec_file="$1"; local name="$2"; local value="$3"
  printf "%s=%s\n" "$name" "$(shell_quote "$value")" >> "$spec_file"
}

run_job_from_spec() {
  local spec_file="$1"
  # shellcheck disable=SC1090
  source "$spec_file"

  local job_start run_dir status exit_code runtime job_name
  local tp ta nodes elems mv vl t2 t3 t4 src_log rc fail_step
  job_start=$(date +%s)
  run_dir="NA"
  status="SMOKE_OK"
  exit_code=0
  runtime="NA"
  job_name="NA"
  tp="NA"; ta="NA"; nodes="NA"; elems="NA"; mv="NA"; vl="NA"; t2="NA"; t3="NA"; t4="NA"

  if [[ ! -f "$vdb_path" ]]; then
    status="SKIPPED_MISSING_VDB"
    exit_code=0
    runtime="NA"
  else
    cp "$BASE_CONFIG" "$cfg"
    sed -i -E "s|^PDB=.*$|PDB=${pdb_u}|" "$cfg"
    sed -i -E "s|^VDB=.*$|VDB=${vdb}|" "$cfg"
    sed -i -E "s|^Res=.*$|Res=${res_f}|" "$cfg"
    sed -i -E "s|^Young=.*$|Young=${young}|" "$cfg"
    sed -i -E "s|^FOLD_TYPE=.*$|FOLD_TYPE=${fold_type}|" "$cfg"
    sed -i -E "s|^FOLD_INDEX=.*$|FOLD_INDEX=${fold_index}|" "$cfg"
    sed -i -E "s|^SOLVER_THREADS=.*$|SOLVER_THREADS=${THREADS}|" "$cfg"
    sed -i -E "s|^BIN=.*$|BIN=${BIN_PATH}|" "$cfg"
    if [[ "$patch_radius_tsv" != "NA" ]]; then
      sed -i -E "s|^cone=.*$|cone=${cone_deg_tsv}|" "$cfg"
    fi

    if (( SMOKE == 0 )); then
      run_dir="$SCRIPT_DIR/runs/$run_tag"
      set +e
      RUN_TAG="$run_tag" "$RUN_SCRIPT" -c "$cfg" -t "$THREADS" --steps "$STEPS_SPEC" --cleanup-checkpoints "$CLEANUP_CHECKPOINTS" --mesh-orientation "$MESH_ORIENTATION_MODE" > "$tmp_log" 2>&1
      rc=$?
      set -e
      if [[ -d "$run_dir" ]]; then
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
      if [[ "$status" == "DONE" && "$tp" == "NA" && "$ta" == "NA" && "$nodes" == "NA" && "$elems" == "NA" && "$mv" == "NA" && "$vl" == "NA" ]]; then
        echo "WARN parse_metric: all structural metrics are NA for job=$job_name log=$src_log"
      fi
      if [[ "$status" == "DONE" && "$run_dir" != "NA" && "$t2" == "NA" && "$t3" == "NA" && "$t4" == "NA" ]]; then
        echo "WARN step_timing: all timing metrics are NA for job=$job_name run_dir=$run_dir"
      fi
    else
      runtime=0
      status="SMOKE_OK"
      job_name="${pdb_u}_F${fold_type}id${fold_index}_R${res_f}_SMOKE"
    fi
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$status" "$exit_code" "$runtime" "$pdb_u" "$vdb" "$res_f" "$young" "$fold_type" "$fold_index" "$THREADS" "$STEPS_SPEC" "$MESH_ORIENTATION_MODE" "$patch_radius_tsv" "$capsid_diameter_tsv" "$cone_deg_tsv" "$run_dir" > "$result_tsv"
  echo "${job_name},${tp},${ta},${nodes},${elems},${mv},${vl},${t2},${t3},${t4}" > "$result_csv"
  printf "%s\n%s\n" "$status" "$exit_code" > "$result_status"
  [[ "$status" == FAILED_STEP_* ]] && return 1
  return 0
}

if [[ -n "$INTERNAL_WORKER_SPEC" ]]; then
  set -e
  run_job_from_spec "$INTERNAL_WORKER_SPEC"
  exit $?
fi

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

job_ids=()
job_specs=()
job_count=0
fail_count=0
skip_count=0
mkdir -p "$work_dir/job_specs" "$work_dir/job_logs" "$work_dir/results"

batch_id="$(date -u +%Y%m%dT%H%M%SZ)_$$"

# Pre-create timing metadata to avoid concurrent header initialization in run_capsim.sh.
mkdir -p "$SCRIPT_DIR/.checkpoints"
if [[ ! -f "$SCRIPT_DIR/.checkpoints/timings.ts" ]]; then
  printf "run_dir\tstep_name\tstatus\telapsed_sec\tsolver_threads\tshear_mode\tnote\n" > "$SCRIPT_DIR/.checkpoints/timings.ts"
fi

for pdb_l in "${pdb_arr[@]}"; do
  pdb_u="$(echo "$pdb_l" | tr '[:lower:]' '[:upper:]')"
  vdb="${pdb_l}_full"
  young="$(young_for_pdb "$pdb_l")"
  patch_radius_tsv="NA"
  capsid_diameter_tsv="NA"
  cone_deg_tsv="NA"
  if [[ -n "$PATCH_RADIUS" ]]; then
    capsid_diameter_tsv="$(diameter_for_pdb "$pdb_l")"
    awk -v r="$PATCH_RADIUS" -v d="$capsid_diameter_tsv" 'BEGIN { exit (2*r <= d ? 0 : 1) }' || {
      echo "--patch-radius $PATCH_RADIUS Å is too large for $pdb_l diameter $capsid_diameter_tsv Å; must be <= $(awk -v d="$capsid_diameter_tsv" 'BEGIN { printf "%.6g", d/2 }') Å"
      exit 1
    }
    cone_deg_tsv="$(cone_for_patch_radius "$PATCH_RADIUS" "$capsid_diameter_tsv")"
    patch_radius_tsv="$PATCH_RADIUS"
  fi
  for fold in "${fold_arr[@]}"; do
    fold_type="${fold%_*}"; fold_index="${fold#*_}"
    for res_i in "${uniq_res[@]}"; do
      res_f=$(printf "%.2f" "$res_i")
      ((++job_count))
      job_id=$(printf "%06d" "$job_count")
      job_key="${pdb_u}_F${fold_type}_${fold_index}_R${res_f}"
      cfg="$work_dir/tmp_configs/${job_key}.sh"
      tmp_log="$work_dir/job_logs/${job_id}_${job_key}.log"
      result_tsv="$work_dir/results/${job_id}.tsv"
      result_csv="$work_dir/results/${job_id}.csv"
      result_status="$work_dir/results/${job_id}.status"
      run_tag="${pdb_u}_F${fold_type}id${fold_index}_R${res_f}_B${batch_id}_J${job_id}"

      vdb_path="$VDB_DIR/${vdb}.vdb"
      if [[ ! -f "$vdb_path" ]]; then
        vdb_path="$SCRIPT_DIR/${vdb}.vdb"
      fi

      spec_file="$work_dir/job_specs/${job_id}.spec"
      : > "$spec_file"
      for var in SCRIPT_DIR BASE_CONFIG RUN_SCRIPT THREADS BIN_PATH STEPS_SPEC CLEANUP_CHECKPOINTS MESH_ORIENTATION_MODE SMOKE cfg tmp_log result_tsv result_csv result_status pdb_u vdb res_f young fold_type fold_index patch_radius_tsv capsid_diameter_tsv cone_deg_tsv run_tag vdb_path; do
        write_spec_var "$spec_file" "$var" "${!var}"
      done
      job_ids+=("$job_id")
      job_specs+=("$spec_file")
    done
  done
done

echo "Planned jobs: $job_count"
if (( PARALLEL_JOBS > 0 )); then
  echo "Parallel job mode enabled: max concurrent jobs = $PARALLEL_JOBS"
else
  echo "Sequential job mode enabled"
fi

run_spec_sequential() {
  local spec_file="$1"
  set +e
  "$0" --_capsim-batch-job-worker "$spec_file"
  local rc=$?
  set -e
  return "$rc"
}

launch_spec_parallel() {
  local spec_file="$1"; local job_id="$2"
  local worker_log="$work_dir/job_logs/${job_id}.worker.nohup.log"
  CAPSIM_BATCH_CHILD=1 nohup "$0" --_capsim-batch-job-worker "$spec_file" > "$worker_log" 2>&1 &
  LAUNCHED_PID=$!
}

if (( PARALLEL_JOBS > 0 )); then
  active_pids=()
  active_ids=()
  completed_jobs=0
  for i in "${!job_specs[@]}"; do
    while (( ${#active_pids[@]} >= PARALLEL_JOBS )); do
      pid="${active_pids[0]}"
      job_id="${active_ids[0]}"
      wait "$pid" || true
      active_pids=("${active_pids[@]:1}")
      active_ids=("${active_ids[@]:1}")
      ((++completed_jobs))
      echo "Completed job $job_id ($completed_jobs/$job_count)"
    done
    job_id="${job_ids[$i]}"
    spec_file="${job_specs[$i]}"
    launch_spec_parallel "$spec_file" "$job_id"
    pid="$LAUNCHED_PID"
    active_pids+=("$pid")
    active_ids+=("$job_id")
    echo "[$((i + 1))/$job_count] spawned PID $pid for job $job_id"
  done

  while (( ${#active_pids[@]} > 0 )); do
    pid="${active_pids[0]}"
    job_id="${active_ids[0]}"
    wait "$pid" || true
    active_pids=("${active_pids[@]:1}")
    active_ids=("${active_ids[@]:1}")
    ((++completed_jobs))
    echo "Completed job $job_id ($completed_jobs/$job_count)"
  done
else
  for i in "${!job_specs[@]}"; do
    echo "[$((i + 1))/$job_count] running job ${job_ids[$i]}"
    run_spec_sequential "${job_specs[$i]}" || true
  done
fi

for job_id in "${job_ids[@]}"; do
  result_tsv="$work_dir/results/${job_id}.tsv"
  result_csv="$work_dir/results/${job_id}.csv"
  result_status="$work_dir/results/${job_id}.status"
  if [[ ! -f "$result_tsv" || ! -f "$result_csv" || ! -f "$result_status" ]]; then
    echo "Missing result for job $job_id"
    ((++fail_count))
    continue
  fi
  cat "$result_tsv" >> "$tsv_file"
  cat "$result_csv" >> "$csv_file"
  status=$(sed -n '1p' "$result_status")
  if [[ "$status" == "SKIPPED_MISSING_VDB" ]]; then
    ((++skip_count))
  elif [[ "$status" == FAILED_STEP_* ]]; then
    ((++fail_count))
  fi
done

echo "Summary TSV: $tsv_file"
echo "Summary CSV: $csv_file"
echo "Finished jobs: $job_count, failed: $fail_count, skipped: $skip_count"

if (( fail_count > 0 )); then
  exit 1
fi
if (( STRICT_SKIPS == 1 && skip_count > 0 )); then
  exit 1
fi
exit 0
