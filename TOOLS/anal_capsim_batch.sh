#!/bin/bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $0 --work_dir DIR [OPTIONS]

Analyze Stage 1 capsim batch run directories and optionally extract completed
simulation result statistics into an accumulated CSV.

Expected run directory format:
  <capsid>_F<fold>id<fold_id>_R<resolution>_<time_stamp>
  <capsid>_F<fold>id<fold_id>_R<resolution>_S<seed>_<time_stamp>

Example run directory:
  3J4U_F5id0_R16.00_S0_20260511T041629Z

Required:
  --work_dir DIR              Directory containing run directories to scan

Reports:
  --tsv PATH                  Optional TSV report output path
  --csv PATH                  Optional CSV report output path

Simulation data extraction:
  --extract_sim_data          For completed runs, create/reuse per-run
                              extracted_sim_data.json and write an accumulated
                              CSV summary
  --extract_sim_data_csv PATH Output path for accumulated extracted data CSV
                              Default: <work_dir>/extracted_sim_data_summary.csv
  --force_extract_sim_data    Regenerate extracted_sim_data.json even if it
                              already exists

Behavior:
  A run is eligible for extraction only when both status and result are
  COMPLETED. The existing result=COMPLETED definition is preserved: non-empty
  octreemesh.post.res, octreemesh.post.msh, and octreemesh_solver.out must all
  exist. Directories named batch_* are ignored without malformed-name warnings.
  Extraction or JSON-read failures are written as NaN in the accumulated CSV.

Accumulated extraction CSV columns:
  run_dir,capsid,fold,fold_id,resolution,Displacement:magnitude:max

Other:
  --strict                    Exit non-zero if malformed entries or extraction
                              failures are encountered
  -h, --help                  Show this help

Examples:
  $0 --work_dir sims_constant_angle
  $0 --work_dir sims_constant_angle --csv runs.csv --tsv runs.tsv
  $0 --work_dir sims_constant_angle --extract_sim_data
  $0 --work_dir sims_constant_angle --extract_sim_data --extract_sim_data_csv sim_displacement_max.csv
  $0 --work_dir sims_constant_angle --extract_sim_data --force_extract_sim_data
USAGE
}
for req in find awk sed sort grep; do
  command -v "$req" >/dev/null 2>&1 || { echo "Missing required tool: $req"; exit 1; }
done

WORK_DIR=""
TSV_OUT=""
CSV_OUT=""
STRICT=0
EXTRACT_SIM_DATA=0
EXTRACT_SIM_DATA_CSV=""
FORCE_EXTRACT_SIM_DATA=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work_dir) WORK_DIR="${2:-}"; shift 2 ;;
    --tsv) TSV_OUT="${2:-}"; shift 2 ;;
    --csv) CSV_OUT="${2:-}"; shift 2 ;;
    --extract_sim_data) EXTRACT_SIM_DATA=1; shift ;;
    --extract_sim_data_csv) EXTRACT_SIM_DATA_CSV="${2:-}"; shift 2 ;;
    --force_extract_sim_data) FORCE_EXTRACT_SIM_DATA=1; shift ;;
    --strict) STRICT=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

[[ -n "$WORK_DIR" ]] || { echo "--work_dir is required"; usage; exit 1; }
[[ -d "$WORK_DIR" ]] || { echo "--work_dir not found: $WORK_DIR"; exit 1; }

WORK_DIR_ABS="$(cd "$WORK_DIR" && pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRACT_SCRIPT="$SCRIPT_DIR/extract_sim_data.py"
if (( EXTRACT_SIM_DATA == 1 )); then
  command -v python3 >/dev/null 2>&1 || { echo "Missing required tool: python3"; exit 1; }
  [[ -f "$EXTRACT_SCRIPT" ]] || { echo "extract_sim_data.py not found: $EXTRACT_SCRIPT"; exit 1; }
  if [[ -z "$EXTRACT_SIM_DATA_CSV" ]]; then
    EXTRACT_SIM_DATA_CSV="$WORK_DIR_ABS/extracted_sim_data_summary.csv"
  fi
fi

# Match both styles seen in runs:
#   <capsid>_F<fold>id<fold_id>_R<resolution>_<timestamp>
#   <capsid>_F<fold>id<fold_id>_R<resolution>_S<seed>_<timestamp>
DIR_RE='^([A-Za-z0-9]+)_F([0-9]+)id([0-9]+)_R([0-9]+(\.[0-9]+)?)(_S[0-9]+)?_([0-9]{8}T[0-9]{6}Z)$'

summary_header=$'run_dir\tcapsid\tfold\tfold_id\tresolution\tstatus\tresult'
rows=()
malformed=0
ignored_batch=0

json_displacement_magnitude_max() {
  local json_path="$1"
  python3 - "$json_path" <<'PY'
import json
import math
import sys
from pathlib import Path

try:
    data = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    value = data["results"]["Displacement"]["magnitude"]["max"]
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("Displacement magnitude max is not a finite number")
except Exception as exc:
    print(f"error reading {sys.argv[1]}: {exc}", file=sys.stderr)
    sys.exit(1)
print(value)
PY
}

write_extracted_sim_data_csv() {
  local extracted=0
  local reused=0
  local failed=0
  local eligible=0

  {
    echo "run_dir,capsid,fold,fold_id,resolution,Displacement:magnitude:max"
    if (( ${#rows[@]} > 0 )); then
      while IFS=$'\t' read -r run_path run_dir capsid fold fold_id resolution status result; do
        if [[ "$status" != "COMPLETED" || "$result" != "COMPLETED" ]]; then
          continue
        fi
        ((eligible+=1)) || true
        local json_path="$run_path/extracted_sim_data.json"
        local value="NaN"
        local extraction_failed=0
        if [[ -f "$json_path" && FORCE_EXTRACT_SIM_DATA -eq 0 ]]; then
          ((reused+=1)) || true
        else
          if (cd "$run_path" && python3 "$EXTRACT_SCRIPT" octreemesh.post.res --format json --output extracted_sim_data.json); then
            ((extracted+=1)) || true
          else
            echo "WARN failed to extract simulation data: $run_dir" >&2
            ((failed+=1)) || true
            extraction_failed=1
          fi
        fi
        if (( extraction_failed == 0 )); then
          if ! value="$(json_displacement_magnitude_max "$json_path")"; then
            echo "WARN failed to read Displacement:magnitude:max: $json_path" >&2
            value="NaN"
            ((failed+=1)) || true
          fi
        fi
        printf "%s,%s,%s,%s,%s,%s\n" "$run_dir" "$capsid" "$fold" "$fold_id" "$resolution" "$value"
      done < <(
        printf '%s\n' "${rows[@]}" \
          | awk -F'\t' '{printf "%s\t%d\t%d\t%.15g\t%s\n", $3, $4, $5, $6 + 0, $0}' \
          | sort -s -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4g \
          | cut -f5-
      )
    fi
  } > "$EXTRACT_SIM_DATA_CSV"

  echo
  echo "# Extract simulation data"
  echo "eligible_completed_runs: $eligible"
  echo "extracted_new_json: $extracted"
  echo "reused_existing_json: $reused"
  echo "failed_extractions_or_reads: $failed"
  echo "Wrote extracted simulation CSV: $EXTRACT_SIM_DATA_CSV"

  if (( STRICT == 1 && failed > 0 )); then
    exit 1
  fi
}

while IFS= read -r run_path; do
  run_name="$(basename "$run_path")"
  if [[ "$run_name" =~ $DIR_RE ]]; then
    capsid="${BASH_REMATCH[1]}"
    fold="${BASH_REMATCH[2]}"
    fold_id="${BASH_REMATCH[3]}"
    resolution="${BASH_REMATCH[4]}"
    status="MISSING_LOG"
    log_path="$run_path/batch_run.log"
    if [[ -f "$log_path" ]]; then
      if grep -q "Selected steps completed successfully!" "$log_path"; then
        status="COMPLETED"
      elif grep -q "Step timing summary" "$log_path"; then
        status="PARTIAL"
      else
        status="FAILED_OR_UNKNOWN"
      fi
    fi

    result="MISSING_RESULT"
    res_file="$run_path/octreemesh.post.res"
    msh_file="$run_path/octreemesh.post.msh"
    solver_out_file="$run_path/octreemesh_solver.out"
    if [[ -s "$res_file" && -s "$msh_file" && -s "$solver_out_file" ]]; then
      result="COMPLETED"
    fi

    rows+=("${run_path}"$'\t'"${run_name}"$'\t'"${capsid}"$'\t'"${fold}"$'\t'"${fold_id}"$'\t'"${resolution}"$'\t'"${status}"$'\t'"${result}")
  else
    if [[ "$run_name" == batch_* ]]; then
      ((ignored_batch+=1)) || true
      continue
    fi
    echo "WARN malformed run directory name: $run_name"
    ((malformed+=1))
  fi
done < <(find "$WORK_DIR_ABS" -mindepth 1 -maxdepth 1 -type d | sort)

echo "# Stage 1 batch directory analysis"
echo "work_dir: $WORK_DIR_ABS"
echo "total_run_dirs: $(find "$WORK_DIR_ABS" -mindepth 1 -maxdepth 1 -type d | wc -l | awk '{print $1}')"
echo "matched: ${#rows[@]}"
echo "malformed: $malformed"
echo "ignored_batch: $ignored_batch"
echo

echo "$summary_header"
if (( ${#rows[@]} > 0 )); then
  printf '%s\n' "${rows[@]}" \
    | awk -F'\t' '{printf "%s\t%d\t%d\t%.15g\t%s\n", $3, $4, $5, $6 + 0, $0}' \
    | sort -s -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4g \
    | cut -f6-
fi

echo
echo "# Unique capsids"
if (( ${#rows[@]} > 0 )); then
  printf '%s\n' "${rows[@]}" | awk -F'\t' '{print $3}' | sort -u | sed 's/^/- /'
else
  echo "- (none)"
fi

echo "# Unique folds"
if (( ${#rows[@]} > 0 )); then
  printf '%s\n' "${rows[@]}" | awk -F'\t' '{print "F"$4"id"$5}' | sort -u | sed 's/^/- /'
else
  echo "- (none)"
fi

echo "# Unique resolutions"
if (( ${#rows[@]} > 0 )); then
  printf '%s\n' "${rows[@]}" | awk -F'\t' '{print $6}' | sort -n -u | sed 's/^/- /'
else
  echo "- (none)"
fi

if [[ -n "$TSV_OUT" ]]; then
  {
    echo "$summary_header"
    if (( ${#rows[@]} > 0 )); then
      printf '%s\n' "${rows[@]}" \
        | awk -F'\t' '{printf "%s\t%d\t%d\t%.15g\t%s\n", $3, $4, $5, $6 + 0, $0}' \
        | sort -s -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4g \
        | cut -f6-
    fi
  } > "$TSV_OUT"
  echo "Wrote TSV report: $TSV_OUT"
fi

if [[ -n "$CSV_OUT" ]]; then
  {
    echo "run_dir,capsid,fold,fold_id,resolution,status,result"
    if (( ${#rows[@]} > 0 )); then
      printf '%s\n' "${rows[@]}" \
        | awk -F'\t' '{printf "%s\t%d\t%d\t%.15g\t%s\n", $3, $4, $5, $6 + 0, $0}' \
        | sort -s -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4g \
        | cut -f6- \
        | awk -F'\t' '{printf "%s,%s,%s,%s,%s,%s,%s\n", $1,$2,$3,$4,$5,$6,$7}'
    fi
  } > "$CSV_OUT"
  echo "Wrote CSV report: $CSV_OUT"
fi

if (( EXTRACT_SIM_DATA == 1 )); then
  write_extracted_sim_data_csv
fi

if (( STRICT == 1 && malformed > 0 )); then
  exit 1
fi
