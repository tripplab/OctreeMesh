#!/bin/bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $0 --work_dir DIR [--tsv OUT.tsv] [--csv OUT.csv] [--strict]

Stage 1: Analyze batch run directory names and report discovered capsids/folds/resolutions.
Expected run directory format:
  <capsid>_F<fold>id<fold_id>_R<resolution>_<time_stamp>
Example:
  3J4U_F5id0_R16.00_S0_20260511T041629Z

Options:
  --work_dir DIR   Directory containing run directories to scan (required)
  --tsv PATH       Optional TSV report output path
  --csv PATH       Optional CSV report output path
  --strict         Exit non-zero if any malformed entries are encountered
  -h, --help       Show this help
USAGE
}

for req in find awk sed sort uniq; do
  command -v "$req" >/dev/null 2>&1 || { echo "Missing required tool: $req"; exit 1; }
done

WORK_DIR=""
TSV_OUT=""
CSV_OUT=""
STRICT=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work_dir) WORK_DIR="${2:-}"; shift 2 ;;
    --tsv) TSV_OUT="${2:-}"; shift 2 ;;
    --csv) CSV_OUT="${2:-}"; shift 2 ;;
    --strict) STRICT=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

[[ -n "$WORK_DIR" ]] || { echo "--work_dir is required"; usage; exit 1; }
[[ -d "$WORK_DIR" ]] || { echo "--work_dir not found: $WORK_DIR"; exit 1; }

WORK_DIR_ABS="$(cd "$WORK_DIR" && pwd)"

# Match both styles seen in runs:
#   <capsid>_F<fold>id<fold_id>_R<resolution>_<timestamp>
#   <capsid>_F<fold>id<fold_id>_R<resolution>_S<seed>_<timestamp>
DIR_RE='^([A-Za-z0-9]+)_F([0-9]+)id([0-9]+)_R([0-9]+(\.[0-9]+)?)(_S[0-9]+)?_([0-9]{8}T[0-9]{6}Z)$'

summary_header=$'run_dir\tcapsid\tfold\tfold_id\tresolution\tstatus\tresult'
rows=()
malformed=0

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

    rows+=("${run_name}"$'\t'"${capsid}"$'\t'"${fold}"$'\t'"${fold_id}"$'\t'"${resolution}"$'\t'"${status}"$'\t'"${result}")
  else
    echo "WARN malformed run directory name: $run_name"
    ((malformed+=1))
  fi
done < <(find "$WORK_DIR_ABS" -mindepth 1 -maxdepth 1 -type d | sort)

echo "# Stage 1 batch directory analysis"
echo "work_dir: $WORK_DIR_ABS"
echo "total_run_dirs: $(find "$WORK_DIR_ABS" -mindepth 1 -maxdepth 1 -type d | wc -l | awk '{print $1}')"
echo "matched: ${#rows[@]}"
echo "malformed: $malformed"
echo

echo "$summary_header"
if (( ${#rows[@]} > 0 )); then
  printf '%s\n' "${rows[@]}" \
    | awk -F'\t' '{printf "%s\t%d\t%d\t%.15g\t%s\n", $2, $3, $4, $5 + 0, $0}' \
    | sort -s -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4g \
    | cut -f5-
fi

echo
echo "# Unique capsids"
if (( ${#rows[@]} > 0 )); then
  printf '%s\n' "${rows[@]}" | awk -F'\t' '{print $2}' | sort -u | sed 's/^/- /'
else
  echo "- (none)"
fi

echo "# Unique folds"
if (( ${#rows[@]} > 0 )); then
  printf '%s\n' "${rows[@]}" | awk -F'\t' '{print "F"$3"id"$4}' | sort -u | sed 's/^/- /'
else
  echo "- (none)"
fi

echo "# Unique resolutions"
if (( ${#rows[@]} > 0 )); then
  printf '%s\n' "${rows[@]}" | awk -F'\t' '{print $5}' | sort -n -u | sed 's/^/- /'
else
  echo "- (none)"
fi

if [[ -n "$TSV_OUT" ]]; then
  {
    echo "$summary_header"
    if (( ${#rows[@]} > 0 )); then
      printf '%s\n' "${rows[@]}" \
        | awk -F'\t' '{printf "%s\t%d\t%d\t%.15g\t%s\n", $2, $3, $4, $5 + 0, $0}' \
        | sort -s -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4g \
        | cut -f5-
    fi
  } > "$TSV_OUT"
  echo "Wrote TSV report: $TSV_OUT"
fi

if [[ -n "$CSV_OUT" ]]; then
  {
    echo "run_dir,capsid,fold,fold_id,resolution,status,result"
    if (( ${#rows[@]} > 0 )); then
      printf '%s\n' "${rows[@]}" \
        | awk -F'\t' '{printf "%s\t%d\t%d\t%.15g\t%s\n", $2, $3, $4, $5 + 0, $0}' \
        | sort -s -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4g \
        | cut -f5- \
        | awk -F'\t' '{printf "%s,%s,%s,%s,%s,%s,%s\n", $1,$2,$3,$4,$5,$6,$7}'
    fi
  } > "$CSV_OUT"
  echo "Wrote CSV report: $CSV_OUT"
fi

if (( STRICT == 1 && malformed > 0 )); then
  exit 1
fi
