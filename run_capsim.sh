#!/bin/bash
# OctreeMesh Pipeline by trippm@bmd [tripplab.com]
# v160326 - Improved version with config file and step selection

# Display usage information
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -c, --config FILE    Configuration file (default: capsim_config.sh)"
    echo "  -s, --steps STEPS     Steps to run (comma-separated or range)"
    echo "                        Available steps:"
    echo "                        1: extract_ATOM"
    echo "                        2: octree_mesh"
    echo "                        3: meshsolver"
    echo "                        4: mesh2pdb"
    echo "                        5: rotate_back"
    echo "                        Examples: 1-5 (all steps), 2,3,4 (specific steps)"
    echo "  --shear               Enable shear force simulation (adds rotation step after mesh)"
    echo "  -t, --threads N       Number of threads for FEM solver (overrides config)"
    echo "  -l, --list            List available fold configurations"
    echo "  -h, --help            Display this help message"
    exit 1
}

# List available fold configurations
list_folds() {
    echo "Available fold configurations:"
    echo ""
    echo "Fold Type 2 (two-fold axis):"
    echo "  Index 0: (0.000, 0.000, 125.100)"
    echo "  Index 1: (23.892, 38.658, 125.100)"
    echo ""
    echo "Fold Type 3 (three-fold axis):"
    echo "  Index 0: (47.784, 0.000, 125.100)"
    echo "  Index 1: (-47.784, 0.000, 125.100)"
    echo ""
    echo "Fold Type 5 (five-fold axis):"
    echo "  Index 0: (0.000, 77.316, 125.100)"
    echo ""
    echo "Custom: Set FOLD_TYPE=custom and define CUSTOM_FX, CUSTOM_FY, CUSTOM_FZ"
}

# Parse step specification
parse_steps() {
    local stepspec=$1
    STEPS_TO_RUN=()
    
    if [[ $stepspec == *-* ]]; then
        # Range format (e.g., 1-5)
        start=${stepspec%-*}
        end=${stepspec#*-}
        for ((i=start; i<=end; i++)); do
            STEPS_TO_RUN+=($i)
        done
    else
        # Comma-separated format
        IFS=',' read -ra STEPS <<< "$stepspec"
        for step in "${STEPS[@]}"; do
            STEPS_TO_RUN+=($step)
        done
    fi
}

# Check if a step should be run
should_run() {
    local step=$1
    for s in "${STEPS_TO_RUN[@]}"; do
        if [[ $s -eq $step ]]; then
            return 0
        fi
    done
    return 1
}

# Get step name
get_step_name() {
    case $1 in
        1) echo "extract_ATOM" ;;
        2) echo "octree_mesh" ;;
        3) echo "meshsolver" ;;
        4) echo "mesh2pdb" ;;
        5) echo "rotate_back" ;;
        s) echo "shear_rotate" ;;
        *) echo "unknown" ;;
    esac
}

# Get epoch seconds with sub-second precision when available
now_epoch() {
    date +%s.%N
}

# Compute elapsed seconds between two epoch values
elapsed_seconds() {
    local start_time=$1
    local end_time=$2
    awk -v s="$start_time" -v e="$end_time" 'BEGIN { printf "%.3f", (e - s) }'
}

# Append timing metadata for cross-run comparisons
record_step_timing() {
    local step_index=$1
    local step_name=$2
    local status=$3
    local elapsed=$4
    local note=${5:-""}

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$RUN_ID" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$PDB" "$Res" "$FOLD_TYPE" "$FOLD_INDEX" \
        "$step_index" "$step_name" "$status" "$elapsed" "$SOLVER_THREADS" "$STEPS_SPEC" "$SHEAR_MODE" "$note" >> "$TIMING_FILE"

    STEP_SUMMARY+=("$step_index|$step_name|$status|$elapsed")
}

# Set fold parameters based on type and index
set_fold_parameters() {
    case $FOLD_TYPE in
        2)
            case $FOLD_INDEX in
                0)
                    Fold=2; inx=0; Fx=00.000; Fy=00.000; Fz=125.100
                    FOLD_DESC="2-fold 0"
                    ;;
                1)
                    Fold=2; inx=1; Fx=23.892; Fy=38.658; Fz=125.100
                    FOLD_DESC="2-fold 1"
                    ;;
                *)
                    echo "Error: Invalid index $FOLD_INDEX for fold type 2"
                    echo "Valid indices: 0, 1"
                    exit 1
                    ;;
            esac
            ;;
        3)
            case $FOLD_INDEX in
                0)
                    Fold=3; inx=0; Fx=47.784; Fy=00.000; Fz=125.100
                    FOLD_DESC="3-fold 0"
                    ;;
                1)
                    Fold=3; inx=1; Fx=-47.784; Fy=00.000; Fz=125.100
                    FOLD_DESC="3-fold 1"
                    ;;
                *)
                    echo "Error: Invalid index $FOLD_INDEX for fold type 3"
                    echo "Valid indices: 0, 1"
                    exit 1
                    ;;
            esac
            ;;
        5)
            case $FOLD_INDEX in
                0)
                    Fold=5; inx=0; Fx=00.000; Fy=77.316; Fz=125.100
                    FOLD_DESC="5-fold"
                    ;;
                *)
                    echo "Error: Invalid index $FOLD_INDEX for fold type 5"
                    echo "Valid indices: 0"
                    exit 1
                    ;;
            esac
            ;;
        custom)
            if [[ -z "$CUSTOM_FX" ]] || [[ -z "$CUSTOM_FY" ]] || [[ -z "$CUSTOM_FZ" ]]; then
                echo "Error: Custom fold selected but CUSTOM_FX, CUSTOM_FY, CUSTOM_FZ not defined"
                exit 1
            fi
            Fold=0; inx=0; Fx=$CUSTOM_FX; Fy=$CUSTOM_FY; Fz=$CUSTOM_FZ
            FOLD_DESC="custom ($Fx, $Fy, $Fz)"
            ;;
        *)
            echo "Error: Invalid FOLD_TYPE '$FOLD_TYPE'"
            echo "Valid options: 2, 3, 5, custom"
            list_folds
            exit 1
            ;;
    esac
}

# Default values
CONFIG_FILE="capsim_config.sh"
STEPS_SPEC="1-5"  # Default: run all steps
SHEAR_MODE=0
CLI_THREADS=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -s|--steps)
            STEPS_SPEC="$2"
            shift 2
            ;;
        --shear)
            SHEAR_MODE=1
            shift
            ;;
        -t|--threads)
            CLI_THREADS="$2"
            shift 2
            ;;
        -l|--list)
            list_folds
            exit 0
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Load configuration
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: Configuration file $CONFIG_FILE not found!"
    echo "Please create it or specify with -c option."
    exit 1
fi

source "$CONFIG_FILE"

# Override threads if provided via CLI
if [[ -n "$CLI_THREADS" ]]; then
    SOLVER_THREADS=$CLI_THREADS
fi

# Export threads for FEM solver
export SOLVER_THREADS

# Set fold parameters based on config
set_fold_parameters

parse_steps "$STEPS_SPEC"

# Validate steps (always 1-5 regardless of shear mode)
for step in "${STEPS_TO_RUN[@]}"; do
    if [[ $step -lt 1 ]] || [[ $step -gt 5 ]]; then
        echo "Error: Step $step is out of range (1-5)"
        exit 1
    fi
done

# Calculate total steps for display
TOTAL_STEPS=5
if [[ $SHEAR_MODE -eq 1 ]]; then
    TOTAL_STEPS=6
fi

# Display run information
echo "=================================================="
echo "OctreeMesh Pipeline"
echo "=================================================="
echo "Configuration: $CONFIG_FILE"
echo "Steps to run: ${STEPS_TO_RUN[@]}"
echo "Shear mode: $([ $SHEAR_MODE -eq 1 ] && echo "Enabled" || echo "Disabled")"
echo "FEM threads: $SOLVER_THREADS"
echo "PDB: $PDB, Resolution: $Res Å"
echo "Indentation: $FOLD_DESC"
if [[ $FOLD_TYPE != "custom" ]]; then
    echo "  (Type $FOLD_TYPE, Index $FOLD_INDEX)"
fi
echo "  Vector: ($Fx, $Fy, $Fz)"
echo "=================================================="
echo ""

# Setup filenames
name="octreemesh"
label="F${Fold}id${inx}"
resn="${PDB}_${label}_${Res}"
CAPR="capsid_rotated.pdb"

GEOM="${name}.geometry.dat"
PROB="${name}.problem.dat"
SOLV="${name}.solver.dat"
MSH1="capside_atoms.msh"
WARN="warnings.log"

OUTF="${name}_solver.out"
POST="${name}.post.res"
MSH2="${name}.post.msh"
PDB_results="${resn}.pdb"
PDB_back="${resn}_back.pdb"

# Create checkpoints directory
mkdir -p .checkpoints

# Timing metadata output (kept across runs for comparisons)
TIMING_FILE=".checkpoints/timings.ts"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
if [[ ! -f "$TIMING_FILE" ]]; then
    printf "run_id\tutc_timestamp\tpdb\tres\tfold_type\tfold_index\tstep_index\tstep_name\tstatus\telapsed_sec\tsolver_threads\tsteps_spec\tshear_mode\tnote\n" > "$TIMING_FILE"
fi
STEP_SUMMARY=()
PIPELINE_START="$(now_epoch)"

# Initialize step counter for display
current_display_step=1

# Step 1: Extract ATOM records
if should_run 1; then
    echo "[Step $current_display_step/$TOTAL_STEPS] extract_ATOM..."
    if [[ ! -f ".checkpoints/step1_done" ]] || [[ ! -f "${VDB}_ATOMS.vdb" ]]; then
        step_start="$(now_epoch)"
        ${BIN}/extract_ATOM "${VDB}.vdb" "${VDB}_ATOMS.vdb"
        step_rc=$?
        step_elapsed="$(elapsed_seconds "$step_start" "$(now_epoch)")"
        if [[ $step_rc -eq 0 ]]; then
            touch .checkpoints/step1_done
            echo "  Time ${step_elapsed}s"
            echo "  ✓ Step $current_display_step complete"
            record_step_timing "1" "extract_ATOM" "done" "$step_elapsed"
        else
            record_step_timing "1" "extract_ATOM" "failed" "$step_elapsed"
            echo "  ✗ Step $current_display_step failed!"
            exit 1
        fi
    else
        echo "  ↻ Step $current_display_step already completed (checkpoint found)"
        record_step_timing "1" "extract_ATOM" "skipped" "0.000" "checkpoint"
    fi
    ((current_display_step++))
fi

# Step 2: Generate mesh and config files
if should_run 2; then
    echo "[Step $current_display_step/$TOTAL_STEPS] octree_mesh..."
    if [[ ! -f ".checkpoints/step2_done" ]] || [[ ! -f "$GEOM" ]]; then
        inp="${VDB}_ATOMS.vdb ${T} ${VDW} ${Res} ${Fold} ${inx} ${Fx} ${Fy} ${Fz} ${PDB} ${cone} ${load_ele} ${Young}"
        step_start="$(now_epoch)"
        ${BIN}/octree_mesh ${inp}
        step_rc=$?
        step_elapsed="$(elapsed_seconds "$step_start" "$(now_epoch)")"
        if [[ $step_rc -eq 0 ]]; then
            touch .checkpoints/step2_done
            echo "  Time ${step_elapsed}s"
            echo "  ✓ Step $current_display_step complete"
            record_step_timing "2" "octree_mesh" "done" "$step_elapsed"
        else
            record_step_timing "2" "octree_mesh" "failed" "$step_elapsed"
            echo "  ✗ Step $current_display_step failed!"
            exit 1
        fi
    else
        echo "  ↻ Step $current_display_step already completed (checkpoint found)"
        record_step_timing "2" "octree_mesh" "skipped" "0.000" "checkpoint"
    fi
    ((current_display_step++))
fi

# Shear rotation step (only if in shear mode and step 2 was just run or we have checkpoint)
if [[ $SHEAR_MODE -eq 1 ]]; then
    # Only run shear if we're going to run solver (step 3) or if explicitly in steps?
    # For simplicity, we run it if shear mode is enabled and step 2 is complete
    if [[ -f ".checkpoints/step2_done" ]] || should_run 2; then
        echo "[Step $current_display_step/$TOTAL_STEPS] shear_rotate..."
        if [[ ! -f ".checkpoints/shear_done" ]]; then
            echo "  Applying shear rotation for shearing force simulation"
            step_start="$(now_epoch)"
            ${BIN}/shear_rotate "${GEOM}" "${GEOM}"
            step_rc=$?
            step_elapsed="$(elapsed_seconds "$step_start" "$(now_epoch)")"
            if [[ $step_rc -eq 0 ]]; then
                touch .checkpoints/shear_done
                echo "  Time ${step_elapsed}s"
                echo "  ✓ Step $current_display_step complete"
                record_step_timing "s" "shear_rotate" "done" "$step_elapsed"
            else
                record_step_timing "s" "shear_rotate" "failed" "$step_elapsed"
                echo "  ✗ Step $current_display_step failed!"
                exit 1
            fi
        else
            echo "  ↻ Step $current_display_step already completed (checkpoint found)"
            record_step_timing "s" "shear_rotate" "skipped" "0.000" "checkpoint"
        fi
        ((current_display_step++))
    fi
fi

# Step 3: Run simulation
if should_run 3; then
    echo "[Step $current_display_step/$TOTAL_STEPS] meshsolver (using $SOLVER_THREADS threads)..."
    if [[ ! -f ".checkpoints/step3_done" ]] || [[ ! -f "$OUTF" ]]; then
        step_start="$(now_epoch)"
        ${BIN}/meshsolver "${name}" "${log_lev}" > "${OUTF}"
        step_rc=$?
        step_elapsed="$(elapsed_seconds "$step_start" "$(now_epoch)")"
        if [[ $step_rc -eq 0 ]]; then
            touch .checkpoints/step3_done
            echo "  Time ${step_elapsed}s"
            echo "  ✓ Step $current_display_step complete"
            record_step_timing "3" "meshsolver" "done" "$step_elapsed"
        else
            record_step_timing "3" "meshsolver" "failed" "$step_elapsed"
            echo "  ✗ Step $current_display_step failed!"
            exit 1
        fi
    else
        echo "  ↻ Step $current_display_step already completed (checkpoint found)"
        record_step_timing "3" "meshsolver" "skipped" "0.000" "checkpoint"
    fi
    ((current_display_step++))
fi

# Step 4: Map results to PDB
if should_run 4; then
    echo "[Step $current_display_step/$TOTAL_STEPS] mesh2pdb..."
    if [[ ! -f ".checkpoints/step4_done" ]] || [[ ! -f "$PDB_results" ]]; then
        step_start="$(now_epoch)"
        ${BIN}/mesh2pdb "${CAPR}" "${MSH2}" "${POST}" "${PDB_results}" "${ref_lev}"
        step_rc=$?
        rm -f "${CAPR}"  # Clean up rotated file
        step_elapsed="$(elapsed_seconds "$step_start" "$(now_epoch)")"
        if [[ $step_rc -eq 0 ]]; then
            touch .checkpoints/step4_done
            echo "  Time ${step_elapsed}s"
            echo "  ✓ Step $current_display_step complete"
            record_step_timing "4" "mesh2pdb" "done" "$step_elapsed"
        else
            record_step_timing "4" "mesh2pdb" "failed" "$step_elapsed"
            echo "  ✗ Step $current_display_step failed!"
            exit 1
        fi
    else
        echo "  ↻ Step $current_display_step already completed (checkpoint found)"
        record_step_timing "4" "mesh2pdb" "skipped" "0.000" "checkpoint"
    fi
    ((current_display_step++))
fi

# Step 5: Rotate back to original orientation
if should_run 5; then
    echo "[Step $current_display_step/$TOTAL_STEPS] rotate_back..."
    if [[ ! -f ".checkpoints/step5_done" ]] || [[ ! -f "$PDB_back" ]]; then
        step_start="$(now_epoch)"
        ${BIN}/apply-matrix.awk ${PDB_results} rotate_Z2F.mtx > ${PDB_back}
        step_rc=$?
        step_elapsed="$(elapsed_seconds "$step_start" "$(now_epoch)")"
        if [[ $step_rc -eq 0 ]]; then
            touch .checkpoints/step5_done
            echo "  Time ${step_elapsed}s"
            echo "  ✓ Step $current_display_step complete"
            record_step_timing "5" "rotate_back" "done" "$step_elapsed"
        else
            record_step_timing "5" "rotate_back" "failed" "$step_elapsed"
            echo "  ✗ Step $current_display_step failed!"
            exit 1
        fi
    else
        echo "  ↻ Step $current_display_step already completed (checkpoint found)"
        record_step_timing "5" "rotate_back" "skipped" "0.000" "checkpoint"
    fi
    ((current_display_step++))
fi

# Check if all requested steps completed
all_completed=true

PIPELINE_ELAPSED="$(elapsed_seconds "$PIPELINE_START" "$(now_epoch)")"
echo ""
echo "Step timing summary (run: $RUN_ID)"
echo "--------------------------------------------------"
printf "%-4s %-15s %-9s %s\n" "Step" "Name" "Status" "Elapsed"
for step_row in "${STEP_SUMMARY[@]}"; do
    IFS='|' read -r summary_step summary_name summary_status summary_elapsed <<< "$step_row"
    printf "%-4s %-15s %-9s %ss\n" "$summary_step" "$summary_name" "$summary_status" "$summary_elapsed"
done
echo "--------------------------------------------------"
echo "Total pipeline wall time: ${PIPELINE_ELAPSED}s"
echo "Timing metadata saved to: $TIMING_FILE"
echo ""

for step in "${STEPS_TO_RUN[@]}"; do
    if [[ ! -f ".checkpoints/step${step}_done" ]]; then
        all_completed=false
        break
    fi
done

# Also check shear checkpoint if in shear mode
if [[ $SHEAR_MODE -eq 1 ]] && [[ ! -f ".checkpoints/shear_done" ]]; then
    all_completed=false
fi

if $all_completed; then
    echo ""
    echo "=================================================="
    echo "Selected steps completed successfully!"
    echo "=================================================="
    
    # Show output files if final steps were run
    if should_run 5 && [[ -f "$PDB_back" ]]; then
        echo "Results in:"
        echo "  - $OUTF (solver output)"
        echo "  - ${resn}.post.res (GID results)"
        echo "  - ${resn}.post.msh (GID mesh)"
        echo "  - $PDB_back (PDB with results for VMD)"
    fi
    
    # Ask about checkpoint cleanup if all steps in full pipeline were run
    if [[ ${#STEPS_TO_RUN[@]} -eq 5 ]] && [[ $all_completed == true ]]; then
        echo ""
        read -p "Clean up checkpoints? (y/n): " -n 1 -r
        echo ""
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            find .checkpoints -mindepth 1 -not -name "$(basename "$TIMING_FILE")" -delete
            echo "Checkpoints cleaned up (timing history preserved in $TIMING_FILE)."
        fi
    fi
else
    echo ""
    echo "=================================================="
    echo "Partial completion. Checkpoints preserved."
    echo "To resume, run with: $0 --steps=${STEPS_SPEC} $([ $SHEAR_MODE -eq 1 ] && echo "--shear")"
    echo "=================================================="
fi
