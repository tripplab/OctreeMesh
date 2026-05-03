# OctreeMesh suite

OCTREEMESH PIPELINE v160326 by trippm@tripplab.com

Some documentation and examples in HOW_TO folder

# OctreeMesh Pipeline - User Guide

# Quick Start

## Installation / compile / run (practical)

git clone -b dev01 --single-branch https://github.com/tripplab/OctreeMesh.git

or

git pull origin dev01

Prereqs:

- make, g++

- gfortran (required by solver build flags)

- OpenMP-capable compiler/runtime (-fopenmp)

shell tools in script path: gunzip, awk

Compile everything

From repo root:

make

This builds each module and copies user-facing binaries into ./bin:

- extract_ATOM

- octree_mesh

- shear_rotate

- meshsolver

- mesh2pdb

Run:

Preferred method is editing capsim_config.sh and running 

chmod u+x run_capsim.sh

./run_capsim.sh



### 1. Configure your simulation

Edit `capsim_config.sh`:

```bash
# Basic settings
PDB=1CWP                 # Capsid identifier
VDB=1cwp_full            # VDB file name (without .vdb)
Res=16.00                # Mesh resolution in angstroms

# Select indentation point (fold type and index)
FOLD_TYPE=2              # Options: 2, 3, 5, or custom
FOLD_INDEX=0             # Index: 0 or 1 (where applicable)

# For custom indentation (if FOLD_TYPE=custom):
#CUSTOM_FX=10.0
#CUSTOM_FY=20.0
#CUSTOM_FZ=30.0

# Other parameters
T=3                      # T-number (1 or 3)
VDW=1                    # Use van der Waals radius (1=yes, 0=no)
cone=15.00               # Boundary condition cone angle
Young=0.020              # Young modulus
SOLVER_THREADS=1         # Number of threads for FEM solver

# Path to executables
BIN=/path/to/OctreeMesh/bin
```

### 2. Available indentation points

View all options:
```bash
./run_capsim.sh --list
```

**Preconfigured points:**
- **2-fold:** Index 0, Index 1
- **3-fold:** Index 0, Index 1
- **5-fold:** Index 0

### 3. Run the pipeline

```bash
# Run complete pipeline
./run_capsim.sh

# Run with shear force enabled
./run_capsim.sh --shear

# Run specific steps only
./run_capsim.sh --steps=2-4        # Steps 2,3,4
./run_capsim.sh --steps=1,3,5      # Steps 1,3,5 only

# Use more threads
./run_capsim.sh --threads=4

# Combine options
./run_capsim.sh --shear --threads=8 --steps=2-5

# Use a different config file
./run_capsim.sh --config=my_experiment.conf
```


## Batch execution with `run_capsim_batch.sh`

Use `run_capsim_batch.sh` to execute many simulations **sequentially** across combinations of PDB, fold token, and resolution.

### What it does
- Builds a full Cartesian product in this exact order: **PDB (outer) × Fold token (middle) × Res (inner)**.
- Runs one simulation at a time (strictly sequential).
- Always self-detaches with `nohup` so jobs continue if your SSH/session ends.
- Creates per-job configs under `runs/batch_<UTC_TS>_<PID>/tmp_configs/`.
- Writes per-job run output to `<RUN_DIR>/batch_run.log`.
- Writes batch-level summaries:
  - `runs/batch_<UTC_TS>_<PID>.tsv` (orchestration/status fields)
  - `runs/batch_<UTC_TS>_<PID>.csv` (requested metrics/timings)

### Required arguments
```bash
./run_capsim_batch.sh   --pdb 1cwp,3j4u   --folds 2_0,2_1,3_0,3_1,5_0   --res 1-16   --threads 30   --bin /path/to/OctreeMesh/bin   --vdb-dir /path/to/vdbs
```

### Optional flags
- `--strict-skips`: if any job is skipped because VDB is missing, final batch exit becomes non-zero.
- `--smoke`: validate/build matrix and summaries without executing simulations.

### Input rules and inference
- Allowed PDB IDs only: `1cwp`, `3j4u`, `3izg`, `4g93` (input case-insensitive).
- PDB is written uppercase into temp config (`1cwp -> 1CWP`).
- VDB is inferred as `${pdb_lower}_full`.
- Hardcoded Young map:
  - `1cwp -> 0.0193`
  - `3j4u -> 0.0096`
  - `3izg -> 0.0245`
  - `4g93 -> 0.0147`
- `--res` supports mixed syntax (example: `1-4,8,10-12`) and each value must be within `[1,16]`.
- `Res` is written as formatted float (`1.00`, `2.00`, ...).
- `--bin` must contain these entries: `extract_ATOM`, `octree_mesh`, `meshsolver`, `mesh2pdb`, and `apply-matrix.awk`.
- Step 5 is named `rotate_back`, but it is executed through `apply-matrix.awk` (there is no separate `rotate_back` executable requirement).
- `--threads` is required and must be integer `[1,30]`.
- Fold tokens must be from: `2_0,2_1,3_0,3_1,5_0`.

### VDB existence check (mandatory)
For each inferred VDB file `<vdb>.vdb`, lookup order is:
1. `--vdb-dir/<vdb>.vdb`
2. repo root `<repo>/<vdb>.vdb`

If not found, job is logged as `SKIPPED_MISSING_VDB` in both TSV and CSV.

### Detach behavior and monitoring
When invoked, the parent process immediately prints:
- spawned PID
- master log path
- exact `tail -f` command

Example:
```bash
./run_capsim_batch.sh --pdb 1cwp --folds 2_0 --res 16 --threads 8 --bin "$PWD/bin" --vdb-dir "$PWD"
# prints:
# Spawned batch PID: <pid>
# Master log: runs/batch_<UTC_TS>_<PID>.log
# Tail command: tail -f runs/batch_<UTC_TS>_<PID>.log
```

### Smoke mode example
```bash
./run_capsim_batch.sh   --pdb 1cwp,3j4u   --folds 2_0,5_0   --res 4-6   --threads 12   --bin /path/to/OctreeMesh/bin   --vdb-dir /path/to/vdbs   --smoke
```
Expected behavior:
- still detaches
- still validates inputs/BIN/VDB checks
- still writes TSV+CSV summaries
- **does not run** `run_capsim.sh`
- non-missing jobs are marked `SMOKE_OK`

### Status and exit semantics
- Normal successful run: `DONE`
- Failed run: `FAILED_STEP_1` ... `FAILED_STEP_5` (or `FAILED_STEP_NA` if step cannot be inferred)
- Missing VDB: `SKIPPED_MISSING_VDB`
- Smoke non-skipped: `SMOKE_OK`

Final batch exit code:
- non-zero if any `FAILED_*`
- additionally non-zero for skipped jobs only when `--strict-skips` is set

### CSV/TSV contents
- TSV columns:
  - `status, exit_code, runtime_sec, pdb, vdb, res, young, fold_type, fold_index, threads, run_dir`
- CSV columns:
  - `job_name, total_proteins, total_atoms, nodes, elements, mesh_volume, volume_loaded, octree_mesh_sec, meshsolver_sec, mesh2pdb_sec`

Notes:
- numeric parsing in CSV normalizes values like `1,234` -> `1234` when possible
- missing fields are written as `NA`
- timings are read from `.checkpoints/timings.ts` for `octree_mesh`, `meshsolver`, `mesh2pdb`

### 4. Pipeline Steps

| Step | Name | Description |
|------|------|-------------|
| 1 | `extract_ATOM` | Extract ATOM records from VDB file |
| 2 | `octree_mesh` | Generate mesh and configuration files |
| 3 | `meshsolver` | Run FEM simulation |
| 4 | `mesh2pdb` | Map results back to PDB structure |
| 5 | `rotate_back` | Rotate back to original orientation |

*With `--shear`: Adds shear rotation step between steps 2 and 3*

### 5. Output Files

After successful completion:
- `*.solver.out` - Solver output log
- `*.post.res` - Results for GID visualization
- `*.post.msh` - Mesh for GID visualization
- `*_back.pdb` - PDB file with results (for VMD)

### 6. Resuming Interrupted Runs

The pipeline creates checkpoints automatically. To resume:

```bash
# Resume from where it left off
./run_capsim.sh --steps=3-5

# Or with same options as original run
./run_capsim.sh --shear --steps=3-5
```

### 7. Command Reference

| Option | Description |
|--------|-------------|
| `-c, --config FILE` | Use custom config file |
| `-s, --steps STEPS` | Steps to run (e.g., 1-5, 2,3,4) |
| `--shear` | Enable shear force simulation |
| `-t, --threads N` | Set number of FEM threads |
| `-l, --list` | List available indentation points |
| `-h, --help` | Show help message |

## Example Workflows

**Basic indentation test:**
```bash
# Edit config: FOLD_TYPE=2, FOLD_INDEX=0
./run_capsim.sh
```

**Shear force simulation with high resolution:**
```bash
# Edit config: Res=10.0
./run_capsim.sh --shear --threads=8
```

**Quick mesh generation only:**
```bash
./run_capsim.sh --steps=1-2
```

**Post-process existing results:**
```bash
./run_capsim.sh --steps=4-5
```

## Troubleshooting

- **"Command not found"**: Check `BIN` path in config file
- **"Invalid index"**: Run `./run_capsim.sh --list` to see valid options
- **Checkpoints not working**: Check write permissions in current directory


# High-level architecture

This repo is a pipeline for viral capsid meshing + FEM solve + back-interpolation to PDB:

extract_ATOM (from CLEANER) filters ATOM/amino-acid records from .vdb.

octree_mesh (from MESHER) builds octree hexahedral mesh and solver input files.

(optional) shear_rotate (from ROTATOR) rotates geometry for shear setups.

meshsolver (renamed Solid binary from SOLVER) runs FEM and writes .post.res.

mesh2pdb (from CREATES_PDB) interpolates element/node results back to atoms in PDB.

apply-matrix.awk rotates output PDB back to original orientation in the example script.

The top-level makefile wires exactly these binaries and names in ./bin.

# What each executable does
1) extract_ATOM (CLEANER/cleaner)

Purpose: reads an input .vdb/PDB-like text file, keeps only ATOM records that belong to amino-acid residues (ALA, ARG, ...), writes filtered output.

Inputs: <input_vdb> <output_vdb> (2 args after program name).

Outputs: filtered atom-only text file used by mesher.

How: parses fixed-width columns, checks record type + residue token, writes matching lines.

Parallelized? No (single-threaded C++ without OpenMP pragmas).

2) octree_mesh (MESHER/output/mesher)

Purpose: builds capsid volumetric mesh and FEM input data.

Inputs (13 args after exe): vdb file, T-number, solution type, resolution, fold/id alignment, vector (x,y,z), virus id, cone angle, load proportion variation, Young’s modulus scaling.

Outputs: mesh/problem/solver files consumed by FEM solver (and other intermediate files).

How: main flow is ReadVdbFile -> refinement/percentages -> local root atom assignment -> octree refine -> set loaded/fixed elements -> print FEM data files.

Parallelized? Compiled with -fopenmp, but I found no OpenMP pragmas in MESHER sources; practically appears single-threaded in this module.

3) shear_rotate (ROTATOR/rotator)

Purpose: rotates a FEM geometry mesh (90° around X) for shear mode.

Inputs: <input.geometry.dat> <output.geometry.dat>.

Outputs: rotated geometry file in same FEM format.

How: reads node/element counts + connectivity, applies rotation matrix to node coordinates, rewrites mesh.

Parallelized? No (single-threaded).

4) meshsolver (copied from SOLVER/gid/problemtypes/Solid.gid/Solid)

Purpose: the structural FEM solve stage used by OctreeMesh.

Inputs: <problem_prefix> [log_level] [log_file], where it expects:

<prefix>.solver.dat
<prefix>.geometry.dat
<prefix>.problem.dat

Outputs: <prefix>.post.res (GiD results), optionally <prefix>.post.mat.

How: loads solver settings + geometry + materials, assembles global system(s), runs selected solver, writes GiD results.

Parallelized? Yes, SOLVER is OpenMP-enabled (-fopenmp) and solver kernels use many #pragma omp parallel for loops.

5) mesh2pdb (CREATES_PDB/output/pdb)

Purpose: interpolates FEM results back onto atoms and writes a PDB with result fields in Occupancy/TempFactor columns.

Inputs: <input_vdb> <mesh.post.msh> <mesh.post.res> <output.pdb> <refinement_level>.

Outputs: result-annotated PDB.

How: translates nodal results, assigns atoms and hexahedra to local octree root, refines search tree, interpolates to atoms, saves.

Parallelized? Compiled with -fopenmp, but like MESHER, no explicit OpenMP pragmas found in this module.

Other executables under SOLVER (not core OctreeMesh pipeline, but present)

SOLVER/build/gcc/Makefile builds libraries, tools, examples, and problem binaries (Solid, Heat, etc.).

# Key tool binaries:

FEMSolver: pipe-based FE system assembly/solve service (FEMData/FEMResult named pipes).

EqnSolver: pipe-based sparse equation solver service (EqnData/EqnResult).

MatSolver: solves A*x=b from MAT-file (A and b in, x out).

Diff: subtracts two .post.res files and writes diff results.

# Inputs/outputs across the full OctreeMesh pipeline

Input structural file: ${VDB}.vdb

extract_ATOM output: ${VDB}_ATOMS.vdb

octree_mesh outputs include

${name}.geometry.dat

${name}.problem.dat

${name}.solver.dat

${name}.post.msh/intermediates

meshsolver output: ${name}.post.res + solver log redirection

mesh2pdb output: ${PDB_results}

apply-matrix.awk output: ${PDB_back}

Final renamed outputs for visualization: ${resn}.post.res, ${resn}.post.msh, ${resn}_back.pdb.



# Parallelization summary

Definitely parallelized: SOLVER core and tool internals (OpenMP loops in solver kernels; OpenMP sections for pipe setup in tool frontends).
Compiled with OpenMP but no explicit pragmas found in these modules: MESHER and CREATES_PDB.
Not parallelized: CLEANER and ROTATOR.

To compile this repo for multithreaded execution, you mainly need to ensure OpenMP is enabled at build time (it is, by default) and then set runtime thread counts for solver binaries.

1) Compile with OpenMP (already configured)

The relevant makefiles already include -fopenmp:

MESHER uses -fopenmp in both compile and link flags. 

CREATES_PDB uses -fopenmp in both compile and link flags. 

SOLVER build flags include -fopenmp for C++ builds in debug/release. 

So a normal build from repo root is enough:

make

Top-level make builds all subprojects and produces bin/* wrappers. 

2) Build parallelism vs runtime parallelism

There are two different “multithread” knobs:

Build parallelism (how fast compilation runs): 

in SOLVER/build/gcc/Makefile, PARALLEL=4 and sub-builds use make -j $(PARALLEL). 

Runtime parallelism (how many threads the solver uses while solving): 

controlled by environment variables like SOLVER_THREADS. 

3) Set runtime solver threads

For solver tools/scripts, increase threads by setting:

export SOLVER_THREADS=N

The provided runner defaults to 1, so change it there or export before execution. 

Example:

export SOLVER_THREADS=8

./bin/meshsolver octreemesh 2

4) Solver code actually uses OpenMP

Solver kernels contain OpenMP loops with num_threads(threads), e.g. conjugate gradient.

Practical recommendation:

Build normally with make (OpenMP already enabled). 

For actual speedup, tune SOLVER_THREADS to your core count (start with physical cores, then benchmark). 

If you want faster compilation, increase PARALLEL in SOLVER/build/gcc/Makefile or call make -j. 

## Verify OpenMP

# 1) How to tell if OpenMP is enabled at build time

For this repo specifically, it is configured in makefiles already:

MESHER compiles and links with -fopenmp.
CREATES_PDB compiles and links with -fopenmp.
SOLVER C++ flags include -fopenmp (debug and release).
Solver code also contains OpenMP pragmas (#pragma omp parallel for), so flags matter and are used.

Practical checks you can run

1) Show compile lines and confirm -fopenmp appears

make V=1

2) If a binary is built, confirm it links OpenMP runtime

ldd ./bin/meshsolver | rg "gomp|omp"

3) If unsure compiler supports OpenMP

echo '#include <omp.h>
#include <stdio.h>
int main(){printf("%d\n", omp_get_max_threads());}' > /tmp/omp_test.c
gcc -fopenmp /tmp/omp_test.c -o /tmp/omp_test && /tmp/omp_test

If (1) includes -fopenmp and (2) shows libgomp (or equivalent), you’re good.

# 2) If OpenMP is not available: micromamba env setup

## Using `environment.yml` with micromamba

From the repository root:

```bash
micromamba env create -f environment.yml
micromamba activate octreemesh
```

Build the project in that environment:

```bash
make
```

If you update dependencies in `environment.yml`, recreate the env:

```bash
micromamba env remove -n octreemesh
micromamba env create -f environment.yml
micromamba activate octreemesh
```

Quick sanity checks after activation:

```bash
which g++
which gfortran
g++ --version
```

To verify OpenMP support from the activated environment:

```bash
echo '#include <omp.h>\n#include <stdio.h>\nint main(){printf("%d\\n", omp_get_max_threads());}' > /tmp/omp_test.c
gcc -fopenmp /tmp/omp_test.c -o /tmp/omp_test && /tmp/omp_test
```


# 3) Extra tip: compile parallelism vs OpenMP runtime

These are different:

Compile parallelism (faster build): make -jN and solver submake PARALLEL=4.
Runtime thread count (OpenMP execution): set env vars like SOLVER_THREADS when running solver tools/scripts.


