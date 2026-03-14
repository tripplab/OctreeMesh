# OctreeMesh

OCTREEMESH v310122
Please take a look at documentation and examples in HOW_TO folder


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
Using HOW_TO/run.sh and top-level makefile naming:
Input structural file: ${VDB}.vdb
extract_ATOM output: ${VDB}_ATOMS.vdb
octree_mesh outputs include ${name}.geometry.dat, ${name}.problem.dat, ${name}.solver.dat, ${name}.post.msh/intermediates
meshsolver output: ${name}.post.res + solver log redirection
mesh2pdb output: ${PDB_results}
apply-matrix.awk output: ${PDB_back}
Final renamed outputs for visualization: ${resn}.post.res, ${resn}.post.msh, ${resn}_back.pdb.

# Installation / compile / run (practical)
Prereqs
make, g++
gfortran (required by solver build flags)
OpenMP-capable compiler/runtime (-fopenmp)
shell tools in script path: gunzip, awk

Compile everything
From repo root:
make
This builds each module and copies user-facing binaries into ./bin:
extract_ATOM
octree_mesh
shear_rotate
meshsolver
mesh2pdb

Run
Preferred method is editing and running HOW_TO/run.sh:
Set parameters (virus ID, fold vector, resolution, BIN path).
chmod u+x run.sh
./run.sh

The sample output in HOW_TO/README documents expected stage-by-stage progress and resulting files. 

# Parallelization summary
Definitely parallelized: SOLVER core and tool internals (OpenMP loops in solver kernels; OpenMP sections for pipe setup in tool frontends).
Compiled with OpenMP but no explicit pragmas found in these modules: MESHER and CREATES_PDB.
Not parallelized: CLEANER and ROTATOR.




