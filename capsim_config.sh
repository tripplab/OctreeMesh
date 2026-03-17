#!/bin/bash
# Configuration file for OctreeMesh pipeline

#--------------- User defined parameters -----------------------------------#
PDB=1CWP  ## capsid identifier
VDB=1cwp_full ## VDB file.vdb with capsid structure from VIPERdb (atoms)
Res=16.00  ## mesh resolution in ang [float]

## Select indentation point by type and index:
FOLD_TYPE=2    # Fold symmetry: 2, 3, or 5
FOLD_INDEX=0   # Index: 0 or 1 (where applicable)
# For custom indentation, set FOLD_TYPE=custom and define:
#CUSTOM_FX=10.000
#CUSTOM_FY=20.000
#CUSTOM_FZ=30.000

T=3 ## Capsid's T-number 1|3
VDW=1 ## use van der Waals radius for atoms
cone=15.00  ## angle in degrees of cone to set boundary conditions
load_ele=0.00  ## variation in the proportion of loaded elements
Young=0.020  ## Young modulus / 10,000
log_lev=2  ## 0: only fatal errors, 1: current process, 2: solver iterations
ref_lev=5  ## octree searching level for interpolation

## Number of threads for FEM solver
SOLVER_THREADS=1

## Set absolute path to executables
BIN=/path/to/OctreeMesh/bin
