##                                                  ##
## OctreeMesh Pipeline by trippm@bmd [tripplab.com] ##
##                 v260122                          ##
##               an example for CCMV                ##

## After the first time you can comment out this line
gunzip -v 1cwp*

#--------------- User defined parameters -----------------------------------#

VDB=1cwp_full ## VDB file with capsid structure from VIPERdb (atoms) [string]
CAPR=1cwp_3F-Z.pdb ##1cwp_full.vdb  ##1cwp_5F-Z.pdb ## PDB file with capsid structure with selected fold aligned to Z (atoms) [string]
Fold=3 ##2 ##5  ## fold to align with Z for nanoindentation (2|3|5) [int]
Res=16.00  ## mesh resolution in ang [float]
PDB=1CWP  ## capsid identifier [string]

T=3 ## 1|3 [int]
VDW=1 ## use van der Waals radius for atoms [int]
inx=0  ## Number of id to align with Y axis (0-11) [int]
cone=15.00  ## angle in degrees of cone to set boundary conditions [float]
load_ele=0.00  ## variation in the proportion of loaded elements [float]
Young=0.020  ## Young modulus / 10,000 [float]
log_lev=2  ## 0: only fatal errors, 1: current process, 2: solver iterations [int]
ref_lev=5  ##  to use in the octree searching of elemnts per atom to interpolate results written to PDB [int]

## Set absolute path to executables
BIN=/path/to/OctreeMesh/bin

#----- No need to modify below this line unless you know what you're doing -------#

name=${Fold}_T${T}_${PDB}

GEOM=${name}.geometry.dat
PROB=${name}.problem.dat
SOLV=${name}.solver.dat
MSH1=capside_mesh.msh
WARN=warnings.log

OUTF=${name}_solver.out
POST=${name}.post.res
MSH2=${name}.post.msh
PDB_results=${name}_${Res}.pdb

  ## Extract all ATOM records from the VIPERdb capsid structure file (virus.vdb)
  inp="${VDB}.vdb ${VDB}_ATOMS.vdb"
  echo "Runing cleanpdb $inp"
  ${BIN}/extract_ATOM ${inp}

  ## Generate the mesh and simulation config files
  inp="${VDB}_ATOMS.vdb ${T} ${VDW} ${Res} ${Fold} ${inx} ${PDB} ${cone} ${load_ele} ${Young}"
  echo "Running biomesh $inp"
  ${BIN}/octree_mesh ${inp}

  ## If simulating a shearing force, the mnesh has to be rotated 90 deg on X
  ##inp="${GEOM} ${GEOM}"
  ##echo "Running meshrotate $inp"
  ##${BIN}/shear_rotate ${inp}

  ## Run the simulation
  inp="${name} $log_lev"
  echo "Running meshsolver $inp"
  ${BIN}/meshsolver ${inp} > ${OUTF}

  ## Transpose the simulation results for the mesh into the capsid structure file (virus.pdb)
  ## suplied as the first argument. Please note that this structure has to be rotated so the
  ## correct axis-fold is aligned with Z to have the right orientation correspondence with the mesh
  inp="${CAPR} $MSH2 $POST ${PDB_results} ${ref_lev}"
  echo "Running mesh2pdb $inp"
  ${BIN}/mesh2pdb ${inp}

if [ -e $OUTF ]
then
	mv $POST ${name}_${Res}.post.res
	mv $MSH2 ${name}_${Res}.post.msh
 echo "Results in $OUTF and ${name}_${Res}.post.res to visualize with GID"
 echo "and mapped into ${name}_${Res}.pdb to visualize with VMD"
 echo "Job done"
else
	echo "Something went wrong"
fi




