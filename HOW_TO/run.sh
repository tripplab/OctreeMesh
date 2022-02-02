##                                                  ##
## OctreeMesh Pipeline by trippm@bmd [tripplab.com] ##
##                 v310122                          ##
##               an example for CCMV                ##

## After the first time you can comment out this line
##gunzip -v 1cwp*

#--------------- User defined parameters -----------------------------------#

PDB=1CWP  ## capsid identifier. As of this version, the code only works with 1CWP, 4G93, 3IZG, 3J4U [string]
VDB=1cwp_full ## VDB file.vdb with capsid structure from VIPERdb (atoms) [string]
Res=16.00  ## mesh resolution in ang [float]

## Set values for (Fold, inx, Fx, Fy, Fz) or uncomment one line for selected fold
## Fold is the fold-label for nanoindentation (2|3|5) [int]
## inx is the id-label of fold to align with Z axis [int]
## (Fx, Fy, Fz) are the vector components [float]

Fold=2; inx=0; Fx=00.000; Fy=00.000; Fz=125.100; ## 2-fold id 0
#Fold=2; inx=1; Fx=23.892; Fy=38.658; Fz=125.100; ## 2-fold id 1
#Fold=3; inx=0; Fx=47.784; Fy=00.000; Fz=125.100; ## 3-fold id 0
#Fold=3; inx=1; Fx=-47.784; Fy=00.000; Fz=125.100; ## 3-fold id 1
#Fold=5; inx=0; Fx=00.000; Fy=77.316; Fz=125.100; ## 5-fold id 0

#Fold=3; inx=2; Fx=-69; Fy=58; Fz=-41; ## random 2
#Fold=3; inx=3; Fx=23; Fy=-79; Fz=100; ## random 3

T=3 ## Capsid's T-number 1|3 [int]
VDW=1 ## use van der Waals radius for atoms [int]
cone=15.00  ## angle in degrees of cone to set boundary conditions [float]
load_ele=0.00  ## variation in the proportion of loaded elements [float]
Young=0.020  ## Young modulus / 10,000 [float]
log_lev=2  ## 0: only fatal errors, 1: current process, 2: solver iterations [int]
ref_lev=5  ##  to use in the octree searching of elemnts per atom to interpolate results written to PDB [int]

## Set absolute path to executables
BIN=/path/to/OctreeMesh/bin

#----- No need to modify below this line unless you know what you're doing -------#

##name=${label}_T${T}_${PDB}
name="octreemesh"
label=${Fold}Fid${inx}
resn=${PDB}_${label}_${Res}
CAPR="capsid_rotated.pdb"

GEOM=${name}.geometry.dat
PROB=${name}.problem.dat
SOLV=${name}.solver.dat
MSH1="capside_atoms.msh"
WARN="warnings.log"

OUTF=${name}_solver.out
POST=${name}.post.res
MSH2=${name}.post.msh
PDB_results=${resn}.pdb
PDB_back=${resn}_back.pdb

#----------------------- Pipeline Starts --------------------------------#

  ## Extract all aminoacid ATOM records from the VIPERdb capsid structure file (virus.vdb)
  inp="${VDB}.vdb ${VDB}_ATOMS.vdb"
  echo "Runing extract_ATOM $inp"
  ${BIN}/extract_ATOM ${inp}

  ## Generate the mesh and simulation config files
  inp="${VDB}_ATOMS.vdb ${T} ${VDW} ${Res} ${Fold} ${inx} ${Fx} ${Fy} ${Fz} ${PDB} ${cone} ${load_ele} ${Young}"
  echo "Running octree_mesh $inp"
  ${BIN}/octree_mesh ${inp}

  ## If simulating a shearing force, the mnesh has to be rotated 90 deg on X
  ##inp="${GEOM} ${GEOM}"
  ##echo "Running shear_rotate $inp"
  ##${BIN}/shear_rotate ${inp}

  ## Run the simulation
  inp="${name} $log_lev"
  echo "Running meshsolver $inp"
  ${BIN}/meshsolver ${inp} > ${OUTF}

  ## Map the simulation results in the mesh into the capsid structure file (virus.pdb)
  inp="${CAPR} $MSH2 $POST ${PDB_results} ${ref_lev}"
  echo "Running mesh2pdb $inp"
  ${BIN}/mesh2pdb ${inp}
  rm ${CAPR}

  ## Rotate capsid structure back to its original orientation (VIPERdb convention)
  inp="${PDB_results} rotate_Z2F.mtx > ${PDB_back}"
  echo "Running apply-matrix.awk $inp"
  ${BIN}/apply-matrix.awk ${PDB_results} rotate_Z2F.mtx > ${PDB_back}

  echo ""
  echo "#############################################################################################"

if [ -e $OUTF ]
then
	mv $POST ${resn}.post.res
	mv $MSH2 ${resn}.post.msh
 echo "Results in $OUTF and ${resn}.post.res to visualize with GID"
 echo "Capsid structure in the VIPERdb convention with nanoindentation results"
 echo "saved in ${PDB_back} to visualize with VMD"
 echo ""
 echo "Job done"
else
	echo "Something went wrong"
fi

  echo "#############################################################################################"



