#include "./mesher/mesher.h"

/**
 *This program generates a capside mesh by using a vdb input file
 *The main idea is generate a mesh composed by spheres representing different kind of 
 *atoms and then determine the octree hexahedral elements intersecting the spherical mesh.
 *The parameters receied by the program are:
 *     1-Name of the vdb input file
 *     2-Virus type (T3 , T4 , T7,etc)
 *     3-Type of solution ( 1 = VDW , 2 = SAS )
 *     4-Resolution of the mesh in Amstrongs
 *     5-Number of fold to align with the Y axis (5,3 or 2)
 *     6-Number of id to align with Y axis (0-11)
 *     7-Virus name (1CWP,4G93,3IZG,3J4U)
 *     8-Cone amplitude in degrees (it is converted to radians by the programm)
 *     9-Variation in the proportion of loaded volume (positive or negative)
 *     10-Young modulus used in the simulation
 */

void CheckParameters( int argc );

int main( int argc , char** argv ){

	CheckParameters( argc );

	double t_ini = clock();
	std::cout << " 0/100 START MESHING " << std::endl;
	Mesher* com = new Mesher( argc , argv );
	com->ReadVdbFile(  );
	com->CalculateRefinementAndPercentage();
	com->SaveVdbOnGiDMesh(  );
  com->ScaleCapside(  );
	std::cout << " 20/100 VDB READED " << std::endl;
	com->AssignAtomsOnLocalRoot( );
	std::cout << " 40/100 ATOMS ASSIGNED " << std::endl;
	com->RefineLocalRoot();
	com->SetAllOctreeNodes();
	std::cout << " 60/100 OCTREE REFINED " << std::endl;
	com->SetLoadedAndFixedElements();	
	std::cout << " 80/100 BOUNDARY CONDITIONS SET  " << std::endl;
	com->PrintDataFilesForFEMT(  );
	std::cout << " 100/100 PROBLEMA GUARDADO " << std::endl;
	std::cout << "  Time " << (clock() - t_ini)/((double)CLOCKS_PER_SEC) << std::endl;
	delete com;
  return 0;
}

/**
 *Checking if number of arguments is correct or not
 *@param[in] argc Number of arguments in program
 */
void CheckParameters( int argc ){
	if(  argc != 11  ){
		std::cout << "\n\n";
		std::cout << "**************************************************************************************" << std::endl;
		std::cout << "********************************MESH TO PDB APPLICATION*******************************" << std::endl;
		std::cout << "Received wrong number of arguments " << std::endl; 
		std::cout << "Arguments must be: " << std::endl;
		std::cout << "  --[ 1][string][input ] Name of vdb file containing the atoms information" << std::endl;
		std::cout << "  --[ 2][int   ][input ] Virus type (3 , 5 , 7, etc)" << std::endl;
		std::cout << "  --[ 3][int   ][input ] Solution type (1=VDW, 2=SAS), actually just working VWD " << std::endl;
		std::cout << "  --[ 4][double][input ] Resolution of mesh in Amstrongs" << std::endl;
		std::cout << "  --[ 5][int   ][input ] Number of fold to align with Y axis (valids are 2,3 and 5) " << std::endl;
		std::cout << "  --[ 6][int   ][input ] Number of fold id to align with Y axis (0-11)" << std::endl;
		std::cout << "  --[ 7][string][input ] Name of virus processed (valids are 1CWP,4G93,3IZG,3J4U)" << std::endl;
		std::cout << "  --[ 8][double][input ] Cone amplitude in degrees" << std::endl;
		std::cout << "  --[ 9][double][input ] Variation in the proportion of loaded volume (positive or negative)" << std::endl;
		std::cout << "  --[10][double][input ] Young modulus divided by 10,000" << std::endl;
		std::cout << "********************************MESH TO PDB APPLICATION*******************************" << std::endl;
		std::cout << "**************************************************************************************" << std::endl;
		std::cout << "\n\n";
		assert( 0 );
	}
}


