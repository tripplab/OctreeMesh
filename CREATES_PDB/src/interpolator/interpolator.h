#pragma once

//C system files
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

//C++ system files
#include <vector>
#include <iostream>
#include <bitset>
#include <cstdio>

//Other libraries
#include "../octree/octree_driver.h"
#include "../results/results.h"


using namespace std;

/**
 *Interpolator class
 *This class contains the methods to interpolate results from GiD files to pdb file
 *@param input_ Is the name of the vdb file containing the mesh information
 *@param vdb_ Is the pointer to the spherical mesh
 *@param mesh_ Is the name of the hexaedral mesh generated
 *@param hex_ Is the class containing the hexaedral mesh
 *@param results_ Is the file that contains the results 
 *@param res_ Is the pointer that contains the results of the capside
 *@param res_nodes_ Is the list of resultas traslated to nodes [ 14 ][ n_nodes ]
 *@param output_ Is the name of the file where the results will be stored
 *@param r_level_ Is the refinement level that the user requires in the mesh
 *@param octree_ Is the octree pointer that the master has created
 *@param local_root_ Is the local cell that each process will use to create its own mesh
 */
class Interpolator{ 

	char* input_;
	Vdb* vdb_;

	char* mesh_;
	Hexahedral* hex_;

	char* results_;
	Results* res_;
	double** res_nodes_;

	char* output_;
	int r_level_;
	OctreeDriver* octree_;
	OctreeCell* local_root_;



	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Interpolator();
		Interpolator( int argc , char** argv );
		~Interpolator();

		//GETS
		void GetAllLocalRootLeaves( OctreeCell_vector& leaves );
		OctreeCell* GetCellOnLocalOctree( key_type* keys );


		//UTILITIES
		void ReadVdbFile(  );
    void ScaleCapside(  );  
		void RefineLocalRoot();
		void AssignAtomsOnLocalRoot( );
		void ScaleHexahedralMesh(  );
		void TranslateAllResultsToNodes(  );
		void AssignHexahedralsOnLocalRoot( );
		void InterpolateResultsToAtomsAndSavePdb();

};







