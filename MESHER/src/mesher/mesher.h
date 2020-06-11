#pragma once

//C system files
#include <math.h>
#include <assert.h>

//C++ system files
#include <vector>
#include <iostream>
#include <cstdio>

//Other libraries
#include "../octree/octree_driver.h"

//defines

using namespace std;

/**
 *Comunicator class
 *This class contains the information requiered to create a mesh based on a octree
 *@param nodes_label Is the value that controls the actual node index setted
 *@param input_ Is the name of the vdb file containing the mesh information
 *@param output_ Is the name of the file where the hexahedral mesh will be stored
 *@param r_level_ Is the refinement level required to generate the capside mesh
 *@param resolution_ Is the required size of the mesh in Amstrongs
 *@param resolution_capside_ Is the equivalent resolution in the capsid mesh mesher
 *@param percentage_ Is the percentage of the octree that will be used in order to hace 
				 the resolution required in the mesh
 *@param octree_ Is the octree pointer that the master has created
 *@param local_root_ Is the local cell that each process will use to create its own mesh
 *@param vdb_ Is the pointer to the spherical mesh
 *@param warnings_ Is the pointer to the warnings file
 *@param msh_ Is the pointer to the hexahedral mesh
 *@param fold_ Is the number of fold to align with the Y axis. Is a number 5,3 or 2
 *@param fold_index_ Is the id of fold to align with the Y axis. Is a number from 0 to 11
 *@param virus Is the name of the virus that is being processed 
 *@param cone_amplitude_ Is the amplitud of the cone to impose conditions
 *@param prop_variation_ Is the change in the proportion of the volume loaded
 *@param young_modulus_ Is the value of young modulus to be used in the simulation
 */
class Mesher{ 

	char* input_;
	char* output_;
	int r_level_;
	int fold_;
	int fold_index_;
	double resolution_;
	double resolution_capside_;
	double percentage_;
	OctreeDriver* octree_;
	OctreeCell* local_root_;
	Vdb* vdb_;
	Mesh* msh_;
	char* virus_;
	double cone_amplitude_;
	double prop_variation_;
	double young_modulus_;
	size_t n_final_nodes_;

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Mesher();
		Mesher( int argc , char** argv );
		~Mesher();

		//GETS
		void GetAllLocalRootLeaves( OctreeCell_vector& leaves );
		OctreeCell* GetCellOnLocalOctree( key_type* keys );
		void GetConeDirection( double* dir );
		void GetConeDirectionFold2( double* dir );
		void GetConeDirectionFold3( double* dir );
		void GetConeDirectionFold5( double* dir );
		size_t GetNNodesInFinalMesh(  );
		size_t GetNElementsInFinalMesh(  );

		//SETS
		void SetAllOctreeNodes( );
		void SetNodesOnCellAndNeighbours( OctreeCell* cell , int* index );
		bool SetNodeInCellAndNeighbors( int node , OctreeCell* cell , double* coord , int index );
		void SetNodeInCell( int node , OctreeCell* cell , double* coord , int index );
		void SetPrintedOnNodeAndNeighbours( int node , OctreeCell* cell );
		void SetPrintedOnFixedNodeAndNeighbours( int node , OctreeCell* cell );
		void SetNodePrintedInCell( int node , OctreeCell* cell );
		void SetFixedNodePrintedInCell( int node , OctreeCell* cell );
		void SetLoadedAndFixedElements(  );
		bool SetLoadedElements( double tot_vol , double prop_vol , double* vertex , 
														double* direction , double cell_prop );
		bool SetFixedElements( double tot_vol , double prop_vol , double* vertex , 
													 double* direction , double cell_prop );

		//UTILITIES
		void ReadVdbFile(  );
		void SaveVdbOnGiDMesh(  );
    void ScaleCapside(  );  
		void RefineLocalRoot();
		void AssignAtomsOnLocalRoot( );
		void CalculateRefinementAndPercentage();
		void RotateCoordianteToFold( double* coord );
		void RotateFold2( double* coord );
		void RotateFold3( double* coord );
		void RotateFold5( double* coord );
		void RotateFold2Id0( double* coord );
		void RotateFold2Id1( double* coord );
		void RotateFold2Id2( double* coord );
		void RotateFold2Id3( double* coord );
		void RotateFold2Id4( double* coord );
		void RotateFold2Id5( double* coord );
		void RotateFold2Id6( double* coord );
		void RotateFold2Id7( double* coord );
		void RotateFold2Id8( double* coord );
		void RotateFold2Id9( double* coord );
		void RotateFold2Id10( double* coord );
		void RotateFold2Id11( double* coord );
		void RotateFold3Id0( double* coord );
		void RotateFold3Id1( double* coord );
		void RotateFold3Id2( double* coord );
		void RotateFold3Id3( double* coord );
		void RotateFold3Id4( double* coord );
		void RotateFold3Id5( double* coord );
		void RotateFold3Id6( double* coord );
		void RotateFold3Id7( double* coord );
		void RotateFold3Id8( double* coord );
		void RotateFold3Id9( double* coord );
		void RotateFold3Id10( double* coord );
		void RotateFold3Id11( double* coord );
		void RotateFold5Id0( double* coord );
		void RotateFold5Id1( double* coord );
		void RotateFold5Id2( double* coord );
		void RotateFold5Id3( double* coord );
		void RotateFold5Id4( double* coord );
		void RotateFold5Id5( double* coord );
		void RotateFold5Id6( double* coord );
		void RotateFold5Id7( double* coord );
		void RotateFold5Id8( double* coord );
		void RotateFold5Id9( double* coord );
		void RotateFold5Id10( double* coord );
		void RotateFold5Id11( double* coord );
		double CalculateOctreeMeshVolume( double* cell_volume );
		double CalculateProportionLoaded(  );
		double CalcDensityForLoadedElements( double masa );
		double InterpolateCapsidResolution(  );
		double InterpolateCapsidResolution1CWP(  );
		double InterpolateCapsidResolution4G93(  );
		double InterpolateCapsidResolution3IZG(  );
		double InterpolateCapsidResolution3J4U(  );

		//SAVING ON FILE
		void PrintDataFilesForFEMT(  );
		void PrintSolverDataFileForFEMT( char* solver , int threads , double tol , 
																		 int max_steps , int preconditioner , 
																		 std::string name );
		void PrintProblemDataFileForFEMT( std::string name );
		void PrintGeometryDataFileForFEMT( std::string name );
};







