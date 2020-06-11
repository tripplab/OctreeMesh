#pragma once

//C system files
#include <math.h>
#include <assert.h>

//C++ system files
#include <vector>
#include <iostream>
#include <string.h>
#include <cstdio>

//Project includes 
#include "../boundary/boundary.h"

/**
 *Tterahedron class
 *This class contains the information from tetrahedron
 *@param nodes_ Is the list of four nodes that form the tetrahedron
 */
class Tetrahedron{

	private:
	  Node* nodes_[ 4 ];

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Tetrahedron();
  	Tetrahedron( Node* nod1 , Node* nod2 , Node* nod3 , Node* nod4 );
  	~Tetrahedron();
		
		//GETS
		size_t GetIndex( int node );

		//SETS

		//UTILITIES
		bool BelongsToEmbededMesh( int* indexes );
		bool BelongsToBodyFittedMesh( int* indexes );

};

/**
 *Tetrahedron mesh
 *This class contains a tetrahedral mesh
 *
 *@param nodes_ Is the list of nodes belonging to the boundary mesh
 *@param elems_ Is the list of tetrahedrons creating the mesh
 *@param scale_ Is the class to mapp between the AABB and the original domain
 */
class Mesh{

	private:
		std::vector<Node*> nodes_;
		std::vector<Tetrahedron*> elems_;
		Scaler* scale_;

	public:
		//CONSTRUCTOR, DESTRUCTOR
		Mesh();
		Mesh( Scaler* scale , int n_nodes );
		~Mesh();

		//GETS
		size_t GetNNodes();
		size_t GetNElements();
		Tetrahedron* GetElement( int i_elem );
		Node* GetNode( int i_node );
		double GetCoord( int node , int dim );
		size_t GetIndex( int i_elem , int i_node );

		//SETS
		void SetElement( Tetrahedron* elem );
		void SetNode( Node* node );
		void SetNode( Node* node , int index );

		//UTILITIES
		bool IsNodeSetted( int node );
		void AddNode( int index , double* coord , bool interfaz , char sign );
		void AddElement( int index0 , int index1 , int index2 , int index3 );
		void GetAllNodesAndElemsOnEmbededMesh(  int* nodes_index , int* elems_index  );
		void GetAllNodesAndElemsOnBodyFittedMesh(  int* nodes_index , int* elems_index  );

		//SAVING MESH ON FILE
		void SaveMeshOnGiDFile( int mpi_rank , char* output_ );
		void SaveEmbededMeshOnGiDFile( int mpi_rank , char* output_ );
		void SaveBodyFittedMeshOnGiDFile( int mpi_rank , char* output_ );


		//DEBUG
		void PrintInfo();
};


/**
 *Class containing a element of type sphere
 *@parameter coords_ Is the spacial position of the sphere
 *@parameter radius_ Is the radius of the sphere
 *@parameter material_ Is the material of the sphere (must be an index from 1-6)
 */
class Sphere{

	private:
		double coords_[ 3 ];
		double radius_;
		int material_;
	public:
		//CONSTRUCTOR, DESTRUCTOR
		Sphere();
		Sphere( double* coords , double radius , int material );
		Sphere( double x_coord , double y_coord , double z_coord , double radius , int material );
		~Sphere();

		//GETS
		double GetCoord( int position );
		void GetCoords( double* coords );
		double GetRadius(  );
		int GetMaterial(  );

		//SETS
		void SetCoord( int position , double coord );
		void SetCoords( double* coords );
		void SetRadius( double radius );
		void SetMaterial( int material );

		//UTILITIES

		//SAVING MESH ON FILE

		//DEBUG
		void PrintInfo();
};

/**
 *Spheres mesh
 *This class contains a mesh composed by spheres
 *@param scale_ Is the class containig all the information required to scale objects
 *@param elements_ Is the list of elements that compose the spherical mesh
 *@param color Is the label that indicates the kind of material of the spheres
 */
class Spheres{

	private:
    Scaler* scale_;
		std::vector<Sphere*> elements_;
		int color_[ 6 ] = { 1/*Carbon*/, 2/*Hydrogen*/ , 3/*Oxygen*/ , 4/*Nitrogen*/ , 5/*phosphorus*/ , 6/*Sulfur*/ };

	public:
		//CONSTRUCTOR, DESTRUCTOR
		Spheres();
		~Spheres();

		//GETS
		size_t GetNElements();
		Sphere* GetElement( size_t index );
		double GetCoord( size_t index , int position );
		void GetCoords( size_t index , double* coords );
		double GetRadius( size_t index );
		int GetMaterial( size_t index );

		//SETS
		void SetCoord( size_t index , int position , double coord );
		void SetCoords( size_t index , double* coord );
		void SetRadius( size_t index , double radius );
		void SetMaterial( size_t index , int material );

		//UTILITIES
		void AddElement( Sphere* element );
    void FillScaler( double percentage );
    void ScaleMesh( );
		double UnscaleCoordinate( int i_pos , double coord );

		//SAVING MESH ON FILE
		void SaveMeshOnGiDFile( char* output_ );

		//DEBUG
		void PrintInfo();
};


