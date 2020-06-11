#pragma once

//C system files
#include <math.h>
#include <assert.h>

//C++ system files
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>

//Other includes
#include "../geometry/geometry.h"

/**
 *Boundary coordinates scaler
 *This class contains the information to can scale proportionally the boundary mesh
 *
 *@param min_coord_ Is the minimum coordinate of the boundary bounding box
 *@param max_coord_ Is the maximum coordinate of the boundary bounding box
 *@param domain_ Is the length of the bounding box in each dimension
 *@param proportion_ Is the proportion to use by the bounding box on the unit cube
 *@param center_ Is the center of the domains
 */
class Scaler{

	private:
		double min_coord_[ 3 ];
		double max_coord_[ 3 ];
		double s_min_coord_[ 3 ];
		double s_max_coord_[ 3 ];
		double domain_[ 3 ];
		double proportion_[ 3 ];
		double center_[ 3 ];

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Scaler();
		Scaler( double* min_coord , double* max_coord , double percentage );
		~Scaler();
		
		//GETS
		double GetMinSizeOfModelBoundingBox();
    double GetLongestDomain( );

		//SETS
		void SetMinCoord( double x_min , double y_min , double z_min );
		void SetMaxCoord( double x_max , double y_max , double z_max );
		void SetProportions( double percentage );

		//UTILITIES
		double ScaleCoord( int i_pos , double coord );
		double UnscaleCoord( int i_pos , double coord );
		bool IsCoordinateInsideDomain( double coord , int dim );

		//CLEANERS

		//DEBUG
		void PrintScaler();
};

/**
 *Node class
 *This class contains the coordinates from a boundary node
 *
 *@param coords_ Is the set of coordinates from the node
 *@param index_ Is the node index in the list of nodes
 *@param interfaz Is a bool value indicating if a node is interfaz betweenj ranks
 *@param sign is a value indicating the node position respect to boundary
 */
class Node{

	private:
  	double coords_[ 3 ];
		size_t index_;
		bool interfaz_;
		char sign_;
	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Node();
		Node( double* coords , size_t index );
		Node( double x_coord , double y_coord , double z_coord , size_t index );
		~Node(  );

		//GETS
		double GetCoord( int i_pos );
		size_t GetIndex();
		char GetSign();

		//SETS
  	void SetCoord( int i_pos , double val );
		void SetIndex( size_t index );
		void SetInterface( bool interfaz );
		void SetSign( char sign );

		//UTILITIES
		bool IntersectsCell(  double* min_coord , double* max_coord  );
};

/**
 *Boundary triangle class
 *This class contains the information from a boundary triangle
 *
 *@param nodes_ Is the list of the three nodes that form the triangle
 */
class Triangle{

	private:
	  Node* nodes_[ 3 ];

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Triangle();
  	Triangle( Node* nod1 , Node* nod2 , Node* nod3 );
  	~Triangle();
		
		//GETS
		Node* GetNode( int i_node );
		double GetCoord( int i_node , int i_coord );
		void GetAllCoordinates( double** T );
		void GetAllCoordinates( double* T0 , double* T1 , double* T2 );
		size_t GetNodeIndex(  int node  );


		//SETS
		void SetNode( Node* node , int i_pos );		

		//UTILITIES
  	bool Intersects( double* bbox_min , double* bbox_max , double tolerance );
  	bool IntersectsBox( double* bbox_min , double* bbox_max );
  	void CalcBoundingBox( double* min_corner , double* max_corner );
		bool IntersectsLine( double* P0 , double* P1 );
		int GetAllIntersectionsWithLine( double* P0 , double* P1 , double* X , int dim );
		bool ContainsCoord( double* coord );
		bool NodeIntersectsCell( int node , double* min_coord , double* max_coord );	


};

/**
 *Boundary mesh of triangles
 *This class contains the boundary mesh representing the domain to mesh and the scaler 
 *to can map the nodes between the original domain and the unit AABB( Axis Aligned Bounding Box )
 *
 *@param elems_ Is the list of triangles belonging to the boundary mesh
 *@param nodes_ Is the list of nodes belonging to the boundary mesh
 *@param scale_ Is the class to mapp between the AABB and the original domain
 */
class Boundary{

	private:
		std::vector<Triangle*> elems_;
		std::vector<Node*> nodes_;
		Scaler* scale_;

	public:
		//CONSTRUCTOR, DESTRUCTOR
		Boundary();
		Boundary( char* name );
		~Boundary();

		//GETS
		size_t GetNNodes();
		size_t GetNElements();
		Triangle* GetElement( int i_elem );
		Node* GetNode( int i_node );
		double GetNodeCoord( size_t i_node , int i_coord );
		double GetMinSizeOfModelBoundingBox();
		Scaler* GetScaler();

		//SETS
		void SetElement( Triangle* elem );
		void SetNode( Node* node );

		//UTILITIES
		void CreatesScalerClass();
		void ScaleNodes();
		double ScaleCoordinate( int i_pos , double coord );
		double UnscaleCoordinate( int i_pos , double coord );
		char IsOuterPoint( double* coord );
		char PointInPolyhedron( double* coord );
		bool IsCoordInsideDomain( double coord , int dim );

		//DEBUG
		void PrintMesh();
		void PrintScaler();
};




















