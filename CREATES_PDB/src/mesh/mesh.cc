#include "mesh.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                   TETRAHEDRON METHODS                                //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Tetrahedron::Tetrahedron(){
}

/**
 *Constructor receiving 4 nodes
 *@param[in] nod1 Is the first node of tetrahedron
 *@param[in] nod2 Is the second node of tetrahedron
 *@param[in] nod3 Is the third node of tetrahedron
 *@param[in] nod4 Is the fourth node of tetrahedron
 */
Tetrahedron::Tetrahedron( Node* nod1 , Node* nod2 , Node* nod3 , Node* nod4 ){
	nodes_[ 0 ] = nod1;
	nodes_[ 1 ] = nod2;
	nodes_[ 2 ] = nod3;
	nodes_[ 3 ] = nod4;
}

/**
 *Defalut destructor
 */
Tetrahedron::~Tetrahedron(){

}

//GETS
/**
 *Getting index from node
 *@param[in] node Is the local position of the node index requested
 *@return A size_t value with the real node index
 */
size_t Tetrahedron::GetIndex( int node ){
	return nodes_[ node ]->GetIndex(  ); 
}

//SETS

//UTILITIES
/**
 *Knowing if a tetrahedron belongs to the embeded mesh or not, a tetrahedron belongs to 
 *the embeded mesh if at less one of its nodes is outside of the model
 *@param[in] indexes Is the array where need to be setted the indexes belonging to embeded mesh
 *@return A bool value indicating if belongs or not to the embeded mesh
 */
bool Tetrahedron::BelongsToEmbededMesh( int* indexes ){
	bool flag = false;
	int sign[ 4 ];
	for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){
		sign[ i_node ] = nodes_[ i_node ]->GetSign();
		if(  (sign[ i_node ] == 1) || (sign[ i_node ] == -5)  ){
			flag = true;
			break;
		}
	}
	if(  flag  ){
		for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){
			indexes[ nodes_[ i_node ]->GetIndex() ] = 1;
		}
	}
	return flag;
}

/**
 *Knowing if a tetrahedron belongs to the embeded mesh or not, a tetrahedron belongs to 
 *the embeded mesh if at less one of its nodes is outside of the model
 *@param[in] indexes Is the array where need to be setted the indexes belonging to embeded mesh
 *@return A bool value indicating if belongs or not to the embeded mesh
 */
bool Tetrahedron::BelongsToBodyFittedMesh( int* indexes ){
	bool flag = false;
	int sign[ 4 ];
	for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){
		sign[ i_node ] = nodes_[ i_node ]->GetSign();
		if(  ( sign[ i_node ] == 0 ) || ( sign[ i_node ] == -1 )  ){
			flag = true;
			break;
		}
	}
	if(  flag  ){
		for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){
			indexes[ nodes_[ i_node ]->GetIndex() ] = 1;
		}
	}
	return flag;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                      MESH METHODS                                    //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Mesh::Mesh(){
	scale_ = NULL;
}

/**
 *Constructor receiving a scaler class
 *@param[in] scale Is a scaler class instanced
 */
Mesh::Mesh( Scaler* scale , int n_nodes ){
	scale_ = scale;
	nodes_.assign( n_nodes , NULL );
}

/** 
 *Default destructor
 */
Mesh::~Mesh(){
	if( scale_ ){
		delete scale_;
		scale_ = NULL;
	}
	elems_.clear();
	nodes_.clear();
}



//GETS
/**
 *Getting number of nodes on mesh
 *@return A size_t value indicating the number of nodes on the mesh
 */
size_t Mesh::GetNNodes(){
	return nodes_.size();
}

/**
 *Getting number of elements on mesh
 *@return A size_t value indicating the amount of elements on mesh
 */
size_t Mesh::GetNElements(){
	return elems_.size();
}

/**
 *Getting element from mesh
 *@param[in] i_elem Is the index element requested
 *@return A tetrahedron pointer to the element requested
 */
Tetrahedron* Mesh::GetElement( int i_elem ){
	return elems_[ i_elem ];
}

/**
 *Getting node from mesh
 *@param[in] i_node Is the index of node requested
 *@return A node pointer to the node requested
 */
Node* Mesh::GetNode( int i_node ){
	return nodes_[ i_node ];
}

/**
 *Getting coordinate from node
 *@param[in] node Is the node index 
 *@param[in] dim Is the dimention where the coordinate is requested 0-X, 1-Y and 2-Z
 *@return A double value with the coordinate requested
 */
double Mesh::GetCoord( int node , int dim ){
	return nodes_[ node ]->GetCoord( dim );
}

/**
 *Getting index from node
 *@param[in] i_eleme Is the element index on the list
 *@param[in] i_node Is the node position on element to obtain its index
 *@return A size_t value with the index requested
 */
size_t Mesh::GetIndex( int i_elem , int i_node ){
	return elems_[ i_elem ]->GetIndex( i_node );
}

//SETS
/**
 *Setting element on mesh
 *@param[in]  elem Is the pointer to the element that will be set
 */
void Mesh::SetElement( Tetrahedron* elem ){
	elems_.push_back( elem );
}

/**
 *Setting node on mesh
 *@param[in] node is the pointer to the node that will be set
 */
void Mesh::SetNode( Node* node ){
	nodes_.push_back( node );
}

/**
 *Setting node on a position of the list
 *@param[in] node Is the pointer to the node that will be setted
 *@param[in] index Is the position where the node will be setted
 */
void Mesh::SetNode( Node* node , int index ){
	nodes_[ index ] = node;
}

//UTILITIES
/**
 *Testing if a mesh node is setted or not
 *@param[in] node Is the node index to be testet to know if has been setted or not
 *@return A bool value indicating if node is setted or not
 */
bool Mesh::IsNodeSetted( int node ){
	return nodes_[ node ];
}

/**
 *Adding a Node on tetrahedral mesh
 *@param[in] index Is the index of the node to be added
 *@param[in] coord Are the coordinates of the node to be added on the mesh
 */
void Mesh::AddNode( int index , double* coord , bool interfaz , char sign ){
	Node* node = new Node( coord , index );
	node->SetInterface( interfaz );
	node->SetSign( sign );
	this->SetNode( node , index );
}

/**
 *Adding element to the mesh
 @param[in] index0 Is the first node of the tetrahedron to be set
 @param[in] index1 Is the first node of the tetrahedron to be set
 @param[in] index2 Is the first node of the tetrahedron to be set
 @param[in] index3 Is the first node of the tetrahedron to be set
 */
void Mesh::AddElement( int index0 , int index1 , int index2 , int index3 ){
	Node* node0 =  this->GetNode( index0 );
	Node* node1 =  this->GetNode( index1 );
	Node* node2 =  this->GetNode( index2 );
	Node* node3 =  this->GetNode( index3 );
	Tetrahedron* tet = new Tetrahedron( node0 , node1 , node2 , node3 );
	this->SetElement( tet );
}

/**
 *Getting all nodes and elements on a embeded mesh, in this method also are actualized the
 *indexes of the elements and nodes belonging to the embeded mesh
 *@param[in] nodes_index Is the array where the node indexes will be actualized
 *@param[in] elems_index Is the array where the elements indexes will be actualized
 */
void Mesh::GetAllNodesAndElemsOnEmbededMesh(  int* nodes_index , int* elems_index  ){
	size_t n_nodes = this->GetNNodes();
	size_t n_elems = this->GetNElements();
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		if(  elems_[ i_elem ]->BelongsToEmbededMesh( nodes_index )  ){
			elems_index[ i_elem ] = 1;
		}
	}
	int cont = 1;
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if(  nodes_index[ i_node ] == 1  ){
			nodes_index[ i_node ] = cont;
			cont++;
		}
	}
	cont = 1;
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		if(  elems_index[ i_elem ] == 1  ){
			elems_index[ i_elem ] = cont;
			cont++;
		}
	}
}

/**
 *Getting all nodes and elements on a body fitted mesh, in this method also are actualized 
 *the indexes of the elements and nodes belonging to the embeded mesh
 *@param[in] nodes_index Is the array where the node indexes will be actualized
 *@param[in] elems_index Is the array where the elements indexes will be actualized
 */
void Mesh::GetAllNodesAndElemsOnBodyFittedMesh(  int* nodes_index , int* elems_index  ){
	size_t n_nodes = this->GetNNodes();
	size_t n_elems = this->GetNElements();
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		if(  elems_[ i_elem ]->BelongsToBodyFittedMesh( nodes_index )  ){
			elems_index[ i_elem ] = 1;
		}
	}
	int cont = 1;
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if(  nodes_index[ i_node ] == 1  ){
			nodes_index[ i_node ] = cont;
			cont++;
		}
	}
	cont = 1;
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		if(  elems_index[ i_elem ] == 1  ){
			elems_index[ i_elem ] = cont;
			cont++;
		}
	}
}

//SAVING MESH ON FILE
/**
 *Saving mesh on GiD File format
 *@param[in] mpi_rank Is the process that will save the mesh
 *@param[in] out_put Is the output file name to save the mesh
 */
void Mesh::SaveMeshOnGiDFile( int mpi_rank , char* output_ ){
	FILE* fp = NULL;
	char name[1000];
	char rank[10];
	sprintf( rank , "%d" , mpi_rank );
	strcpy( name , rank );
	strcat( name , "_tet_" );
	strcat( name , output_ );

	fp = fopen( name , "w" );	
	if(  !fp  ){
		std::cout << " Error, output file not opened " << std::endl;
		std::cout << name << std::endl;
		assert( fp );
	}	

	fprintf( fp , "MESH \"tet\" dimension 3 ElemType Tetrahedra Nnode 4\n" );
	fprintf( fp , "# color 96 96 96\n" );
	fprintf( fp , "Coordinates\n" );
	fprintf( fp , "# node number coordinate_x coordinate_y coordinate_z  \n" );
	size_t n_nodes = this->GetNNodes();
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		fprintf( fp , "%d ", (int)(i_node+1) );
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			double coord = this->GetCoord( (int)i_node , i_dim );
			fprintf( fp , "%lf ", scale_->UnscaleCoord( i_dim , coord ) );
		}
		fprintf( fp, "\n" ); 
	}

	fprintf( fp , "end coordinates\n" );
	fprintf( fp , "Elements\n" );
	fprintf( fp , "# element node_1 node_2 node_3 node_4 material_number\n" );

	size_t n_elems = this->GetNElements();

	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		fprintf( fp , "%d ", (int)(i_elem+1) );
		for(  int i_node = 0 ;  i_node < 4  ;  i_node++  ){
			size_t index = this->GetIndex( i_elem , i_node ) +1;
			fprintf( fp , "%d ", (int)index );
		}
		fprintf( fp , "%d\n", 100 );
	}
	fprintf( fp , "end elements\n" );
	fclose( fp );
}

/**
 *Saving embeded mesh on GiD File format
 *@param[in] mpi_rank Is the process that will save the mesh
 *@param[in] out_put Is the output file name to save the mesh
 */
void Mesh::SaveEmbededMeshOnGiDFile( int mpi_rank , char* output_ ){
	size_t n_nodes = this->GetNNodes();
	size_t n_elems = this->GetNElements();
	int* nodes_index = new int [ n_nodes ];
	int* elems_index = new int [ n_elems ];
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		nodes_index[ i_node ] = 0;
	}
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		elems_index[ i_elem ] = 0;
	}
	this->GetAllNodesAndElemsOnEmbededMesh(  nodes_index , elems_index  );
	FILE* fp = NULL;
	char name[1000];
	char rank[10];
	sprintf( rank , "%d" , mpi_rank );
	strcpy( name , rank );
	strcat( name , "_tet_embeded_" );
	strcat( name , output_ );

	fp = fopen( name , "w" );	
	if(  !fp  ){
		std::cout << " Error, output file not opened " << std::endl;
		std::cout << name << std::endl;
		assert( fp );
	}	

	fprintf( fp , "MESH \"tet\" dimension 3 ElemType Tetrahedra Nnode 4\n" );
	fprintf( fp , "# color 96 96 96\n" );
	fprintf( fp , "Coordinates\n" );
	fprintf( fp , "# node number coordinate_x coordinate_y coordinate_z  \n" );
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if( nodes_index[ i_node ] > 0 ){
			fprintf( fp , "%d ", nodes_index[ i_node ] );
			for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
				double coord = this->GetCoord( (int)i_node , i_dim );
				fprintf( fp , "%lf ", scale_->UnscaleCoord( i_dim , coord ) );
			}
			fprintf( fp, "\n" ); 
		}
	}

	fprintf( fp , "end coordinates\n" );
	fprintf( fp , "Elements\n" );
	fprintf( fp , "# element node_1 node_2 node_3 node_4 material_number\n" );

	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		if(  elems_index[ i_elem ] > 0  ){
			fprintf( fp , "%d ", elems_index[ i_elem ] );
			for(  int i_node = 0 ;  i_node < 4  ;  i_node++  ){
				size_t index = this->GetIndex( i_elem , i_node );
				fprintf( fp , "%d ", nodes_index[ index ] );
			}
			fprintf( fp , "%d\n", 100 );
		}
	}

	fprintf( fp , "end elements\n" );
	fclose( fp );



	delete[] nodes_index;
	delete[] elems_index;
}

/**
 *Saving body fittedmesh on GiD File format
 *@param[in] mpi_rank Is the process that will save the mesh
 *@param[in] out_put Is the output file name to save the mesh
 */
void Mesh::SaveBodyFittedMeshOnGiDFile( int mpi_rank , char* output_ ){
	size_t n_nodes = this->GetNNodes();
	size_t n_elems = this->GetNElements();
	int* nodes_index = new int [ n_nodes ];
	int* elems_index = new int [ n_elems ];
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		nodes_index[ i_node ] = 0;
	}
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		elems_index[ i_elem ] = 0;
	}
	this->GetAllNodesAndElemsOnBodyFittedMesh(  nodes_index , elems_index  );
	FILE* fp = NULL;
	char name[1000];
	char rank[10];
	sprintf( rank , "%d" , mpi_rank );
	strcpy( name , rank );
	strcat( name , "_tet_body_" );
	strcat( name , output_ );

	fp = fopen( name , "w" );	
	if(  !fp  ){
		std::cout << " Error, output file not opened " << std::endl;
		std::cout << name << std::endl;
		assert( fp );
	}	

	fprintf( fp , "MESH \"tet\" dimension 3 ElemType Tetrahedra Nnode 4\n" );
	fprintf( fp , "# color 96 96 96\n" );
	fprintf( fp , "Coordinates\n" );
	fprintf( fp , "# node number coordinate_x coordinate_y coordinate_z  \n" );
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if( nodes_index[ i_node ] > 0 ){
			fprintf( fp , "%d ", nodes_index[ i_node ] );
			for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
				double coord = this->GetCoord( (int)i_node , i_dim );
				fprintf( fp , "%lf ", scale_->UnscaleCoord( i_dim , coord ) );
			}
			fprintf( fp, "\n" ); 
		}
	}

	fprintf( fp , "end coordinates\n" );
	fprintf( fp , "Elements\n" );
	fprintf( fp , "# element node_1 node_2 node_3 node_4 material_number\n" );

	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		if(  elems_index[ i_elem ] > 0  ){
			fprintf( fp , "%d ", elems_index[ i_elem ] );
			for(  int i_node = 0 ;  i_node < 4  ;  i_node++  ){
				size_t index = this->GetIndex( i_elem , i_node );
				fprintf( fp , "%d ", nodes_index[ index ] );
			}
			fprintf( fp , "%d\n", 100 );
		}
	}
	fprintf( fp , "end elements\n" );
	fclose( fp );



	delete[] nodes_index;
	delete[] elems_index;
}


//DEBUG
void Mesh::PrintInfo(){
	size_t n_nodes = this->GetNNodes();
	size_t n_elems = this->GetNElements();
	std::cout<<" Nodes = "<<n_nodes<<std::endl;
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if(  nodes_[ i_node ] != NULL  ){
			std::cout<<i_node<<"-"<<this->GetCoord( i_node , 0 )<<" "<<this->GetCoord( i_node , 1 )<<" "<<this->GetCoord( i_node , 2 )<<std::endl;
		}else{
			std::cout<<i_node<<" Nodes not setted"<<std::endl;
			//break;
		}
	}
	std::cout<<" Elements = "<<n_elems<<std::endl;
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		std::cout<<i_elem<<"-"<<this->GetIndex( i_elem , 0 )<<" "<<this->GetIndex( i_elem , 1 )<<" "<<this->GetIndex( i_elem , 2 )<<" "<<this->GetIndex( i_elem , 3 )<<std::endl;
	}
	scale_->PrintScaler();
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                    SPHERE METHODS                                    //

//CONSTRUCTOR, DESTRUCTOR
/**
 *Default constructor
 */
Sphere::Sphere(){
	coords_[ 0 ] = 0.0;
	coords_[ 1 ] = 0.0;
	coords_[ 2 ] = 0.0;
	radius_ = 0.0;
	material_ = 0;
}

/**
 *Second constructor
 *@param[in] coords Is the array with the sphere coordinates
 *@param[in] radius Is the radius of the sphere
 *@param[in] material Is the material index of the sphere
 */
Sphere::Sphere( double* coords , double radius , int material ){
	coords_[ 0 ] = coords[ 0 ];
	coords_[ 1 ] = coords[ 1 ];
	coords_[ 2 ] = coords[ 2 ];
	radius_ = radius;
	material_ = material;
}

/**
 *Third constructor
 *@param[in] x_coord Is the x_coordinate
 *@param[in] y_coord Is the y_coordinate
 *@param[in] z_coord Is the z_coordinate
 *@param[in] radius Is the radius of the sphere
 *@param[in] material Is the material index of the sphere
 */
Sphere::Sphere( double x_coord , double y_coord , double z_coord , double radius , int material ){
	coords_[ 0 ] = x_coord;
	coords_[ 1 ] = y_coord;
	coords_[ 2 ] = z_coord;
	radius_ = radius;
	material_ = material;
}

/**
 *Default destructor
 */
Sphere::~Sphere(){
	coords_[ 0 ] = 0.0;
	coords_[ 1 ] = 0.0;
	coords_[ 2 ] = 0.0;
	radius_ = 0.0;
	material_ = 0;
}

//GETS
/**
 *Getting coordinate from sphere
 *@param[in] position Indicates the dimention where the coordinate is requested X,Y,Z
 *@return A double value with the coordinate requested
 */
double Sphere::GetCoord( int position ){
	assert( ( position >= 0 ) && ( position < 3 ) );
	return coords_[ position ];
}

/**
 *Geting coordinates from sphere
 *@param[out] coords Is the array where coordinates will be stored
 */
void Sphere::GetCoords( double* coords ){
	for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
		coords[ i_pos ] = coords_[ i_pos ];
	}
}

/**
 *Getting radius from sphere
 *@return A double value with the radius
 */
double Sphere::GetRadius(  ){
	return radius_;
}

/**
 *Getting material from sphere
 *@return A int value with the material index
 */
int Sphere::GetMaterial(  ){
	return material_;
}

//SETS
/**
 *Setting coordinate on sphere
 *@param[in] position Is the dimention where the coordinate will be set
 *@param[in] coord Is the coordinate to be set
 */
void Sphere::SetCoord( int position , double coord ){
	assert( ( position >= 0 ) && ( position < 3 ) );
	coords_[ position ] = coord;
}

/**
 *Setting all coordinates on sphere
 */
void Sphere::SetCoords( double* coords ){
	for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
		coords_[ i_pos ] = coords[ i_pos ];
	}
}

/**
 *Setting radius on sphere
 *@param radius Is the radius to be set on the sphere
 */
void Sphere::SetRadius( double radius ){
	radius_ = radius;
}

/**
 *Setting material on sphere
 */
void Sphere::SetMaterial( int material ){
	material_ =  material;
}

//UTILITIES

//DEBUG
void Sphere::PrintInfo(){
	std::cout << "Printing sphere information " << std::endl;
	std::cout << "   X coordinate : " << coords_[ 0 ] << std::endl;
	std::cout << "   Y coordinate : " << coords_[ 1 ] << std::endl;
	std::cout << "   Z coordinate : " << coords_[ 2 ] << std::endl;
	std::cout << "   Radius       : " << radius_      << std::endl;
	std::cout << "   Material     : " << material_    << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                   SPHERES METHODS                                    //

//CONSTRUCTOR, DESTRUCTOR
/**
 *Default constructor
 */
Spheres::Spheres(){
  scale_ = NULL;
	elements_.clear();
}

/**
 *Default destructor
 */
Spheres::~Spheres(){
  scale_ = NULL;
	elements_.clear();
}

//GETS
/**
 *Getting number of elements on the mesh
 *@return A size_t value with the amount of elements
 */
size_t Spheres::GetNElements(){
	return elements_.size(  );
}

/**
 *Getting element from spherical mesh
 *@param[in] index Is the index of the element requested
 *@return A pointer to the Sphere class requested
 */
Sphere* Spheres::GetElement( size_t index ){
	assert( ( index >= 0 )  &&  ( index < this->GetNElements() ) );
	return elements_[ index ];

}

/**
 *Getting coordinate from sphere on mesh
 *@param[in] index Is the sphere index on the elements array
 *@param[in] position Is the dimention of the requested coordinate
 *@return A double value with the coordinate requested
 */
double Spheres::GetCoord( size_t index , int position ){
	assert( ( index >= 0 )  &&  ( index < this->GetNElements() ) );
	assert( ( position >= 0 )  &&  ( position < 3 ) );
	return elements_[ index ]->GetCoord( position );
}

/**
 *Getting coordinates from sphere on mesh
 *@param[in] index Is the sphere index on the elements array
 *@param[out] coords Is the array to store the sphere coordinates
 */
void Spheres::GetCoords( size_t index , double* coords ){
	assert( ( index >= 0 )  &&  ( index < this->GetNElements() ) );
	elements_[ index ]->GetCoords( coords );
}

/**
 *Getting radius from sphere
 *@param[in] index Is he sphere index on the elements array
 *@return A double value with the radius of the sphere
 */
double Spheres::GetRadius( size_t index ){
	assert( ( index >= 0 )  &&  ( index < this->GetNElements() ) );
	return elements_[ index ]->GetRadius();
}

/**
 *Getting material from sphere
 *@param[in] indes Is the sphere index on the elements array
 *@return A int value with the material index of the sphere
 */
int Spheres::GetMaterial( size_t index ){
	assert( ( index >= 0 )  &&  ( index < this->GetNElements() ) );
	return elements_[ index ]->GetMaterial();
}

//SETS
/**
 *Setting coordinate on sphere
 *@param[in] index Is the sphere index on the spheres array
 *@param[in] position Is the dimention where the coordinate will be set
 *@param[in] coord  Is the value of the coordinate to be set
 */
void Spheres::SetCoord( size_t index , int position , double coord ){
	assert( ( index >= 0 )  &&  ( index < this->GetNElements() ) );
	assert( ( position >= 0 )  &&  ( position < 3 ) );
	elements_[ index ]->SetCoord( position , coord );
}

/**
 *Setting coordinates on sphere
 *@param[in] index Is the sphere index on the spheres array
 *@param[in] coord Is the array with the coordinates to be set
 */
void Spheres::SetCoords( size_t index , double* coord ){
	assert( ( index >= 0 )  &&  ( index < this->GetNElements() ) );
	elements_[ index ]->SetCoords( coord );
}

/**
 *Setting radius on sphere
 *@param[in] index Is the sphere index on the spheres array
 *@param[in] radius Is the radius value to be set
 */
void Spheres::SetRadius( size_t index , double radius ){
	assert( ( index >= 0 )  &&  ( index < this->GetNElements() ) );
	elements_[ index ]->SetRadius( radius );
}

/**
 *Setting material on sphere
 *@param[in] index Is the sphere index on the spheres array
 *@param[in] material Is the material value to be set
 */
void Spheres::SetMaterial( size_t index , int material ){
	assert( ( index >= 0 )  &&  ( index < this->GetNElements() ) );
	elements_[ index ]->SetMaterial( material );
}

//UTILITIES
/**
 *Adding element to the list of elements
 *@param[in] element Is the pointer to the element to be added
 */
void Spheres::AddElement( Sphere* element ){
	elements_.push_back( element );
}

/**
 *This method fills the scaler information in order to use it to scale and unscale
 */
void Spheres::FillScaler( ){
	double min_coord[ 3 ] = {  1e20 ,  1e20 ,  1e20 };
	double max_coord[ 3 ] = { -1e20 , -1e20 , -1e20 };
	size_t n_elems = this->GetNElements();
	//Obtaining minimum and maximum
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
    double radius = this->GetRadius( i_elem ); 
		for(  int i_coord = 0  ;  i_coord < 3  ;  i_coord++  ){
			double coord = this->GetCoord( i_elem , i_coord );
			if(  ( coord - radius ) < min_coord[ i_coord ]  ){
				min_coord[ i_coord ] = ( coord - radius );
			}
			if(  ( coord + radius ) > max_coord[ i_coord ]  ){
				max_coord[ i_coord ] = ( coord + radius );
			}
		}
	}
	scale_ = new Scaler( min_coord , max_coord );	
}

/**
 *This method scales the coordinates of the spherical mesh
 */
void Spheres::ScaleMesh( ){
  double radius_proportion = scale_->GetLongestDomain() / 0.99;

	size_t n_elements = this->GetNElements();
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		Sphere* aux = this->GetElement( i_elem );
		for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
			aux->SetCoord(  i_pos , scale_->ScaleCoord( i_pos , aux->GetCoord( i_pos ) )  );
		}
    aux->SetRadius( aux->GetRadius() / radius_proportion );
	}
}

/**
 *Unscaling coordinate
 *@param[in] i_pos Is the domain where the coordinate is located
 *@param[in] coord Is the coordinate value to be unscaled
 *@return A double value with he coordinate unscaled
 */
double Spheres::UnscaleCoordinate( int i_pos , double coord ){
	double value;
	value = scale_->UnscaleCoord( i_pos , coord );
	return value;
}

/**
 *Rotating the mesh in order to have fold 5 aligned with the Y axis
 *@param[in] fold Is the index (0-11) indiccating the fold to use in the rotation
 */
void Spheres::Rotate( ){
	double teta = -1.0187486844;
	double eta = -0.3141592654;
	//double eta = 0.0;
	

	double x_rotation[ 3 ][ 3 ];
	x_rotation[ 0 ][ 0 ] = 1;	x_rotation[ 0 ][ 1 ] = 0         ;	x_rotation[ 0 ][ 2 ] = 0          ;
	x_rotation[ 1 ][ 0 ] = 0;	x_rotation[ 1 ][ 1 ] = cos(teta);	x_rotation[ 1 ][ 2 ] = -sin(teta);
	x_rotation[ 2 ][ 0 ] = 0;	x_rotation[ 2 ][ 1 ] = sin(teta);	x_rotation[ 2 ][ 2 ] = cos(teta) ;

	double y_rotation[ 3 ][ 3 ];
	y_rotation[ 0 ][ 0 ] = cos(eta)  ;	y_rotation[ 0 ][ 1 ] = 0;	y_rotation[ 0 ][ 2 ] = sin(eta)          ;
	y_rotation[ 1 ][ 0 ] = 0          ;	y_rotation[ 1 ][ 1 ] = 1;	y_rotation[ 1 ][ 2 ] = 0;
	y_rotation[ 2 ][ 0 ] = -sin(eta) ;	y_rotation[ 2 ][ 1 ] = 0;	y_rotation[ 2 ][ 2 ] = cos(eta) ;

	size_t n_elements = this->GetNElements();
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		double coords[ 3 ];
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			coords[ i_dim ] = this->GetCoord( i_elem , i_dim );
		}
		double new_coord[ 3 ];
		for(  int i_row = 0  ;  i_row < 3  ;  i_row++  ){
			new_coord[ i_row ] = 0.0;
			for(  int i_col = 0  ;  i_col < 3  ;  i_col++  ){
				new_coord[ i_row ] += ( x_rotation[ i_row ][ i_col ] * coords[ i_col ] );
			}
		}
		double new_new_coord[ 3 ];
		for(  int i_row = 0  ;  i_row < 3  ;  i_row++  ){
			new_new_coord[ i_row ] = 0.0;
			for(  int i_col = 0  ;  i_col < 3  ;  i_col++  ){
				new_new_coord[ i_row ] += ( y_rotation[ i_row ][ i_col ] * new_coord[ i_col ] );
			}
		}
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			this->SetCoord( i_elem , i_dim , new_new_coord[ i_dim ] );
		}
	}
}

/**
 *Scaling coordinate
 *@param[in] i_pos Is the index of the dimention where the coordinate is scaled 0-X, 1-Y, 2-Z
 *@param[in] coord Is the coordinate to be scaled
 *@return A double value with the coordinate scaled
 */
double Spheres::ScaleCoord( int i_pos , double coord ){
	double coordinate = scale_->ScaleCoord( i_pos , coord );
	return coordinate;
}

//SAVING MESH ON FILE
/**
 *Saving mesh as a GiD mesh
 *@param[in] output_ Is the name of the file to be printed
 */
void Spheres::SaveMeshOnGiDFile( char* output_ ){
	FILE* fp = NULL;
	fp = fopen( output_ , "w" );	
	if(  !fp  ){
		std::cout << " Error, output file not opened: " << output_ << std::endl;
		assert( fp );
	}	

	fprintf( fp , "MESH \"capside\" dimension 3 ElemType Sphere Nnode 1\n" );
	fprintf( fp , "# color 96 96 96\n" );
	fprintf( fp , "Coordinates\n" );
	fprintf( fp , "# node_number coordinate_x coordinate_y coordinate_z  \n" );
	size_t n_elems = this->GetNElements( );
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		fprintf( fp , "%d ", (int)( i_elem + 1 ) );
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			double coord = elements_[ i_elem ]->GetCoord( i_dim );
			fprintf( fp , "%lf ", coord );
		}
		fprintf( fp, "\n" ); 
	}
	fprintf( fp , "end coordinates\n" );
	fprintf( fp , "Elements\n" );
	fprintf( fp , "# element node_1 radius material_number\n" );

	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		fprintf( fp , "%d %d ", (int)( i_elem + 1 ) , (int)( i_elem + 1 ) );
		fprintf( fp , "%lf ", elements_[ i_elem ]->GetRadius( ) );
		fprintf( fp , "%d\n", elements_[ i_elem ]->GetMaterial( ) );
	}
	fprintf( fp , "end elements\n" );
	fclose( fp );
}

//DEBUG
/**
 *Printing mesh information
 */
void Spheres::PrintInfo(){
	printf( "Printing spherical mesh information\n" );
	printf( "Number of elements : %ld\n" , this->GetNElements() );
	for(  size_t i_element = 0  ;  i_element < this->GetNElements()  ;  i_element++  ){
		elements_[ i_element ]->PrintInfo();
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                   HEXAHEDRON METHODS                                 //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Hexahedron::Hexahedron(){
	for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
		nodes_[ i_node ] =  NULL;
	}
}

/**
 *Constructor receiving the list of nodes
 *@param[in] nod1-8 Is the list of nodes forming the element
 */
Hexahedron::Hexahedron( Node* nod1 , Node* nod2 , Node* nod3 , Node* nod4 , Node* nod5 , Node* nod6 , Node* nod7 , Node* nod8 ){
	nodes_[ 0 ] = nod1;
	nodes_[ 1 ] = nod2;
	nodes_[ 2 ] = nod3;
	nodes_[ 3 ] = nod4;
	nodes_[ 4 ] = nod5;
	nodes_[ 5 ] = nod6;
	nodes_[ 6 ] = nod7;
	nodes_[ 7 ] = nod8;
}
 
/**
 *Default destructor
 */
Hexahedron::~Hexahedron(){
	for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
		nodes_[ i_node ] =  NULL;
	}
}
		
//GETS
/**
 *Getting real index of node
 */
size_t Hexahedron::GetIndex( int node ){
	return nodes_[ node ]->GetIndex();
}

/**
 *Getting coordinate fron node
 *@param[in] node Is the node index on the local numeration
 *@param[in] dim Is the dimention where the coordinate is requested 0-X,1-Y, 2-Z
 */
double Hexahedron::GetCoord( int node , int dim ){
	assert( ( node >= 0 ) && ( node < 8 ) );
	assert( ( dim >= 0 )  && ( dim < 3 ) );
	return ( nodes_[ node ]->GetCoord( dim ) );
}

/**
 *Getting the min coord of the bounding box of the hexahedron
 *@param[out] min_coord Pointer where the mimimum coordinate in X,Y and Z will be saved
 */
void Hexahedron::GetMinPoint( double* min_coord ){
	min_coord[ 0 ] = 1e20;
	min_coord[ 1 ] = 1e20;
	min_coord[ 2 ] = 1e20;
	for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			double coord = nodes_[ i_node ]->GetCoord( i_dim );
			if(  coord < min_coord[ i_dim ]  ){
				min_coord[ i_dim ] = coord; 
			}
		}
	}
}

/**
 *Getting the max coord of the bounding box of the hexahedron
 *@param[out] max_coord Pointer where the maximum coordinate in X,Y and Z will be saved
 */
void Hexahedron::GetMaxPoint( double* max_coord ){
	max_coord[ 0 ] = -1e20;
	max_coord[ 1 ] = -1e20;
	max_coord[ 2 ] = -1e20;
	for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			double coord = nodes_[ i_node ]->GetCoord( i_dim );
			if(  coord > max_coord[ i_dim ]  ){
				max_coord[ i_dim ] = coord; 
			}
		}
	}
}

/**
 *Getting normalized coordinates
 *@param[in] X Is the X coordinate in the cartesian space
 *@param[in] Y Is the Y coordinate in the cartesian space
 *@param[in] Z Is the Z coordinate in the cartesian space
 *@param[in] psi Is the reference to the psi coordinate in the normalized space
 *@param[in] eta Is the reference to the eta coordinate in the normalized space
 *@param[in] mu Is the reference to the mu coordinate in the normalized space
 */
void Hexahedron::GetNormalizedCoordinates( double X , double Y , double Z , double* psi , double* eta, double* mu ){
	int MaxIter = 10000, iter = 0;
	double epsilon = 1e-12;
	double delta = 1e-12;
	double tolx = 2.0*epsilon;
	double tolf = 2.0*delta;
	(*psi) = 0.0;
	(*eta) = 0.0;
	(*mu) = 0.0;
	double** J;
	double** invJ;
	J = new double*[ 3 ];
	invJ = new double*[ 3 ];
	for(  int i_row = 0  ;  i_row < 3  ;  i_row++  ){
		J[ i_row ] = new double[ 3 ];
		invJ[ i_row ] = new double[ 3 ];
	}
	double Fx[ 3 ] , dx[ 3 ];
	while( iter < MaxIter ){
		J[ 0 ][ 0 ] = InterpolationDPsi( (*psi) , (*eta) , (*mu) , 0 );
		J[ 0 ][ 1 ] = InterpolationDEta( (*psi) , (*eta) , (*mu) , 0 );
		J[ 0 ][ 2 ] = InterpolationDMu( (*psi) , (*eta) , (*mu) , 0 );
		J[ 1 ][ 0 ] = InterpolationDPsi( (*psi) , (*eta) , (*mu) , 1 );
		J[ 1 ][ 1 ] = InterpolationDEta( (*psi) , (*eta) , (*mu) , 1 );
		J[ 1 ][ 2 ] = InterpolationDMu( (*psi) , (*eta) , (*mu) , 1 );
		J[ 2 ][ 0 ] = InterpolationDPsi( (*psi) , (*eta) , (*mu) , 2 );
		J[ 2 ][ 1 ] = InterpolationDEta( (*psi) , (*eta) , (*mu) , 2 );
		J[ 2 ][ 2 ] = InterpolationDMu( (*psi) , (*eta) , (*mu) , 2 );
		Fx[ 0 ] =  InterpolationFunction( (*psi) , (*eta) , (*mu) , 0 , X );
		Fx[ 1 ] = InterpolationFunction( (*psi) , (*eta) , (*mu) , 1 , Y );
		Fx[ 2 ] = InterpolationFunction( (*psi) , (*eta) , (*mu) , 2 , Z );
		InvertJacobiano( J , invJ );
		for(  int i_row = 0  ;  i_row < 3  ;  i_row++  ){
			dx[ i_row ] = 0.0;		
			for(  int i_col = 0  ;  i_col < 3  ;  i_col++  ){
				dx[ i_row ] += invJ[ i_row ][ i_col ] * Fx[ i_col ];
			}
		}
		(*psi) -= dx[ 0 ];
		(*eta) -= dx[ 1 ];
		(*mu) -= dx[ 2 ];
		if( (*psi) > 1.0 ){
			(*psi) = 1.0;
		}
		if( (*psi) < -1.0 ){
			(*psi) = -1.0;
		}
		if( (*eta) > 1.0 ){
			(*eta) = 1.0;
		}
		if( (*eta) < -1.0 ){
			(*eta) = -1.0;
		}
		if( (*mu) > 1.0 ){
			(*mu) = 1.0;
		}
		if( (*mu) < -1.0 ){
			(*mu) = -1.0;
		}
		tolx = sqrt( dx[ 0 ] * dx[ 0 ] + dx[ 1 ] * dx[ 1 ] + dx[ 2 ] * dx[ 2 ] );
		tolf = sqrt( Fx[ 0 ] * Fx[ 0 ] + Fx[ 1 ] * Fx[ 1 ] + Fx[ 2 ] * Fx[ 2 ] );
		if(  ( tolx < epsilon ) || ( tolf < delta )  ){
			break;
		}
		iter++;
	}
	for(  int i_row = 0  ;  i_row < 3  ;  i_row++  ){
		delete[] J[ i_row ];
		delete[] invJ[ i_row ];
	}
	delete[] J;
	delete[] invJ;
}
//SETS

//UTILITIES
/**
 *Getting the evaluation of the derivates of the form functions
 *param[in] id Is the index of the form function dferivative to be evaluated
 *@param[in] dim Is the index of the derivative variable 0-PSI,1-ETA,2-MU
 *@param[in] psi Is the value of psi coordinate on the normalized space
 *@param[in] eta Is the value of eta coordinate on the normalized space
 *@param[in] mu Is the value of mu coordinate on the normalized space
 *@return A double value with the value of the derivative
 */
double Hexahedron::DerivateFunction( int id , int dim , double psi , double eta , double mu ){
	assert(  ( id >= 0 ) && ( id < 8 )  );
	assert(  ( dim >= 0 ) && ( dim < 3 )  );
	double value;
	switch( id ){
		case 0:
			if(  dim == 0  ){
				value = ( ( -0.125 ) * ( 1-eta ) * ( 1-mu ) );
			}else{
				if(  dim == 1  ){
					value = ( ( -0.125 ) * ( 1-psi ) * ( 1-mu ) );
				}else{
					value = ( ( -0.125 ) * ( 1-psi ) * ( 1-eta ) );
				}
			}
			break;	
		case 1:
			if(  dim == 0  ){
				value = ( ( 0.125 ) * ( 1-eta ) * ( 1-mu ) );
			}else{
				if(  dim == 1  ){
					value = ( ( -0.125 ) * ( 1+psi ) * ( 1-mu ) );
				}else{
					value = ( ( -0.125 ) * ( 1+psi ) * ( 1-eta ) );
				}
			}
			break;
		case 2:
			if(  dim == 0  ){
				value = ( ( 0.125 ) * ( 1+eta ) * ( 1-mu ) );
			}else{
				if(  dim == 1  ){
					value = ( ( 0.125 ) * ( 1+psi ) * ( 1-mu ) );
				}else{
					value = ( ( -0.125 ) * ( 1+psi ) * ( 1+eta ) );
				}
			}
			break;
		case 3:
			if(  dim == 0  ){
				value = ( ( -0.125 ) * ( 1+eta ) * ( 1-mu ) );
			}else{
				if(  dim == 1  ){
					value = ( ( 0.125 ) * ( 1-psi ) * ( 1-mu ) );
				}else{
					value = ( ( -0.125 ) * ( 1-psi ) * ( 1+eta ) );
				}
			}
			break;
		case 4:
			if(  dim == 0  ){
				value = ( ( -0.125 ) * ( 1-eta ) * ( 1+mu ) );
			}else{
				if(  dim == 1  ){
					value = ( ( -0.125 ) * ( 1-psi ) * ( 1+mu ) );
				}else{
					value = ( ( 0.125 ) * ( 1-psi ) * ( 1-eta ) );
				}
			}
			break;
		case 5:
			if(  dim == 0  ){
				value = ( ( 0.125 ) * ( 1-eta ) * ( 1+mu ) );
			}else{
				if(  dim == 1  ){
					value = ( ( -0.125 ) * ( 1+psi ) * ( 1+mu ) );
				}else{
					value = ( ( 0.125 ) * ( 1+psi ) * ( 1-eta ) );
				}
			}
			break;
		case 6:
			if(  dim == 0  ){
				value = ( ( 0.125 ) * ( 1+eta ) * ( 1+mu ) );
			}else{
				if(  dim == 1  ){
					value = ( ( 0.125 ) * ( 1+psi ) * ( 1+mu ) );
				}else{
					value = ( ( 0.125 ) * ( 1+psi ) * ( 1+eta ) );
				}
			}
			break;
		case 7:
			if(  dim == 0  ){
				value = ( ( -0.125 ) * ( 1+eta ) * ( 1+mu ) );
			}else{
				if(  dim == 1  ){
					value = ( ( 0.125 ) * ( 1-psi ) * ( 1+mu ) );
				}else{
					value = ( ( 0.125 ) * ( 1-psi ) * ( 1+eta ) );
				}
			}
			break;
	}
	return value;
}

/**
 *Getting the jacobian for the element
 *@param[out] jacobiano Is the array where the jacobian matrix will be saved
 */
void Hexahedron::CalculateJacobian( double** jacobiano , int pi ){
	//Inicializing jacobian in zeros
	for(  int i_row = 0  ;  i_row < 3  ;  i_row++  ){
		for(  int i_col = 0  ;  i_col < 3  ;  i_col++  ){
			jacobiano[ i_row ][ i_col ] = 0.0;
		}
	}
	double value =0.5773502692; 
	double psi[ 8 ] = { -value,value,value,-value,-value,value,value,-value };
	double eta[ 8 ] = { -value,-value,value,value,-value,-value,value,value };
	double mu [ 8 ] = { -value,-value,-value,-value,value,value,value,value };
	//begining to fill information
	for(  int i_func = 0  ;  i_func < 8  ;  i_func++  ){
		for(  int i_row = 0  ;  i_row < 3  ;  i_row++  ){
			for(  int i_col = 0  ;  i_col < 3  ;  i_col++  ){
				jacobiano[ i_row ][ i_col ] += 	DerivateFunction( i_func , i_row , psi[ pi ] , eta[ pi ] , mu[ pi ] )*nodes_[ i_func ]->GetCoord( i_col );
			}
		}
	}
}

/**
 *Evaluating a form function givven the normalized coordinates
 *@param[in] id is the index of the form function
 *@param[in] psi Is the psi normalized coordinate
 *@param[in] eta Is the eta normalized coordinate
 *@param[in] mu Is the mu normalized coordinate
 *@return A double value with the evaluation of the form function
 */
double Hexahedron::FormFunction( int id , double psi , double eta , double mu ){
	assert( ( id >= 0 ) && ( id < 8 ) );
	double value;
	switch( id ){
		case 0:
			value = ( 0.125 )*( 1-psi )*( 1-eta )*( 1-mu );
			break;
		case 1:
			value = ( 0.125 )*( 1+psi )*( 1-eta )*( 1-mu );
			break;
		case 2:
			value = ( 0.125 )*( 1+psi )*( 1+eta )*( 1-mu );
			break;
		case 3:
			value = ( 0.125 )*( 1-psi )*( 1+eta )*( 1-mu );
			break;
		case 4:
			value = ( 0.125 )*( 1-psi )*( 1-eta )*( 1+mu );
			break;
		case 5:
			value = ( 0.125 )*( 1+psi )*( 1-eta )*( 1+mu );
			break;
		case 6:
			value = ( 0.125 )*( 1+psi )*( 1+eta )*( 1+mu );
			break;
		case 7:
			value = ( 0.125 )*( 1-psi )*( 1+eta )*( 1+mu );
			break;
	}
	return value;
}

/**
 *Inverting the jacobian matrix
 *@param[in] jac Is the matrix to be inverted
 *@param[out] inv Is the inverted matrix
 */
void Hexahedron::InvertJacobiano( double** jac , double** inv ){
	//Calculo de determinante
	double det = 0.0;
	double jacT[ 3 ][ 3 ];
	det += jac[ 0 ][ 0 ]*( jac[ 1 ][ 1 ] * jac[ 2 ][ 2 ] - jac[ 1 ][ 2 ]*jac[ 2 ][ 1 ] );
	det -= jac[ 0 ][ 1 ]*( jac[ 1 ][ 0 ] * jac[ 2 ][ 2 ] - jac[ 2 ][ 0 ]*jac[ 1 ][ 2 ] );
	det += jac[ 0 ][ 2 ]*( jac[ 1 ][ 0 ] * jac[ 2 ][ 1 ] - jac[ 2 ][ 0 ]*jac[ 1 ][ 1 ] );

	//Matriz transpuesta
	for(  int i_row = 0  ;  i_row < 3  ;  i_row++  ){
		for(  int i_col = 0  ;  i_col < 3  ;  i_col++  ){
			jacT[ i_row ][ i_col ] = jac[ i_col ][ i_row ];
		}
	}
	inv[ 0 ][ 0 ] =  ( ( jacT[ 1 ][ 1 ]*jacT[ 2 ][ 2 ] )-( jacT[ 2 ][ 1 ]*jacT[ 1 ][ 2 ] ) )/det; 
	inv[ 0 ][ 1 ] = -( ( jacT[ 1 ][ 0 ]*jacT[ 2 ][ 2 ] )-( jacT[ 2 ][ 0 ]*jacT[ 1 ][ 2 ] ) )/det;
	inv[ 0 ][ 2 ] =  ( ( jacT[ 1 ][ 0 ]*jacT[ 2 ][ 1 ] )-( jacT[ 2 ][ 0 ]*jacT[ 1 ][ 1 ] ) )/det; 
	inv[ 1 ][ 0 ] = -( ( jacT[ 0 ][ 1 ]*jacT[ 2 ][ 2 ] )-( jacT[ 2 ][ 1 ]*jacT[ 0 ][ 2 ] ) )/det; 
	inv[ 1 ][ 1 ] =  ( ( jacT[ 0 ][ 0 ]*jacT[ 2 ][ 2 ] )-( jacT[ 2 ][ 0 ]*jacT[ 0 ][ 2 ] ) )/det; 
	inv[ 1 ][ 2 ] = -( ( jacT[ 0 ][ 0 ]*jacT[ 2 ][ 1 ] )-( jacT[ 2 ][ 0 ]*jacT[ 0 ][ 1 ] ) )/det; 
	inv[ 2 ][ 0 ] =  ( ( jacT[ 0 ][ 1 ]*jacT[ 1 ][ 2 ] )-( jacT[ 1 ][ 1 ]*jacT[ 0 ][ 2 ] ) )/det; 
	inv[ 2 ][ 1 ] = -( ( jacT[ 0 ][ 0 ]*jacT[ 1 ][ 2 ] )-( jacT[ 1 ][ 0 ]*jacT[ 0 ][ 2 ] ) )/det; 
	inv[ 2 ][ 2 ] =  ( ( jacT[ 0 ][ 0 ]*jacT[ 1 ][ 1 ] )-( jacT[ 1 ][ 0 ]*jacT[ 0 ][ 1 ] ) )/det; 
}

/**
 *This is the evaluation of the interpolation function given a psi, eta and mu
 *f(psi,eta,mu,X) = sum( Ni*Xi ) - X
 *@param[in] psi Is the psi vale of the cordinate normalized
 *@param[in] eta Is the eta vale of the cordinate normalized
 *@param[in] mu Is the mu vale of the cordinate normalized
 *@param[in] dim Is the index that indicates the coordinates to be used in the interpolation 0-X,1-Y, 2-Z
 *@param[in] X Is the coordinate that is being interpolated
 */
double Hexahedron::InterpolationFunction( double psi , double eta , double mu , int dim , double X ){
  return (0.125*((1-psi)*(1-eta)*(1-mu)*nodes_[ 0 ]->GetCoord( dim )+(1+psi)*(1-eta)*(1-mu)*nodes_[ 1 ]->GetCoord( dim )+
								 (1+psi)*(1+eta)*(1-mu)*nodes_[ 2 ]->GetCoord( dim )+(1-psi)*(1+eta)*(1-mu)*nodes_[ 3 ]->GetCoord( dim )+
								 (1-psi)*(1-eta)*(1+mu)*nodes_[ 4 ]->GetCoord( dim )+(1+psi)*(1-eta)*(1+mu)*nodes_[ 5 ]->GetCoord( dim )+
								 (1+psi)*(1+eta)*(1+mu)*nodes_[ 6 ]->GetCoord( dim )+(1-psi)*(1+eta)*(1+mu)*nodes_[ 7 ]->GetCoord( dim ))-X);
}

/**
 *Derivate of the interpolation function respect to psi.
 *df(psi,eta,mu,X)/dpsi = dsum( Ni*Xi ) - X/dpsi
 *@param[in] psi Is the psi value of the cordinate normalized
 *@param[in] eta Is the eta value of the cordinate normalized
 *@param[in] mu Is the mu value of the cordinate normalized
 *@param[in] dim Is the index that indicates the coordinates to be used in the interpolation 0-X,1-Y, 2-Z
 */
double Hexahedron::InterpolationDPsi( double psi , double eta , double mu , int dim ){
	return (0.125*(-(1-eta)*(1-mu)*nodes_[ 0 ]->GetCoord( dim )+(1-eta)*(1-mu)*nodes_[ 1 ]->GetCoord( dim )+
									(1+eta)*(1-mu)*nodes_[ 2 ]->GetCoord( dim )-(1+eta)*(1-mu)*nodes_[ 3 ]->GetCoord( dim )-
                  (1-eta)*(1+mu)*nodes_[ 4 ]->GetCoord( dim )+(1-eta)*(1+mu)*nodes_[ 5 ]->GetCoord( dim )+
									(1+eta)*(1+mu)*nodes_[ 6 ]->GetCoord( dim )-(1+eta)*(1+mu)*nodes_[ 7 ]->GetCoord( dim )) );
}

/**
 *Derivate of the interpolation function respect to eta.
 *df(psi,eta,mu,X)/deta = dsum( Ni*Xi ) - X/deta
 *@param[in] psi Is the psi value of the cordinate normalized
 *@param[in] eta Is the eta value of the cordinate normalized
 *@param[in] mu Is the mu value of the cordinate normalized
 *@param[in] dim Is the index that indicates the coordinates to be used in the interpolation 0-X,1-Y, 2-Z
 */
double Hexahedron::InterpolationDEta( double psi , double eta , double mu , int dim ){
	return ( 0.125*(-(1-psi)*(1-mu)*nodes_[ 0 ]->GetCoord( dim )-(1+psi)*(1-mu)*nodes_[ 1 ]->GetCoord( dim )+
									 (1+psi)*(1-mu)*nodes_[ 2 ]->GetCoord( dim )+(1-psi)*(1-mu)*nodes_[ 3 ]->GetCoord( dim )-
                   (1-psi)*(1+mu)*nodes_[ 4 ]->GetCoord( dim )-(1+psi)*(1+mu)*nodes_[ 5 ]->GetCoord( dim )+
									 (1+psi)*(1+mu)*nodes_[ 6 ]->GetCoord( dim )+(1-psi)*(1+mu)*nodes_[ 7 ]->GetCoord( dim )) );
}

/**
 *Derivate of the interpolation function respect to mu.
 *df(psi,eta,mu,X)/dmu = dsum( Ni*Xi ) - X/dmu
 *@param[in] psi Is the psi value of the cordinate normalized
 *@param[in] eta Is the eta value of the cordinate normalized
 *@param[in] mu Is the mu value of the cordinate normalized
 *@param[in] dim Is the index that indicates the coordinates to be used in the interpolation 0-X,1-Y, 2-Z
 */
double Hexahedron::InterpolationDMu( double psi , double eta , double mu , int dim ){
	return ( 0.125*(-(1-psi)*(1-eta)*nodes_[ 0 ]->GetCoord( dim )-(1+psi)*(1-eta)*nodes_[ 1 ]->GetCoord( dim )-
									 (1+psi)*(1+eta)*nodes_[ 2 ]->GetCoord( dim )+(1-psi)*(1+eta)*nodes_[ 3 ]->GetCoord( dim )+
                   (1-psi)*(1-eta)*nodes_[ 4 ]->GetCoord( dim )+(1+psi)*(1-eta)*nodes_[ 5 ]->GetCoord( dim )+
									 (1+psi)*(1+eta)*nodes_[ 6 ]->GetCoord( dim )+(1-psi)*(1+eta)*nodes_[ 7 ]->GetCoord( dim )) );
}

/**
 *Testing if point is contained int the hexahedron
 *@param[in] coord Is the coordinate to be tested
 *@return A bool value indicating if intersects or not
 */
bool Hexahedron::ContainPoint( double* coord ){
	double nodes[ 8 ][ 3 ];
	for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			nodes[ i_node ][ i_dim ] = nodes_[ i_node ]->GetCoord( i_dim );
		}
	}
	double vec1[ 6 ][ 3 ];
	double vec2[ 6 ][ 3 ];
	double normal[ 6 ][ 3 ];
	double test[ 6 ][ 3 ];
	double PP[ 6 ];
	//Calculating vec1
	CalcVector( vec1[ 0 ] , nodes[ 0 ] , nodes[ 3 ] );
	CalcVector( vec1[ 1 ] , nodes[ 0 ] , nodes[ 4 ] );
	CalcVector( vec1[ 2 ] , nodes[ 1 ] , nodes[ 5 ] );
	CalcVector( vec1[ 3 ] , nodes[ 2 ] , nodes[ 6 ] );
	CalcVector( vec1[ 4 ] , nodes[ 0 ] , nodes[ 1 ] );
	CalcVector( vec1[ 5 ] , nodes[ 4 ] , nodes[ 7 ] );

	//Calculating vec2
	CalcVector( vec2[ 0 ] , nodes[ 0 ] , nodes[ 4 ] );
	CalcVector( vec2[ 1 ] , nodes[ 0 ] , nodes[ 1 ] );
	CalcVector( vec2[ 2 ] , nodes[ 1 ] , nodes[ 2 ] );
	CalcVector( vec2[ 3 ] , nodes[ 2 ] , nodes[ 3 ] );
	CalcVector( vec2[ 4 ] , nodes[ 0 ] , nodes[ 3 ] );
	CalcVector( vec2[ 5 ] , nodes[ 4 ] , nodes[ 5 ] ); 

	//Calculating test vector
	CalcVector( test[ 0 ] , nodes[ 0 ] , coord );
	CalcVector( test[ 1 ] , nodes[ 0 ] , coord );
	CalcVector( test[ 2 ] , nodes[ 1 ] , coord );
	CalcVector( test[ 3 ] , nodes[ 2 ] , coord );
	CalcVector( test[ 4 ] , nodes[ 0 ] , coord );
	CalcVector( test[ 5 ] , nodes[ 4 ] , coord ); 

	//Calculating Normal, normalizing test and normal vector
	int count = 0;
	for(  int i_face = 0  ;  i_face < 6  ;  i_face++  ){
		Cross( normal[ i_face ] , vec1[ i_face ] , vec2[ i_face ] );
		Normalize( normal[ i_face ] );
		Normalize( test[ i_face ] );
		PP[ i_face ] = Dot( normal[ i_face ] , test[ i_face ] );
		if(  PP[ i_face ] >= 0.0  ){
			count++;
		}
	}
	//Counting amount of positive dot products and returning result
	if(  count == 6  ){
		return true;
	}
	return false;
}

/**
 *Calculating the distande from a point outside of an hexahedron to the hexahedron by 
 *using the dot product of the normal to each face and the vector to the coordinate
 *@param[in] coord Is the coordinate from the one that the distance is calculated
 *@return A double value with the distance calculated
 */
double Hexahedron::CalcDistance( double* coord ){
	double nodes[ 8 ][ 3 ];
	for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			nodes[ i_node ][ i_dim ] = nodes_[ i_node ]->GetCoord( i_dim );
		}
	}
	double vec1[ 6 ][ 3 ];
	double vec2[ 6 ][ 3 ];
	double normal[ 6 ][ 3 ];
	double test[ 6 ][ 3 ];
	double PP[ 6 ];
	//Calculating vec1
	CalcVector( vec1[ 0 ] , nodes[ 0 ] , nodes[ 3 ] );
	CalcVector( vec1[ 1 ] , nodes[ 0 ] , nodes[ 4 ] );
	CalcVector( vec1[ 2 ] , nodes[ 1 ] , nodes[ 5 ] );
	CalcVector( vec1[ 3 ] , nodes[ 2 ] , nodes[ 6 ] );
	CalcVector( vec1[ 4 ] , nodes[ 0 ] , nodes[ 1 ] );
	CalcVector( vec1[ 5 ] , nodes[ 4 ] , nodes[ 7 ] );

	//Calculating vec2
	CalcVector( vec2[ 0 ] , nodes[ 0 ] , nodes[ 4 ] );
	CalcVector( vec2[ 1 ] , nodes[ 0 ] , nodes[ 1 ] );
	CalcVector( vec2[ 2 ] , nodes[ 1 ] , nodes[ 2 ] );
	CalcVector( vec2[ 3 ] , nodes[ 2 ] , nodes[ 3 ] );
	CalcVector( vec2[ 4 ] , nodes[ 0 ] , nodes[ 3 ] );
	CalcVector( vec2[ 5 ] , nodes[ 4 ] , nodes[ 5 ] ); 

	//Calculating test vector
	CalcVector( test[ 0 ] , nodes[ 0 ] , coord );
	CalcVector( test[ 1 ] , nodes[ 0 ] , coord );
	CalcVector( test[ 2 ] , nodes[ 1 ] , coord );
	CalcVector( test[ 3 ] , nodes[ 2 ] , coord );
	CalcVector( test[ 4 ] , nodes[ 0 ] , coord );
	CalcVector( test[ 5 ] , nodes[ 4 ] , coord ); 

	//Calculating Normal, normalizing test and normal vector
	double dist = 0.0;
	for(  int i_face = 0  ;  i_face < 6  ;  i_face++  ){
		Cross( normal[ i_face ] , vec1[ i_face ] , vec2[ i_face ] );
		Normalize( normal[ i_face ] );
		Normalize( test[ i_face ] );
		PP[ i_face ] = Dot( normal[ i_face ] , test[ i_face ] );
		if(  PP[ i_face ] < 0.0  ){
			dist += fabs( PP[ i_face ] );
		}
	}
	return dist;
}

//DEBUG

//////////////////////////////////////////////////////////////////////////////////////////
//                                   HEXAHEDRAL METHODS                                 //
//CONSTRUCTOR, DESTRUCTOR
/**
 *Default constructor
 */
Hexahedral::Hexahedral(){
	nodes_.clear();
	elems_.clear();
	scale_ = NULL;
}

/**
 *Constructor receiving the scaler and the name of the file containing the mesh
 *@param[in] scale Is the pointer to the scaler used in to move the mesh to the omain [0,1]
 *@param[in] name Is the name of the file containing the mesh
 */
Hexahedral::Hexahedral( char* name ){
	scale_ = NULL;
	FILE *fin = NULL;
	fin = fopen( name , "r" );
	if( !fin ){
		std::cout << "Mesh file not opened " << std::endl;
		assert( 0 ); 
	}
	while( !feof( fin ) ){
		char aux[100];
		fscanf( fin , "%s" , aux );
		if(  !strcmp( aux , "Coordinates" )  ){
			//Reading nodes
			size_t index = 0;
			while( 1 ){
				double coords[ 3 ];
				fscanf( fin , "%s" , aux );
				if(  !strcmp( aux , "End" )  ){
					break;
				}
				fscanf( fin , "%s" , aux );
				coords[ 0 ] = atof( aux );
				fscanf( fin , "%s" , aux );
				coords[ 1 ] = atof( aux );
				fscanf( fin , "%s" , aux );
				coords[ 2 ] = atof( aux );
				this->AddNode( index , coords );
				index++;
			}
			index = 0;
			fscanf( fin , "%s" , aux );
			fscanf( fin , "%s" , aux );
			//Reading elements
			while( 1 ){
				size_t nodes[ 8 ];
				fscanf( fin , "%s" , aux );
				if(  !strcmp( aux , "End" )  ){
					break;
				}
				for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
					fscanf( fin , "%s" , aux );
					nodes[ i_node ] = atoi( aux );
					nodes[ i_node ]--;	
				}
				this->AddElement( nodes[ 0 ] , nodes[ 1 ] , nodes[ 2 ] , nodes[ 3 ] , nodes[ 4 ] , 
													nodes[ 5 ] , nodes[ 6 ] , nodes[ 7 ]  );
			}
			break;
		}
	}
}

/**
 *Default destructor
 */
Hexahedral::~Hexahedral(){
	nodes_.clear();
	elems_.clear();
	if( scale_ ){
		delete scale_;
	}
}


//GETS
/**
 *Getting number of nodes in the mesh
 *@return A size_t value witht the amount of nodes in the mesh
 */
size_t Hexahedral::GetNNodes(){
	return nodes_.size();
}

/**
 *Getting number of elements on the mesh
 *@return A size_t value with the amount of elements in the mesh
 */
size_t Hexahedral::GetNElements(){
	return elems_.size();
}

/**
 *Getting the pointer of an element in the mesh
 *@param[in] i_elem Is the index of the requested element
 *@return A pointer to an hexahedron
 */
Hexahedron* Hexahedral::GetElement( size_t i_elem ){
	assert( ( i_elem >= 0 )  && ( i_elem < elems_.size() ) );
	return elems_[ i_elem ];
}

/**
 *Getting node from mesh 
 *@param[in] i_node Is the index of the node requested
 *@return A Node pointer 
 */
Node* Hexahedral::GetNode( size_t i_node ){
	assert( ( i_node > 0 )  && ( i_node < nodes_.size() ) );
	return nodes_[ i_node ];
}

/**
 *This function gets one of the coordinates of a node in the mesh
 *@param[in] node Is the index of the node that its coordinate is requested
 *@param[in] dimention Is the index of the coordinate requested (0-X,1-Y,2-Z)
 *@return A double value with the coordinate requested
 */
double Hexahedral::GetCoord( size_t node , int dimention ){
	assert( ( node >= 0 ) && ( node < this->GetNNodes() ) );
	assert( ( dimention >= 0 ) && ( dimention < 3 ) );
	return nodes_[ node ]->GetCoord( dimention );
}

/**
 *Getting the real index of nothe from an element
 *@param[in] element Is the index of the element on the list of elements
 *@param[in] node Is the local labeling of the node requested
 *@return A size_t value with the index requested
 */
size_t Hexahedral::GetIndex( size_t element , int node ){
	assert( ( element >= 0 ) && ( element < this->GetNElements() ) );
	assert( ( node >= 0 ) && ( node < 8 ) );
	return elems_[ element ]->GetIndex( node );
}


//SETS
/**
 *Setting coordinate on one of the nodes in the mesh
 *@param[in] index Is the index of the node where the coordinate will be set
 *@param[in] dimention Is the dimention where the coordinate will be set 0-X,1-Y,2-Z
 *@param[in] coord Is the value to be set
 */
void Hexahedral::SetCoord( size_t index , int dimention , double coord ){
	assert( ( index >= 0 ) && ( index < this->GetNNodes() ) );
	assert( ( dimention >= 0 ) && ( dimention < 3 )  );
	nodes_[ index ]->SetCoord( dimention , coord );
}


//UTILITIES
/**
 *Adding a new node to the mesh
 *@param[in] index Is the node index to be added
 *@param[in] coord Is the vector of coordinates of the node
 */
void Hexahedral::AddNode( size_t index , double* coord ){
	Node* node =  new Node( coord , index );
	nodes_.push_back( node );
}

/**
 *Adding a new element on the mesh
 *@param[in] index0-7 Is the list of indexes of the element
 */
void Hexahedral::AddElement( size_t index0 , size_t index1 , size_t index2 , size_t index3 , 
										 size_t index4 , size_t index5 , size_t index6 , size_t index7 ){
	Hexahedron* hex =  new Hexahedron( nodes_[ index0 ] , nodes_[ index1 ] , nodes_[ index2 ] , 
																		 nodes_[ index3 ] , nodes_[ index4 ] , nodes_[ index5 ] , 
																		 nodes_[ index6 ] , nodes_[ index7 ] );
	elems_.push_back( hex );
}

/**
 *Scaling hexahedral mesh to domain [0,1] 
 */
void Hexahedral::ScaleMesh(  ){
  this->FillScaler( );
	size_t n_nodes = this->GetNNodes(  );
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			this->SetCoord( i_node , i_dim , scale_->ScaleCoord( i_dim , this->GetCoord( i_node , i_dim ) ) );
		} 
	}
}

/**
 *Filling the information on the scaler 
 */
void Hexahedral::FillScaler( ){
	double min_coord[ 3 ] = {  1e20 ,  1e20 ,  1e20 };
	double max_coord[ 3 ] = { -1e20 , -1e20 , -1e20 };
	size_t n_nodes = this->GetNNodes();
	//Obtaining minimum and maximum
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		for(  int i_coord = 0  ;  i_coord < 3  ;  i_coord++  ){
			double coord = this->GetCoord( i_node , i_coord );
			if(  coord < min_coord[ i_coord ]  ){
				min_coord[ i_coord ] = coord;
			}
			if(  coord > max_coord[ i_coord ]  ){
				max_coord[ i_coord ] = coord;
			}
		}
	}
	scale_ = new Scaler( min_coord , max_coord );	
}

/**
 *Unscaling hexahedral mesh
 */
void Hexahedral::UnscaleMesh(  ){


	assert( 0 );//NOT IMPLEMENTED
}

/**
 *This function rotates the hexahedral mesh in order to fit it with the spheres mesh
 */
void Hexahedral::Rotate( ){
	double teta = 1.0187486844;
	double eta = 0.3141592654;
	//double eta = 0.0;
	

	double x_rotation[ 3 ][ 3 ];
	x_rotation[ 0 ][ 0 ] = 1;	x_rotation[ 0 ][ 1 ] = 0         ;	x_rotation[ 0 ][ 2 ] = 0          ;
	x_rotation[ 1 ][ 0 ] = 0;	x_rotation[ 1 ][ 1 ] = cos(teta);	x_rotation[ 1 ][ 2 ] = -sin(teta);
	x_rotation[ 2 ][ 0 ] = 0;	x_rotation[ 2 ][ 1 ] = sin(teta);	x_rotation[ 2 ][ 2 ] = cos(teta) ;

	double y_rotation[ 3 ][ 3 ];
	y_rotation[ 0 ][ 0 ] = cos(eta)  ;	y_rotation[ 0 ][ 1 ] = 0;	y_rotation[ 0 ][ 2 ] = sin(eta)          ;
	y_rotation[ 1 ][ 0 ] = 0          ;	y_rotation[ 1 ][ 1 ] = 1;	y_rotation[ 1 ][ 2 ] = 0;
	y_rotation[ 2 ][ 0 ] = -sin(eta) ;	y_rotation[ 2 ][ 1 ] = 0;	y_rotation[ 2 ][ 2 ] = cos(eta) ;

	size_t n_elements = this->GetNNodes();
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		double coords[ 3 ];
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			coords[ i_dim ] = this->GetCoord( i_elem , i_dim );
		}
		double new_coord[ 3 ];
		for(  int i_row = 0  ;  i_row < 3  ;  i_row++  ){
			new_coord[ i_row ] = 0.0;
			for(  int i_col = 0  ;  i_col < 3  ;  i_col++  ){
				new_coord[ i_row ] += ( y_rotation[ i_row ][ i_col ] * coords[ i_col ] );
			}
		}
		double new_new_coord[ 3 ];
		for(  int i_row = 0  ;  i_row < 3  ;  i_row++  ){
			new_new_coord[ i_row ] = 0.0;
			for(  int i_col = 0  ;  i_col < 3  ;  i_col++  ){
				new_new_coord[ i_row ] += ( x_rotation[ i_row ][ i_col ] * new_coord[ i_col ] );
			}
		}
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			this->SetCoord( i_elem , i_dim , new_new_coord[ i_dim ] );
		}
	}
}
//SAVING MESH ON FILE
/**
 *Saving hexahedral mesh on GiD file
 */
void Hexahedral::SaveMeshOnGiDFile( int mpi_rank , char* output_ ){
	FILE* fp = NULL;
	char name[1000];
	char rank[10];
	sprintf( rank , "%d" , mpi_rank );
	strcpy( name , rank );
	strcat( name , "_hex_" );
	strcat( name , output_ );

	fp = fopen( name , "w" );	
	if(  !fp  ){
		std::cout << " Error, output file not opened " << std::endl;
		std::cout << name << std::endl;
		assert( fp );
	}	

	fprintf( fp , "MESH \"hex\" dimension 3 ElemType Hexahedra Nnode 8\n" );
	fprintf( fp , "# color 96 96 96\n" );
	fprintf( fp , "Coordinates\n" );
	fprintf( fp , "# node_number coordinate_x coordinate_y coordinate_z  \n" );
	size_t n_nodes = this->GetNNodes();
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		fprintf( fp , "%d ", (int)(i_node+1) );
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			double coord = this->GetCoord( (int)i_node , i_dim );
			fprintf( fp , "%lf ", coord );
		}
		fprintf( fp, "\n" ); 
	}

	fprintf( fp , "End coordinates\n" );
	fprintf( fp , "Elements\n" );

	size_t n_elems = this->GetNElements();

	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		fprintf( fp , "%d ", (int)(i_elem+1) );
		for(  int i_node = 0 ;  i_node < 8  ;  i_node++  ){
			size_t index = this->GetIndex( i_elem , i_node ) +1;
			fprintf( fp , "%d ", (int)index );
		}
		fprintf( fp , "\n" );
	}
	fprintf( fp , "End elements\n" );
	fclose( fp );	
}


//DEBUG











