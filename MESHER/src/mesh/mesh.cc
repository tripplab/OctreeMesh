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
 *@param[in] percentage Is the percentage of the octre to be used in order to reach the 
 *           desired resolution in the mesh
 */
void Spheres::FillScaler( double percentage ){
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
	scale_ = new Scaler( min_coord , max_coord , percentage );	
	
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
















