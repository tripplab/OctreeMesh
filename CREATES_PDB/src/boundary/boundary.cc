#include "boundary.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                   SCALER METHODS                                     //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Scaler::Scaler(){
}

/**
 *Constructor receiving min and max coords
 *@param[in] min_coord Is the minimum coordinates of boundary mesh
 *@param[in] max_coord Is the maximum coordinates of boundary mesh
 */
Scaler::Scaler( double* min_coord , double* max_coord ){
	for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
		min_coord_[ i_pos ] = min_coord[ i_pos ];
		max_coord_[ i_pos ] = max_coord[ i_pos ];
		center_[ i_pos ] = ( max_coord_[ i_pos ] + min_coord_[ i_pos ] ) * 0.5;
	}	
	this->SetProportions();
}

/**
 *Default destructor
 */
Scaler::~Scaler(){

}

//GETS
/**
 *Getting min size of model bounding box
 */
double Scaler::GetMinSizeOfModelBoundingBox(){
	double value = 0.0;
	if( proportion_[ 0 ] < proportion_[ 1 ] ){
		if(  proportion_[ 0 ] < proportion_[ 2 ]  ){
			value = proportion_[ 0 ];
		}else{
			value = proportion_[ 2 ];
		}
	}else{
		if(  proportion_[ 1 ] < proportion_[ 2 ]  ){
			value = proportion_[ 1 ];
		}else{
			value = proportion_[ 2 ];
		}
	}
	return value;
}

/**
 *Getting the longest domain of the capside
 *@return A double value with the longest domain of the capside
 */
double Scaler::GetLongestDomain( ){
  double value = -1e20;
  if( domain_[ 0 ] > value   ){
    value = domain_[ 0 ]; 
  }
  if( domain_[ 1 ] > value   ){
    value = domain_[ 1 ]; 
  }
  if( domain_[ 2 ] > value   ){
    value = domain_[ 2 ]; 
  }
	return value;
}

//SETS
/**
 *Set minimum coordinate
 *
 *@param[in] x_min Is the minimum coordinate in the x axis
 *@param[in] y_min Is the minimum coordinate in the y axis
 *@param[in] z_min Is the minimum coordinate in the z axis
 */
void Scaler::SetMinCoord( double x_min , double y_min , double z_min ){
	min_coord_[ 0 ] = x_min;
	min_coord_[ 1 ] = y_min;
	min_coord_[ 2 ] = z_min;
}

/**
 *Set maximum coordinates
 *
 *@param[in] x_max Is the maximum coordinate in the x axis
 *@param[in] y_max Is the maximum coordinate in the y axis
 *@param[in] z_max Is the maximum coordinate in the z axis
 */
void Scaler::SetMaxCoord( double x_max , double y_max , double z_max ){
	max_coord_[ 0 ] = x_max;
	max_coord_[ 1 ] = y_max;
	max_coord_[ 2 ] = z_max;
}

/**
 *Set proportions to scale on every dimension
 *
 *This function establishes the boundary domains in every dimention and calculates a 
 *proportion to scale in order to maintain the topology of the boundary mesh
 */
void Scaler::SetProportions(){
	double max = -1e20;
	int index = -1;
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		domain_[ i_dim ] = max_coord_[ i_dim ] - min_coord_[ i_dim ];
		if(  domain_[ i_dim ] >= max  ){
			max = domain_[ i_dim ]; 
			index = i_dim;
		}
	}
	assert( index >= 0 );
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		if(  i_dim == index  ){
			proportion_[ i_dim ] = 0.99;
		}else{
			proportion_[ i_dim ] = domain_[ i_dim ] / domain_[ index ]*0.99;
		}
	}
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		double size = proportion_[ i_dim ] * 0.5; 
		s_min_coord_[ i_dim ] = 0.5 - size;
		s_max_coord_[ i_dim ] = 0.5 + size;
	}
}

//UTILITIES
/**
 *Coordinate scaler
 *This method takes a coordinate and it is scaled proportionally in order to have the 
 *boundary on the domain [0,1]. new = ( ( ( ( coord - C ) / D ) * P ) + 0.5 )
 *@param[in] i_pos Indicates the domain (X,Y,Z) to use in the scaler
 *@param[in] coord Is the coordinate to be scaled
 *@return The scaled coordinate in a double value
 */
double Scaler::ScaleCoord( int i_pos , double coord ){
	return ( ( ( ( coord - center_[ i_pos ] ) / domain_[ i_pos ] ) * proportion_[ i_pos ] ) + 0.5 );
}

/**
 *Coordinate unscaler
 *This method takes a coordinate and it is unscaled in order to return the coordinate to
 *the original boundary domain new = ( ( ( ( coord - 0.5 ) / P ) * D ) + C )
 *@param[in] i_pos Indicates the domain (X,Y,Z) to use in the unscaler
 *@param[in] coord Is the coordinate to be unscaled
 *@return The coordinate unscaled in a double value
 */
double Scaler::UnscaleCoord( int i_pos , double coord ){
	return ( ( ( ( coord - 0.5 ) / proportion_[ i_pos ] ) * domain_[ i_pos ] ) + center_[ i_pos ] );
}

/**
 *Testing if coordinate is inside domain
 *@param[in] coord Is the coordinate to be tested
 *@param[in] dim Is the axis over the coordinate is tested
 *@return A bool value indicating if intersects or not 
 */
bool Scaler::IsCoordinateInsideDomain( double coord , int dim ){
	if(  coord > s_max_coord_[ dim ] || coord < s_min_coord_[ dim ]  ){
		return false;
	}
	return true;
}

//CLEANERS

//DEBUG
/**
 *Printing information from class scaler
 */
void Scaler::PrintScaler(){
	std::cout << " Printing Scaler Info " << std::endl;
	std::cout << "Min        - " <<  min_coord_[ 0 ] << " " <<  min_coord_[ 1 ] << " " <<  min_coord_[ 2 ] << std::endl;
	std::cout << "Max        - " <<  max_coord_[ 0 ] << " " <<  max_coord_[ 1 ] << " " <<  max_coord_[ 2 ] << std::endl;
	std::cout << "Domain     - " <<     domain_[ 0 ] << " " <<     domain_[ 1 ] << " " <<     domain_[ 2 ] << std::endl;
	std::cout << "Proportion - " << proportion_[ 0 ] << " " << proportion_[ 1 ] << " " << proportion_[ 2 ] << std::endl;
	std::cout << "Center     - " <<     center_[ 0 ] << " " <<     center_[ 1 ] << " " <<     center_[ 2 ] << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                   NODE METHODS                                       //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Node::Node(  ){
}

/**
 *Constructor receiving a pointer to the coordinates 
 *@param[in] coords Is the pointer to the coordinates of the node
 *@param[in] index Id the node index on the mesh
 */
Node::Node( double* coords , size_t index ){
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		coords_[ i_dim ] = coords[ i_dim ];
	}
	index_ = index;
}

/**
 *Constructor receibing three values of the node coordinates
 *@param[in] x_coord Is the x coordinate of the node
 *@param[in] y_coord Is the y coordinate of the node
 *@param[in] z_coord Is the z coordinate of the node
 *@param[in] index Id the node index on the mesh
 */
Node::Node( double x_coord , double y_coord , double z_coord , size_t index ){
	coords_[ 0 ] = x_coord;
	coords_[ 1 ] = y_coord;
	coords_[ 2 ] = z_coord;
  index_ = index;
}

/**
 *Default destructor
 */
Node::~Node(  ){

}

//GETS
/**
 *Get coordinate of node
 *
 *@param[in] i_pos Is the dimension of the coordinate to return
 *@return A double value containing the coordinate
 */
double Node::GetCoord( int i_pos ){
	assert( i_pos >= 0 );
	assert( i_pos <= 2 );
	return coords_[ i_pos ];
}

/**
 *Get node index
 *@return A size_t value that indicates the node index on the mesh
 */
size_t Node::GetIndex(){
	return index_;
}

/**
 *Getting sign from node
 *@return A char value with the node sign
 */
char Node::GetSign(){
	return sign_;
}

//SETS
/**
 *Set coordinate of node
 *
 *@param[in] i_pos Is the position where the coordinate will be located
 *@param[in] val Is the value of the coordinate
 */
void Node::SetCoord( int i_pos , double val ){
	assert( i_pos >= 0 );
	assert( i_pos <= 2 );
	coords_[ i_pos ] = val;
}

/**
 *Set node index on mesh
 *@param[in] index Is the index to set in the node
 */
void Node::SetIndex( size_t index ){
	index_ = index;
}

/**
 *Setting if node is on interfaz or not
 */
void Node::SetInterface( bool interfaz ){
	interfaz_ = interfaz;
}

/**
 *Setting node sign respect to boundary
 */
void Node::SetSign( char sign ){
	sign_ = sign;
}

//UTILITIES
bool Node::IntersectsCell( double* min_coord , double* max_coord ){
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		if(  ( coords_[ i_dim ] < min_coord[ i_dim ] ) || ( coords_[ i_dim ] > max_coord[ i_dim ] )  ){
			return false;
		}
	}
	return true;
}


//////////////////////////////////////////////////////////////////////////////////////////
//                                   TRIANGLE METHODS                                   //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Triangle::Triangle(){

}

/**
 *Constructor receiving the pointer to the nodes
 *
 *@param[in] nod1 Pointer to the first node of the triangle
 *@param[in] nod2 Pointer to the second node of the triangle
 *@param[in] nod3 Pointer to the third node of the triangle
 */
Triangle::Triangle( Node* nod1 , Node* nod2 , Node* nod3 ){
	nodes_[ 0 ] = nod1;
	nodes_[ 1 ] = nod2;
	nodes_[ 2 ] = nod3;
}

/**
 *Default destructor
 */
Triangle::~Triangle(){
}

//GETS
/**
 *Getting node from triangle
 *@param[in] i_node Indicates the node index solicited
 *@return A node pointer to the node solicited
 */
Node* Triangle::GetNode( int i_node ){
	return nodes_[ i_node ];
}

/**
 *Getting coordinate from node
 *@param[in] i_node Is the node used to obtain coordinate
 *@param[in] i_coord Is the coordinate position
 *@return A double value with the coordinate
 */
double Triangle::GetCoord( int i_node , int i_coord ){
	return nodes_[ i_node ]->GetCoord( i_coord );
}

/**
 *Getting all coordinates from a triangle
 *@param[out] T Is the matrix where the coordinates will be stored
 */
void Triangle::GetAllCoordinates( double** T ){
	for(  int i_node = 0  ;  i_node < 3  ;  i_node++  ){
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			T[ i_node ][ i_dim ] = this->GetCoord( i_node , i_dim );
		}
	}
} 

/**
 *Getting all coordinates from triangle
 *@param[out] T0 Is the first node of the triangle
 *@param[out] T1 Is the second node of the triangle
 *@param[out] T2 Is the third node of the triangle
 */
void Triangle::GetAllCoordinates( double* T0 , double* T1 , double* T2 ){
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		T0[ i_dim ] = this->GetCoord( 0 , i_dim );
		T1[ i_dim ] = this->GetCoord( 1 , i_dim );
		T2[ i_dim ] = this->GetCoord( 2 , i_dim );
	}
}

/**
 *Getting global node index
 *This method gets the real node index on the boundary mesh
 *@param[in] node Is the number of node of the triangle to obtain its index
 *@return A int value with the node index 
 */
size_t Triangle::GetNodeIndex(  int node  ){
	return nodes_[ node ]->GetIndex();
}

//SETS
/**
 *Setting node 
 *@param[in] node Is the pointer to the node that will be added
 *@param[in] i_pos Is the node position on triangle that will be added
 */
void Triangle::SetNode( Node* node , int i_pos ){
	nodes_[ i_pos ] = node;
}

//UTILITIES
/**
 *Method to test triangle intersection using kratos
 */
bool Triangle::Intersects( double* bbox_min , double* bbox_max , double tolerance ){
	assert(0);
	return false;
}

/**
 *Intersection between the triangle and box
 *In this method the box is received with a tolerance aplied in the boundaries to avoid 
 *numerical mistakes
 *@param[in] bbox_min Is the lower, left frontal coordinate of the box 
 *@param[in] bbox_max Is the upper, right back coordinate of the box
 *@param[in] tolerance Is the epsilon value during the operations
 */
bool Triangle::IntersectsBox( double* bbox_min , double* bbox_max ){
	//Gettinn data to perform the test
	double box_center[ 3 ], box_half_size[ 3 ];
	for(  int i_dim = 0  ;  i_dim < 3  ; i_dim++  ){
		box_center[ i_dim ] = 0.5 * ( bbox_min[ i_dim ] + bbox_max[ i_dim ] );
		box_half_size[ i_dim ] = 0.5 * ( bbox_max[ i_dim ] - bbox_min[ i_dim ] );
	}
	//getting triangle vertexes
	double* triverts[ 3 ];
	triverts[ 0 ] = new double[ 3 ];
	triverts[ 1 ] = new double[ 3 ];
	triverts[ 2 ] = new double[ 3 ];
	for(  int i_node = 0  ;  i_node < 3  ;  i_node++  ){
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			triverts[ i_node ][ i_dim ] = this->GetCoord( i_node , i_dim );
		}
	}
	//performing the test
	bool intersection = TriangleIntersectsBox( box_center , box_half_size , triverts );	

	delete triverts[ 0 ];	
	delete triverts[ 1 ];	
	delete triverts[ 2 ];	

	return intersection;
}

/**
 *Obtaining bounding box of triangle
 *This method is used to reduce the calculation time when the intersection is tested
 *
 *@param[out] min_corner Is the lower, left, frontal corner of the triangle bounding box 
 *@param[out] max_corner Is the upper, right, back corner of the triangle bounding box 
 */
void Triangle::CalcBoundingBox( double* min_corner , double* max_corner ){
	min_corner[ 0 ] = this->GetCoord( 0 , 0 );		max_corner[ 0 ] = this->GetCoord( 0 , 0 );
	min_corner[ 1 ] = this->GetCoord( 0 , 1 );		max_corner[ 1 ] = this->GetCoord( 0 , 1 );
	min_corner[ 2 ] = this->GetCoord( 0 , 2 );		max_corner[ 2 ] = this->GetCoord( 0 , 2 );
	for(  int i_node = 1  ;  i_node < 3  ;  i_node++  ){
		for(  int i_dim = 0  ;  i_dim < 3  ; i_dim++  ){
			double test_coord = this->GetCoord( i_node , i_dim );
			if(  test_coord < min_corner[ i_dim ]  ){
				min_corner[ i_dim ] = test_coord;
			}
			if(  test_coord > max_corner[ i_dim ]  ){
				max_corner[ i_dim ] = test_coord;
			}
		}
	}
}

/**
 *Testing if triangle intersects a line
 *@param[in] P0 Is the first point of the line
 *@param[in] P1 Is the second point of the line
 *@return a bool value indicating if the line interscts the triangle
 */
bool Triangle::IntersectsLine( double* P0 , double* P1 ){
	double T0[ 3 ], T1[ 3 ], T2[ 3 ];
	for(  int i_coord = 0  ;  i_coord < 3  ;  i_coord++  ){
		T0[ i_coord ] = this->GetCoord( 0 , i_coord );
		T1[ i_coord ] = this->GetCoord( 1 , i_coord );
		T2[ i_coord ] = this->GetCoord( 2 , i_coord );
	}
	return LineIntersectsTriangle( P0 , P1 , T0 , T1 , T2 );
}

/**
 *Getting all intersections of a line with a triangle
 *This method calculates all the intersection between a line and a triangle. If the line 
 *is complanar to the triangle then could exist at most two intersections.
 *Tolerance is stablished because the maximum refinement level allowed is 27, and to
 *determine if a line intersects a triangular face is calculated the signed distance of 
 *the nodes of the line to the plane that contains the triangular face, and given the
 *maximum refinement level, the minimum distance that can exist between the face and the 
 *node is 1/2^27, then obtaining the squared this is 1/2^54, so, the tolerance is 1/2^55
 *because islower than the minimum distance multiplied.  
 *@param[in] P0 Is the initial point of the line
 *@param[in] P1 Is the final point of the line
 *@param[out] X Is the vector where the intersections will be stored
 *@param[in] dim Is the parallel axis to the ray traced
 *@return 0-If there is not intersection
 *@return 1-If there is one intersection
 *@return 2-If there are two intersection
 */
int Triangle::GetAllIntersectionsWithLine( double* P0 , double* P1 , double* X , int dim ){
	long double tol = ( 1.0 / (long double)( 1UL<<55 ) );
	double T0[ 3 ], T1[ 3 ] , T2[ 3 ];
	this->GetAllCoordinates( T0 , T1 , T2 );
	double dist1 , dist2 , val_test , eps ;
	dist1 = SignedDistanceFromFaceToPoint( P0 , T0 , T1 , T2 );
	dist2 = SignedDistanceFromFaceToPoint( P1 , T0 , T1 , T2 );
	val_test = dist1 * dist2;
	if(  ((long double) fabs( val_test ) ) < tol  ){
		//Line is over the plane that contains the triangular face, need to be determined if 
		//line intersects the triangular face and in which points
		if(  !FindStartAndEndOfIntersection( P0 , P1 , T0 , T1 , T2 , X , dim )  ){
			//The line and triangular face are coplanars and in this case the line does not 
			//intersects the triangular face
			return 0;
		}
		return 2;
	}else{
		if(  val_test < 0.0  ){
			//There is one intersection
			eps = ( dist1 ) / ( dist1 - dist2 );
			if(  !FindIntersectionLineTriangle( P0 , P1 , T0 , T1 , T2 , X , eps , dim )  ){
				//Line intersects plane containing triangle but does not intersects the triangle
				return 0;
			}
			return 1;
		}else{
			//Line does not intersects
			return 0;
		}
	}
}

/**
 *testing if a triangle contains a ccordinate
 *@param[in] coord Is the coordinate to be tested
 *@return A bool value indicating if coord is contained on triangle
 */
bool Triangle::ContainsCoord( double* coord ){
	double T0[ 3 ], T1[ 3 ] , T2[ 3 ];
	this->GetAllCoordinates( T0 , T1 , T2 );	
	//testing if plane contains coord
	double AB[ 3 ] , AC[ 3 ] , N[ 3 ] , NN[ 3 ] , D , BAR[ 3 ];
	CalcVector( AB , T0 , T1 );
	CalcVector( AC , T0 , T2 );
	Cross( N , AB , AC );
	NN[ 0 ] = N[ 0 ];
	NN[ 1 ] = N[ 1 ];
	NN[ 2 ] = N[ 2 ];
	Normalize( NN );
	D = Dot( NN , T0 );
	if(  !( fabs( Dot( coord , NN ) + D ) < 1e-10 )  ){
		return false;
	}
	//testing if triangle contains coord
	CalculateBaricentricCoordinates( coord , T0 , T1 , T2 , BAR );
	if(  ( BAR[ 0 ] >= -1e-10 ) && ( BAR[ 1 ] >= -1e-10 ) && ( BAR[ 2 ] >= -1e-10 )  ){
		return true;
	}
	return false;
}

/**
 *Testing if node intersects box using separating axis therorem
 *@param[in] node Is the node to be tested
 *@param[in] min_coord Is the min coordinate of the AABB
 *@param[in] max_coord Is the maximum coordinate of the AABB
 *@return A bool value indicating if node lies inside AABB or not
 */
bool Triangle::NodeIntersectsCell( int node , double* min_coord , double* max_coord ){
	return nodes_[ node ]->IntersectsCell(  min_coord , max_coord  );
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                   BOUNDARY METHODS                                   //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Boundary::Boundary(){
	scale_ = NULL;
}

/**
 *Constructor receiving a boundary filename
 *@param[in] name Is the name of the input file
 */ 
Boundary::Boundary( char* name ){
	std::ifstream file;
	file.open( name );
	if(  !file.is_open()  ){
		std::cout << "File could not be opened---PROGRAM FINISHED" << name << std::endl;
		assert( 0 );
	}
	std::string word;
	while(  file >> word  ){
		if(  word == "Coordinates"  ){
			//Reading nodes information
			file >> word;
			do{
				Node* node;
				double coord[ 3 ];
				for(  size_t i_pos = 0  ;  i_pos < 3  ;  i_pos++){
					file >> word;
					coord[ i_pos ] = std::stod( word );
				}
				node = new Node( coord , GetNNodes() );
				this->SetNode( node );
				file >> word;
			}while(  word != "End" );
			file >> word;
		}
		if(  word == "Elements"  ){
			//reading elements information
			file >> word;
			do{
				Node* node[ 3 ];
				int index;
				for(  size_t i_pos = 0  ;  i_pos < 3  ;  i_pos++){
					file >> word;
					index = stoi( word );
					node[ i_pos ] = this->GetNode( index - 1 );
				}
				Triangle* tri = new Triangle( node[ 0 ] , node[ 1 ] , node[ 2 ] );
				this->SetElement( tri );
				file >> word;
			}while(  word != "End"  );
			file >> word;
		}
	}
	file.close();
	this->CreatesScalerClass();
}

/**
 *Default destructor
 */
Boundary::~Boundary(){
	if( scale_ ){
		delete scale_;
	}
	elems_.clear();
	nodes_.clear();
}


//GETS
/**
 *Getting number of nodes on boundary mesh
 *@return A size_t value indicating the amount of nodes contained on the boundary
 */
size_t Boundary::GetNNodes(){
	return nodes_.size();
}

/**
 *Getting number of triangles on boundary
 *@return A size_t value indicating the amount of triangles contained on the boundary
 */
size_t Boundary::GetNElements(){
	return elems_.size();
}

/**
 *Getting node from boundary
 *@param[in] i_node Is the node index to request
 *@return A Node pointer to the node requested
 */
Node* Boundary::GetNode( int i_node ){
	return nodes_[ i_node ];
}

/**
 *Getting element from boundary mesh
 *@param[in] i_elem Is the index of requested element
 *@return A Triangle pointer to the element
 */
Triangle* Boundary::GetElement( int i_elem ){
	return elems_[ i_elem ];
}

/**
 *Getting coordinate from node on boundary mesh
 *@param[in] i_node Is the node index
 *@param[in] i_coord Is the coordinate index requested
 *@return A double value with the coordinate solicited
 */
double Boundary::GetNodeCoord( size_t i_node , int i_coord ){
	return nodes_[ i_node ]->GetCoord( i_coord );
}

/**
 *Getting minsize of bounding box from model
 *@return A double value with the minimum size of the model boundingbox
 */
double Boundary::GetMinSizeOfModelBoundingBox(){
	return scale_->GetMinSizeOfModelBoundingBox();
}

/**
 *Getting scaler from boundary class
 *@return A pointer to the scaler class contained on boundary class
 */
Scaler* Boundary::GetScaler(){
	return scale_;
}

//SETS
/**
 *Adding a element to the boundary mesh
 *@param[in] elem Is a Tringle pointer containing the new element information
 */
void Boundary::SetElement( Triangle* elem ){
	elems_.push_back( elem );
}

/**
 *Adding a node to boundary mesh
 *@param[in] node Is a Node pointer containing the new node information
 */
void Boundary::SetNode( Node* node ){
	nodes_.push_back( node );
}

//UTILITIES
/**
 *This method calculates minimum ans maximum from boundary mesh and instance a Scaler class
 */
void Boundary::CreatesScalerClass(){
	double min_coord[ 3 ] = {  1e20 ,  1e20 ,  1e20 };
	double max_coord[ 3 ] = { -1e20 , -1e20 , -1e20 };
	size_t n_nodes = this->GetNNodes();
	//Obtaining minimum and maximum
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		for(  int i_coord = 0  ;  i_coord < 3  ;  i_coord++  ){
			double coord;
			coord = this->GetNodeCoord( i_node , i_coord );
			if(  coord < min_coord[ i_coord ]  ){
				min_coord[ i_coord ] = coord;
			}
			if(  coord > max_coord[ i_coord ]  ){
				max_coord[ i_coord ] = coord;
			}
		}
	}
	//Instancin scaler
	scale_ = new Scaler( min_coord , max_coord );
}

/**
 *Scaling boundary nodes
 *This method scale the boundary nodes to the domain [0,1]
 */
void Boundary::ScaleNodes(){
	int n_nodes = (int)this->GetNNodes();
	for(  int i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		Node* node = this->GetNode( i_node );
		for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
			node->SetCoord(  i_pos , scale_->ScaleCoord( i_pos , node->GetCoord( i_pos ) )  );
		}
	}
}

/**
 *Scaling coordinate
 *@param[in] i_pos Is the dimension to scale
 *@param[in] coord Is the coordinate to be scaled
 *@return A double value with the coordinate scaled
 */
double Boundary::ScaleCoordinate( int i_pos , double coord ){
	return scale_->ScaleCoord( i_pos , coord );
}

/**
 *Unscaling coordinate
 *@param[in] i_pos Is the dimension to unscale
 *@param[in] coord Is the coordinate to be unscaled
 *@return A double value with the coordinate unscaled
 */
double Boundary::UnscaleCoordinate( int i_pos , double coord ){
	return scale_->UnscaleCoord( i_pos , coord );
}

/**
 *Testing if coordinate is inside, or outside boundary
 *@param[in] coord Is the coordinate to be tested
 *@return A char value, -1 means that coordinate is inside boundary, 1 means the opposite
 */
char Boundary::IsOuterPoint( double* coord ){
	int cont = 0;
	size_t n_triangles = this->GetNElements();
	double P1[ 3 ] = {-5.0,-5.0,-5.0};
	GenerateOuterOctreePoint( P1 );
	for(  size_t i_tri = 0  ;  i_tri < n_triangles  ;  i_tri++  ){
		Triangle* tri = elems_[ i_tri ];
		if(  tri->ContainsCoord( coord )  ){
			return 0;
		}
		if(  tri->IntersectsLine( coord , P1 )  ){
			cont++;
		}
	}
	//if the number of intersection with the boundary is even, then the point is outside
	if( cont%2 == 0 ){
		return 1;
	}else{	//if the number of intersection with the boundary is odd, then the point is inside
		return -1;
	}	
}

/**
 *Performing test to know the relative position of point with respect to a triangle mesh
 *@param[in] coord Is the point to be tested
 *@return a char value indicating if a point is inside or outside a polyhedron defined
 *by triangles
 */
char Boundary::PointInPolyhedron( double* coord ){
	bool flag = true;
	char result = 0;
	while( flag ){
		char test1 = this->IsOuterPoint( coord );
		char test2 = this->IsOuterPoint( coord );
		char test3 = this->IsOuterPoint( coord );
		char test4 = this->IsOuterPoint( coord );
		char test5 = this->IsOuterPoint( coord );
		if(  (test1==test2) && (test2==test3) && (test3==test4) && (test4==test5)  ){
			result = test1;
			flag = false;
		}
	}
	return result;
}

/**
 *Testing if coordinate is inside domain
 *@param[in] coord Is the coordinate to be tested
 *@param[in] dim Is the axis over the coordinate is tested
 *@return A bool value indicating if intersects or not 
 */
bool Boundary::IsCoordInsideDomain( double coord , int dim ){
	return scale_->IsCoordinateInsideDomain( coord , dim );
}

//DEBUG
void Boundary::PrintMesh(){
	std::cout << " Printing boundary mesh information " << std::endl;
	size_t n_nodes = this->GetNNodes();
	size_t n_elements = this->GetNElements();
	std::cout << "Nodes = " << n_nodes << std::endl;
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		Node* node = this->GetNode( (int)i_node );
		std::cout << i_node << " - " << node->GetCoord( 0 ) << " " <<  node->GetCoord( 1 ) << " " << node->GetCoord( 2 ) << std::endl;
	}
	std::cout << "Elements " << n_elements << std::endl;
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		std::cout << i_elem << " - ";
		Triangle *tri = this->GetElement( (int)i_elem );
		for(  int i_node = 0  ;  i_node < 3  ;  i_node++  ){
			Node* node = tri->GetNode( i_node );
			std::cout << node->GetIndex() << " ";
		}
		std::cout << std::endl;
	}
}

void Boundary::PrintScaler(){
	scale_->PrintScaler();
}











