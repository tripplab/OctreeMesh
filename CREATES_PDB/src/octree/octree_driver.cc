#include "octree_driver.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                          DATAINSIDEOCTREECELL METHODS                                //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
DataInsideOctreeCell::DataInsideOctreeCell(){
	is_mesh_ = -1;
	int_protein_.clear();
	int_aminoacid_.clear();
	int_atoms_.clear();
	bnd_protein_.clear();
	bnd_aminoacid_.clear();
	bnd_atoms_.clear();
	for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
		points_[ i_node ] = NULL;
		r_index_[ i_node ] = -1;
		printed_[ i_node ] = false;
	}

}

/**
 *Default destructor
 */
DataInsideOctreeCell::~DataInsideOctreeCell(){
}

//GETS
/**
 *Getting info to know if cell belongs to mesh
 *@return A char vlue indicating if cell belongs or not to the final mesh
 */
char DataInsideOctreeCell::_GetIsMesh( ){
	return is_mesh_;
}

/**
 *Getting number of atoms in the interior of the cell
 *@return A size_t value with the amount of atoms on cell
 */
size_t DataInsideOctreeCell::_GetNInteriorAtoms( ){
	return int_atoms_.size();
}

/**
 *Getting number of atoms on the boundary of the cell
 *@return A size_t value with the amount of atoms on cell
 */
size_t DataInsideOctreeCell::_GetNBoundaryAtoms( ){
	return bnd_atoms_.size();
}

/**
 *Getting an atom on interior of the cell
 *@param[in] id Is the atom index
 *@return A Atom pointer with the atom requested
 */
Atom* DataInsideOctreeCell::_GetInteriorAtom( size_t id ){
	assert( ( id >= 0 ) && ( id < int_atoms_.size() ) );
	return int_atoms_[ id ];	
}

/**
 *Getting an atom on boundary of the cell
 *@param[in] id Is the atom index
 *@return A Atom pointer with the atom requested
 */
Atom* DataInsideOctreeCell::_GetBoundaryAtom( size_t id ){
	assert( ( id >= 0 ) && ( id < bnd_atoms_.size() ) );
	return bnd_atoms_[ id ];	
}

/**
 *Getting index from an interior protein on the cell
 *@param[in] id Is the id of the protein on the proteins list
 *@return A size_t value with the protein index on the vdb file
 */
size_t DataInsideOctreeCell::_GetInteriorProteinIndex( size_t id ){
	assert( ( id >= 0 ) && ( id < int_protein_.size() ) );
	return int_protein_[ id ];
}

/**
 *Getting index from an boundary protein on the cell
 *@param[in] id Is the id of the protein on the proteins list
 *@return A size_t value with the protein index on the vdb file
 */
size_t DataInsideOctreeCell::_GetBoundaryProteinIndex( size_t id ){
	assert( ( id >= 0 ) && ( id < bnd_protein_.size() ) );
	return bnd_protein_[ id ];
}

/**
 *Getting index from an interior aminoacid on the cell
 *@param[in] id Is the id of the aminoacid on the aminoacid list
 *@return A size_t value with the aminoacid index on the vdb file
 */
size_t DataInsideOctreeCell::_GetInteriorAminoacidIndex( size_t id ){
	assert( ( id >= 0 ) && ( id < int_aminoacid_.size() ) );
	return int_aminoacid_[ id ];
}

/**
 *Getting real index of node
 *@param[in] node Is the local position of node to know its index
 *@return a int value with the real index
 */
int DataInsideOctreeCell::_GetIndex( int node ){
	assert(  ( node >=0 )  &&  ( node < 8 )  );
	return r_index_[ node ];
}

/**
 *Getting the amount of hexahedrons on the cell
 *@return A size_T value with the amount of hexahedrons in the mesh
 */
size_t DataInsideOctreeCell::_GetNHexahedrons( ){
	return hex_list_.size();

}

/**
 *Getting hexahedron contained in cell
 *@param[in] index Is the local index of the hexahedron contained in the cell
 *@return A hexahedronpointer with the hexahedron requested
 */
Hexahedron* DataInsideOctreeCell::_GetHexahedron( size_t index ){
	assert( (index >=0) && (index < hex_list_.size()) );
	return hex_list_[ index ];
}

/**
 *Getting index from an boundary aminoacid on the cell
 *@param[in] id Is the id of the aminoacid on the aminoacid list
 *@return A size_t value with the aminoacid index on the vdb file
 */
size_t DataInsideOctreeCell::_GetBoundaryAminoacidIndex( size_t id ){
	assert( ( id >= 0 ) && ( id < bnd_aminoacid_.size() ) );
	return bnd_aminoacid_[ id ];
}

//SETS
/**
 *Adding an atom to the interior of actual cell.
 *@param[in] atom Is the pointer to the atom added to the cell
 *@param[in] prot_id Is the protein index of the atom added
 *@param[in] prot_id Is the aminoacid index of the atom added
 */
void DataInsideOctreeCell::_AddInteriorAtom( Atom* atom , size_t prot_id , size_t amino_id ){
	is_mesh_ = 1;
	int_protein_.push_back( prot_id );
	int_aminoacid_.push_back( amino_id );
	int_atoms_.push_back( atom );
}

/**
 *Adding an atom to the boundary of actual cell.
 *@param[in] atom Is the pointer to the atom added to the cell
 *@param[in] prot_id Is the protein index of the atom added
 *@param[in] prot_id Is the aminoacid index of the atom added
 */
void DataInsideOctreeCell::_AddBoundaryAtom( Atom* atom , size_t prot_id , size_t amino_id ){
	is_mesh_ = 1;
	bnd_protein_.push_back( prot_id );
	bnd_aminoacid_.push_back( amino_id );
	bnd_atoms_.push_back( atom );
}

/**
 *Setting node on cell
 *@param[in] node Is the node index to be set
 *@param[in] coord Is the coordinate to be set
 */
void DataInsideOctreeCell::_SetNode( int node , double* coord ){
	assert( ( node >= 0) && ( node < 8 ) );

	points_[ node ] = coord;
}

/**
 *Setting real node index
 *@param[in] node Is the node position on the cell
 *@param[in] index Is the real index to be set on the node
 */
void DataInsideOctreeCell::_SetNodeIndex( int node , int index ){
	r_index_[ node ] = index;
}

/**
 *Setting node printeds
 *@param[in] node Is the local index of the node
 */
void DataInsideOctreeCell::_SetPrinted( int node ){
	assert(  ( node >= 0 ) && ( node < 8 )  );
	printed_[ node ] = true;
}

/**
 *Adding an hexahedron on the local_root
 */
void DataInsideOctreeCell::_AddHexahedron( Hexahedron* hex ){
	hex_list_.push_back( hex );
}

//UTILITIES
/**
 *Checking if node has been set
 *@param[in] Is the node that will be known if is set or not
 *@return A bool value indicating if node is set or nor
 */
bool DataInsideOctreeCell::_IsSetNode( int node ){
	return ( points_[ node ] != NULL );

}

/**
 *knowing if node was printed to file or not
 *@param[in] node Is the node index to be tested
 *@return A bool value indicating if node was printed or not
 */
bool DataInsideOctreeCell::_NodeWasPrinted(  int node  ){
	assert(  ( node >=0 )  &&  ( node < 8 )  );
	return printed_[ node ];
}

//DEBUG
/**
 *Printing data cell info
 */
void DataInsideOctreeCell::_PrintInfo(){
	std::cout << "Cell belongs to mesh : " << is_mesh_ << std::endl;
	if(  is_mesh_ == 1  ){
		std::cout << "Number of interior atoms : " << int_atoms_.size() << std::endl;
		std::cout << "Number of boundary atoms : " << bnd_atoms_.size() << std::endl;
		std::cout << "Nodes : " << std::endl;
		for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
			if(  points_[ i_node ]  ){
				std::cout << " X: " << points_[ i_node ][ 0 ] << " Y: " << points_[ i_node ][ 1 ] 
				<< " Z: " << points_[ i_node ][ 2 ] << " id: " <<r_index_[ i_node ] 
				<< " printed: " << printed_[ i_node ] << std::endl;
			}
		}
	}else{
		std::cout << " Cell does not contain any atom " << std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//                          				OCTREEDRIVER METHODS                                //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
OctreeDriver::OctreeDriver(){
  octree_bin = new Octree();
}

/**
 *Default destructor
 */
OctreeDriver::~OctreeDriver(){
  delete octree_bin;
}

//GETS
/**
 *Getting the octree binary
 *@return A octree pointer containing the octree 
 */
Octree* OctreeDriver::GetOctree()const{
	return octree_bin;
}

/**
 *Getting root cell on octree
 *@return A OctreeCell pointer to the root cell of octree
 */
OctreeCell* OctreeDriver::GetRoot(){
	return octree_bin->pGetRoot();
}

/**
 *Getting cell giving a coordinate
 *@param[in] coord
 *@return A OctreeCell pointer to the cell that contains coordinate coord
 */
OctreeCell* OctreeDriver::GetOctreeCell( const double *coord )const{
	return octree_bin->pGetCellNormalized( coord );
}

/**
 *Getting cell giving a key
 *@param[in] keys Is the key pointer to look for a cell
 *@return A OctreeCell pointer thar contains point related to keys
 */ 
OctreeCell* OctreeDriver::GetOctreeCell( key_type* keys )const{
	return octree_bin->pGetCell( keys );
}

/**
 *Getting leaves from octree
 *This method gets the list of leaves on the octree
 *@param[out] all_leaves Is the reference to the pointer containing all the leaves
 *@return A int value with the amount of leaves on the list 
 */
int OctreeDriver::CalcAllLeavesVector(OctreeCell_vector*& all_leaves){
  int num_leaves;
  if( !all_leaves ){
		all_leaves = new OctreeCell_vector();
	}else{
		all_leaves->clear();
	}
  octree_bin->GetAllLeavesVector( *all_leaves );
  num_leaves = (int)( all_leaves->size() );
  return num_leaves;
}

/**
 *Getting leaves inside a bounding box
 *This function calls to the same function on the OctreeDriver class, and gets all leaves
 *that intersect with a given bounding box
 *@param[in] coord1 Is the lower, left and frontal node of the bounding box
 *@param[in] coord2 Is the upper, right and back node of the bounding box
 *@param[out] leaves Is the reference to the leaves vector
 */
void OctreeDriver::GetLeavesInBoundingBox( const double* coord1 , const double* coord2 , OctreeCell_vector& leaves )const{
 	return octree_bin->GetLeavesInBoundingBoxNormalized( coord1 , coord2 , leaves );
}

/**
 *Getting a neighbour from a cell
 *@param[in] cell Is the pointer to the input cell
 *@param[in] idir Is the direction to search for a neighbour
 *@return A OctreeCell pointer to the neighbour of cell (NULL if neighbour does not exist)
 */
OctreeCell* OctreeDriver::GetNeighbour( const OctreeCell* cell , const int idir )const{    
 	OctreeCell* ret = octree_bin->pGetNeighbourCell( cell , idir );
 	return ret;
}

/**
 *Get left neighbour from cell
 *@param[in] cell Is the cell that its neighbour will be look for
 *@return A OctreeCell pointer to the neighbour cell (NULL if neighbour does not exist)
 */
OctreeCell* OctreeDriver::GetLeftNeighbour( const OctreeCell* cell )const{
 	OctreeCell* ret = octree_bin->pGetLeftCell( cell );
 	return ret;
}

/**
 *Get right neighbour from cell
 *@param[in] cell Is the cell that its neighbour will be look for
 *@return A OctreeCell pointer to the neighbour cell (NULL if neighbour does not exist)
 */
OctreeCell* OctreeDriver::GetRightNeighbour( const OctreeCell* cell )const{
 	OctreeCell* ret = octree_bin->pGetRightCell( cell );
 	return ret;      
}

/**
 *Get top neighbour from cell
 *@param[in] cell Is the cell that its neighbour will be look for
 *@return A OctreeCell pointer to the neighbour cell (NULL if neighbour does not exist)
 */
OctreeCell* OctreeDriver::GetTopNeighbour( const OctreeCell* cell )const{
 	OctreeCell* ret = octree_bin->pGetTopCell( cell );
 	return ret;      
}

/**
 *Get bottom neighbour from cell
 *@param[in] cell Is the cell that its neighbour will be look for
 *@return A OctreeCell pointer to the neighbour cell (NULL if neighbour does not exist)
 */
OctreeCell* OctreeDriver::GetBottomNeighbour( const OctreeCell* cell )const{
 	OctreeCell* ret = octree_bin->pGetBottomCell( cell );
 	return ret;     
}

/**
 *Get front neighbour from cell
 *@param[in] cell Is the cell that its neighbour will be look for
 *@return A OctreeCell pointer to the neighbour cell (NULL if neighbour does not exist)
 */
OctreeCell* OctreeDriver::GetFrontNeighbour( const OctreeCell* cell )const{
 	OctreeCell* ret = octree_bin->pGetFrontCell( cell );
 	return ret;
}

/**
 *Get back neighbour from cell
 *@param[in] cell Is the cell that its neighbour will be look for
 *@return A OctreeCell pointer to the neighbour cell (NULL if neighbour does not exist)
 */
OctreeCell* OctreeDriver::GetBackNeighbour( const OctreeCell* cell )const{
 	OctreeCell* ret = octree_bin->pGetBackCell( cell );
 	return ret;
}

/**
 *Getting coordinate from a position on cell
 *@param[in] cell Is the cell pointer in which the coordinate will be calculated
 *@param[in] ipos Indicates the position where the coordinate will be calculated
 *@param[out] coord_point Is the coordinate calculated on ipos
 */
void OctreeDriver::GetCoordOctreePosition( const OctreeCell* cell , int ipos , double* coord_point )const{
  key_type keys[ 3 ];
  cell->GetKey( ipos , keys );
  for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
    coord_point[ i_dim ] = GetCoordinate( keys[ i_dim ] );
  }
}

/**
 *Gettin coordinate given a key
 *@param[out] key Is the key value from which the coordinate is calculated
 *@return A double value with the coordinate solicited
 */
double OctreeDriver::GetCoordinate( Octree::key_type key )const{
	return GetOctree()->GetCoordinateNormalized( key );
}


//SETS


//UTILITIES
/**
 *Check if octree leaves are balanced
 *This method checks if for each leave its neighbours have at most one level of difference
 *in the refinement level
 *@return A bool value indicating if the octree is balanced
 */
bool OctreeDriver::CheckIsBalanced()const{
	return octree_bin->CheckConstrain2To1();
}

/**
 *Balancing octree
 *This method makes octree cells balanced, that is, for each leave its neighbours have at 
 *most one level of difference in the refinement level.
 *@return Zero if the balance operation was performed correctly 
 */
int OctreeDriver::BalanceOctree(){
  octree_bin->Constrain2To1();
  assert( octree_bin->CheckConstrain2To1() );
  return 0;
}

/**
 *Subdividing an OctreeCell
 *This method divides a cell into eight childs
 *@param[in] cell Is the pointer to the cell that will be divided
 *@return Zero indicating that there was not error when subdividing
 */
int OctreeDriver::SubdivideCell( OctreeCell* cell )const{
  assert( cell->GetLevel() > cell->MIN_LEVEL );
  if( !IsCellEmpty( cell ) ){
		//this is executed if exists triangles intersecting cell
    cell->SubdivideCell();
    //TransferTrianglesToSons
  }else{
		//this is executed if there are not triangles intersecting cell
    cell->SubdivideCell();
  }
  return 0;
}

/**
 *Refines octree until reach a size
 *@param[in] size Is the cell size required in the octree
 *@return Zero if the refinement was performed correctly 
 */
int OctreeDriver::RefineOctreeWithSize( const double size ){     
	int fail = octree_bin->RefineWithUniformSizeNormalized( size );     
 	return fail;
}

/**
 *Knowing if cell is empty
 *@param[in] cell Is the pointer to the cell that is tested
 *@return true If cell does not intersects any triangle, false if intersects at less one triangle 
 */
bool OctreeDriver::IsCellEmpty( const OctreeCell* cell )const{
 	if( !cell->pGetObjects() ){
		return true;
	}
 	return !cell->pGetObjects()->size();    
}

/**
 *Calculates the cell size
 *@param[in] cell Is the pointer to the cell that will be calculated its size
 *@return The size of coordinate in a double value 
 */
double OctreeDriver::CalcSize( const OctreeCell* cell )const{
	return GetOctree()->CalcSizeNormalized( cell );
}


//DEBUG
/**
 *Printing information of octree
 *This method prints the amount of leaves on the octree
 */
void OctreeDriver::PrintInfo(){
	octree_bin->PrintInfo();
}

/**
 *Getting if keys are on the maximum on the different dimentions
 *If this is the case then is returned a value indicating the position on it
 *		CASE----X----Y----Z
 *		  0 ----0----0----0 
 *		  1 ----0----0----1
 *		  2 ----0----1----0
 *		  3 ----0----1----1
 *		  4 ----1----0----0
 *		  5 ----1----0----1
 *		  6 ----1----1----0
 *		  7 ----1----1----1
 *@param[in] keys Is the key of the node
 *@return A int value indicating the position on where the node is, or cero if it is inside the octree  
 */
int GetKindOfCaseOnLocalRoot( key_type* root_max_key , key_type* keys ){

	if(  ( keys[ 0 ] != root_max_key[ 0 ]) && ( keys[ 1 ] != root_max_key[ 1 ]) && ( keys[ 2 ] == root_max_key[ 2 ])   )
		return 1;
	if(  ( keys[ 0 ] != root_max_key[ 0 ]) && ( keys[ 1 ] == root_max_key[ 1 ]) && ( keys[ 2 ] != root_max_key[ 2 ])   )
		return 2;
	if(  ( keys[ 0 ] != root_max_key[ 0 ]) && ( keys[ 1 ] == root_max_key[ 1 ]) && ( keys[ 2 ] == root_max_key[ 2 ])   )
		return 3;
	if(  ( keys[ 0 ] == root_max_key[ 0 ]) && ( keys[ 1 ] != root_max_key[ 1 ]) && ( keys[ 2 ] != root_max_key[ 2 ])   )
		return 4;
	if(  ( keys[ 0 ] == root_max_key[ 0 ]) && ( keys[ 1 ] != root_max_key[ 1 ]) && ( keys[ 2 ] == root_max_key[ 2 ])   )
		return 5;
	if(  ( keys[ 0 ] == root_max_key[ 0 ]) && ( keys[ 1 ] == root_max_key[ 1 ]) && ( keys[ 2 ] != root_max_key[ 2 ])   )
		return 6;
	if(  ( keys[ 0 ] == root_max_key[ 0 ]) && ( keys[ 1 ] == root_max_key[ 1 ]) && ( keys[ 2 ] == root_max_key[ 2 ])   )
		return 7;
	return 0;
}

/**
 *This function substract 1 from the keys depending on the kind of case that is the key
 */
void PerturbateKeysDependingOnKind( key_type* keys , int kind ){
	switch( kind ){
		case 0:
			break;
		case 1:
			keys[ 2 ] -= 1;
			break;
		case 2:
			keys[ 1 ] -= 1;
			break;
		case 3:
			keys[ 1 ] -= 1;
			keys[ 2 ] -= 1;
			break;
		case 4:
			keys[ 0 ] -= 1;
			break;
		case 5:
			keys[ 0 ] -= 1;
			keys[ 2 ] -= 1;
			break;
		case 6:
			keys[ 0 ] -= 1;
			keys[ 1 ] -= 1;
			break;
		case 7:
			keys[ 0 ] -= 1;
			keys[ 1 ] -= 1;
			keys[ 2 ] -= 1;
			break;
	}
}

void UnPerturbateKeysDependingOnKind( key_type* keys , int kind ){
	switch( kind ){
		case 0:
			break;
		case 1:
			keys[ 2 ] += 1;
			break;
		case 2:
			keys[ 1 ] += 1;
			break;
		case 3:
			keys[ 1 ] += 1;
			keys[ 2 ] += 1;
			break;
		case 4:
			keys[ 0 ] += 1;
			break;
		case 5:
			keys[ 0 ] += 1;
			keys[ 2 ] += 1;
			break;
		case 6:
			keys[ 0 ] += 1;
			keys[ 1 ] += 1;
			break;
		case 7:
			keys[ 0 ] += 1;
			keys[ 1 ] += 1;
			keys[ 2 ] += 1;
			break;
	}
}

void GetNodesToAddByRayCasting( int position , int dim , int* nodes ){
	switch( dim ){
		case 0:
			GetNodesToAddByRayCastingX( position , nodes );
			break;
		case 1:
			GetNodesToAddByRayCastingY( position , nodes );
			break;
		case 2:
			GetNodesToAddByRayCastingZ( position , nodes );
			break;
		default:
			assert( 0 );
			break;
	}
}

void GetNodesToAddByRayCastingX( int position , int* nodes ){
	switch(position){
		case 0:
			nodes[ 0 ] = 0;
			nodes[ 1 ] = 9;
			nodes[ 2 ] = 1;
			break;
		case 3:
			nodes[ 0 ] = 3;
			nodes[ 1 ] = 11;
			nodes[ 2 ] = 2;
			break;
		case 4:
			nodes[ 0 ] = 4;
			nodes[ 1 ] = 17;
			nodes[ 2 ] = 5;
			break;
		case 7:
			nodes[ 0 ] = 7;
			nodes[ 1 ] = 19;
			nodes[ 2 ] = 6;
			break;
		case 12:
			nodes[ 0 ] = 12;
			nodes[ 1 ] = 21;
			nodes[ 2 ] = 10;
			break;
		case 13:
			nodes[ 0 ] = 13;
			nodes[ 1 ] = 22;
			nodes[ 2 ] = 14;
			break;
		case 16:
			nodes[ 0 ] = 16;
			nodes[ 1 ] = 24;
			nodes[ 2 ] = 15;
			break;
		case 20:
			nodes[ 0 ] = 20;
			nodes[ 1 ] = 26;
			nodes[ 2 ] = 18;
			break;
		case 25:
			nodes[ 0 ] = 25;
			nodes[ 1 ] = 8;
			nodes[ 2 ] = 23;
			break;
		default:
			assert( 0 );
			break;
	}
}


void GetNodesToAddByRayCastingY( int position , int* nodes ){
	switch(position){
		case 0:
			nodes[ 0 ] = 0;
			nodes[ 1 ] = 12;
			nodes[ 2 ] = 3;
			break;
		case 1:
			nodes[ 0 ] = 1;
			nodes[ 1 ] = 10;
			nodes[ 2 ] = 2;
			break;
		case 4:
			nodes[ 0 ] = 4;
			nodes[ 1 ] = 20;
			nodes[ 2 ] = 7;
			break;
		case 5:
			nodes[ 0 ] = 5;
			nodes[ 1 ] = 18;
			nodes[ 2 ] = 6;
			break;
		case 9:
			nodes[ 0 ] = 9;
			nodes[ 1 ] = 21;
			nodes[ 2 ] = 11;
			break;
		case 13:
			nodes[ 0 ] = 13;
			nodes[ 1 ] = 25;
			nodes[ 2 ] = 16;
			break;
		case 14:
			nodes[ 0 ] = 14;
			nodes[ 1 ] = 23;
			nodes[ 2 ] = 15;
			break;
		case 17:
			nodes[ 0 ] = 17;
			nodes[ 1 ] = 26;
			nodes[ 2 ] = 19;
			break;
		case 22:
			nodes[ 0 ] = 22;
			nodes[ 1 ] = 8;
			nodes[ 2 ] = 24;
			break;
		default:
			assert( 0 );
			break;
	}
}


void GetNodesToAddByRayCastingZ( int position , int* nodes ){
	switch(position){
		case 0:
			nodes[ 0 ] = 0;
			nodes[ 1 ] = 13;
			nodes[ 2 ] = 4;
			break;
		case 1:
			nodes[ 0 ] = 1;
			nodes[ 1 ] = 14;
			nodes[ 2 ] = 5;
			break;
		case 2:
			nodes[ 0 ] = 2;
			nodes[ 1 ] = 15;
			nodes[ 2 ] = 6;
			break;
		case 3:
			nodes[ 0 ] = 3;
			nodes[ 1 ] = 16;
			nodes[ 2 ] = 7;
			break;
		case 9:
			nodes[ 0 ] = 9;
			nodes[ 1 ] = 22;
			nodes[ 2 ] = 17;
			break;
		case 10:
			nodes[ 0 ] = 10;
			nodes[ 1 ] = 23;
			nodes[ 2 ] = 18;
			break;
		case 11:
			nodes[ 0 ] = 11;
			nodes[ 1 ] = 24;
			nodes[ 2 ] = 19;
			break;
		case 12:
			nodes[ 0 ] = 12;
			nodes[ 1 ] = 25;
			nodes[ 2 ] = 20;
			break;
		case 21:
			nodes[ 0 ] = 21;
			nodes[ 1 ] = 8;
			nodes[ 2 ] = 26;
			break;
		default:
			assert( 0 );
			break;
	}
}

int GetPositionOverRay( std::vector<double> points , key_type key ){

	//se determina la posicion de 'key' sobre el rayo, regresa -1 si es un punto
	//contenido en 'points'
	int flag = -1;
	int n_points = (int)points.size();
	double coord = ( (double)key ) / ( (double)(1<<ROOT_LEVEL) );
	//Se verifica si el key es un punto de interseccion
	for(  int i_point = 0  ;  i_point < n_points  ;  i_point++  ){
		if(  points[ i_point ] == coord  ){
			return flag;
		}
	}
	//Si no es un nodo entonces se verifica dentro de que rango esta
	for(  int i_point = 0  ;  i_point < ( n_points - 1 )  ;  i_point++  ){
		if(  ( coord > points[ i_point ] ) && ( coord < points[ i_point + 1 ] )  ){
			flag = i_point;
			break;
		}
	}
	return flag;
}

/**
 *Determining if there are two equal signs and if it hapens, then the node has to be set
 *with that sign
 *@param[in] signs Is the list of three signs (X,Y,Z)
 *@param[out] real_sign Is the value where the real sign will be set if it is possible
 *@return A bool value indicating if node sign can be set or not.
 */
bool NodeSignCanBeSet( int* signs , int* real_sign ){
	int cont_1 = 0;
	int cont_0 = 0;
	int cont1  = 0;
	int contN  = 0;
	switch( signs[ 0 ] ){
		case -1:
			cont_1++;
			break;
		case 0:
			cont_0++;
			break;
		case 1:
			cont1++;
			break;
		case -5:
			contN++;
			break;
	}

	switch( signs[ 1 ] ){
		case -1:
			cont_1++;
			break;
		case 0:
			cont_0++;
			break;
		case 1:
			cont1++;
			break;
		case -5:
			contN++;
			break;
	}

	switch( signs[ 2 ] ){
		case -1:
			cont_1++;
			break;
		case 0:
			cont_0++;
			break;
		case 1:
			cont1++;
			break;
		case -5:
			contN++;
			break;
	}
	if(  cont_1 == 2  ){
		(*real_sign) = -1;
		return true;
	}
	if(  cont_0 == 2  ){
		(*real_sign) = 0;
		return true;
	}
	if(  cont1 == 2  ){
		(*real_sign) = 1;
		return true;
	}
	if(  contN == 2  ){
		if(  cont_1 == 1  ){
			(*real_sign) = -1;
			return true;
		}
		if(  cont_0 == 1  ){
			(*real_sign) = 0;
			return true;
		}
		if(  cont1 == 1  ){
			(*real_sign) = 1;
			return true;
		}
		assert( 0 );
	}
	return false;
}

/**
 *Getting index of center node from a face
 *@param[in] src_dir Is the source direction from where the index is requested
 *@return A int value with the index requested 
 */
int GetCenterIndexFromFace( int src_dir ){

	switch( src_dir ){
		case 0:
			return 23;
			break;
		case 1:
			return 25;
			break;
		case 2:
			return 24;
			break;
		case 3:
			return 22;
			break;
		case 4:
			return 26;
			break;
		case 5:
			return 21;
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting all linear indexes from a cell face
 *@param[in] src_dir Is the source direction from the face
 *@param[out] index Is the array where the indexes will be setted
 */
void GetLinearIndexFromFace( int src_dir , int* index ){
	switch( src_dir ){
		case 0:
			index[ 0 ] = 1;
			index[ 1 ] = 2;
			index[ 2 ] = 5;
			index[ 3 ] = 6;
			break;
		case 1:
			index[ 0 ] = 0;
			index[ 1 ] = 3;
			index[ 2 ] = 4;
			index[ 3 ] = 7;
			break;
		case 2:
			index[ 0 ] = 3;
			index[ 1 ] = 2;
			index[ 2 ] = 6;
			index[ 3 ] = 7;
			break;
		case 3:
			index[ 0 ] = 0;
			index[ 1 ] = 1;
			index[ 2 ] = 4;
			index[ 3 ] = 5;
			break;
		case 4:
			index[ 0 ] = 4;
			index[ 1 ] = 5;
			index[ 2 ] = 6;
			index[ 3 ] = 7;
			break;
		case 5:
			index[ 0 ] = 0;
			index[ 1 ] = 1;
			index[ 2 ] = 2;
			index[ 3 ] = 3;
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting all quadratic indexes from a cell face
 *@param[in] src_dir Is the source direction from the face
 *@param[out] index Is the array where the indexes will be setted
 */
void GetQuadraticIndexFromFace(int src_dir , int* index ){
	switch( src_dir ){
		case 0:
			index[ 0 ] = 10;
			index[ 1 ] = 14;
			index[ 2 ] = 15;
			index[ 3 ] = 18;
			break;
		case 1:
			index[ 0 ] = 12;
			index[ 1 ] = 13;
			index[ 2 ] = 16;
			index[ 3 ] = 20;
			break;
		case 2:
			index[ 0 ] = 11;
			index[ 1 ] = 15;
			index[ 2 ] = 16;
			index[ 3 ] = 19;
			break;
		case 3:
			index[ 0 ] = 9;
			index[ 1 ] = 13;
			index[ 2 ] = 14;
			index[ 3 ] = 17;
			break;
		case 4:
			index[ 0 ] = 17;
			index[ 1 ] = 18;
			index[ 2 ] = 19;
			index[ 3 ] = 20;
			break;
		case 5:
			index[ 0 ] = 9;
			index[ 1 ] = 10;
			index[ 2 ] = 11;
			index[ 3 ] = 12;
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting keys from all neighbours by face when the neighbour is more refined
 *@param[in] direction Is the direcion from where the neighbours are searched
 *@param[in] keys0 Is the neighbour key calculated
 *@param[out] keys1 Is the 2 neighbour key from face
 *@param[out] keys2 Is the 3 neighbour key from face
 *@param[out] keys3 Is the 4 neighbour key from face
 *@param[in] lvel Is the refinement level on face
 */
void GetAllNeighbourKeysFromFaceLessRefined( int direction , key_type* keys0 , key_type* keys1 , 
																						 key_type* keys2 , key_type* keys3 , int level ){
	switch( direction ){
		case 0:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys1[ 2 ] = keys0[ 2 ];
			keys2[ 0 ] = keys0[ 0 ];
			keys2[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys2[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			keys3[ 0 ] = keys0[ 0 ];
			keys3[ 1 ] = keys0[ 1 ];
			keys3[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			break;
		case 1:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys1[ 2 ] = keys0[ 2 ];
			keys2[ 0 ] = keys0[ 0 ];
			keys2[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys2[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			keys3[ 0 ] = keys0[ 0 ];
			keys3[ 1 ] = keys0[ 1 ];
			keys3[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			break;
		case 2:
			keys1[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys1[ 1 ] = keys0[ 1 ];
			keys1[ 2 ] = keys0[ 2 ];
			keys2[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys2[ 1 ] = keys0[ 1 ];
			keys2[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			keys3[ 0 ] = keys0[ 0 ];
			keys3[ 1 ] = keys0[ 1 ];
			keys3[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			break;
		case 3:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ];
			keys1[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			keys2[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys2[ 1 ] = keys0[ 1 ];
			keys2[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			keys3[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys3[ 1 ] = keys0[ 1 ];
			keys3[ 2 ] = keys0[ 2 ];
			break;
		case 4:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys1[ 2 ] = keys0[ 2 ];
			keys2[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys2[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys2[ 2 ] = keys0[ 2 ];
			keys3[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys3[ 1 ] = keys0[ 1 ];
			keys3[ 2 ] = keys0[ 2 ];
			break;
		case 5:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys1[ 2 ] = keys0[ 2 ];
			keys2[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys2[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys2[ 2 ] = keys0[ 2 ];
			keys3[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys3[ 1 ] = keys0[ 1 ];
			keys3[ 2 ] = keys0[ 2 ];
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting linear indexes from edge
 *@param[in] src_dir Is the source direction from the face
 *@param[out] index Is the array where the indexes will be setted
 */
void GetLinearIndexFromEdge( int src_dir , int* index ){
	switch( src_dir ){
		case 6:
			index[ 0 ] = 2;
			index[ 1 ] = 6;
			break;
		case 7:
			index[ 0 ] = 3;
			index[ 1 ] = 7;
			break;
		case 8:
			index[ 0 ] = 1;
			index[ 1 ] = 5;
			break;
		case 9:
			index[ 0 ] = 0;
			index[ 1 ] = 4;
			break;
		case 10:
			index[ 0 ] = 5;
			index[ 1 ] = 6;
			break;
		case 11:
			index[ 0 ] = 4;
			index[ 1 ] = 7;
			break;
		case 12:
			index[ 0 ] = 1;
			index[ 1 ] = 2;
			break;
		case 13:
			index[ 0 ] = 0;
			index[ 1 ] = 3;
			break;
		case 14:
			index[ 0 ] = 6;
			index[ 1 ] = 7;
			break;
		case 15:
			index[ 0 ] = 4;
			index[ 1 ] = 5;
			break;
		case 16:
			index[ 0 ] = 2;
			index[ 1 ] = 3;
			break;
		case 17:
			index[ 0 ] = 0;
			index[ 1 ] = 1;
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting center index from edge
 *@param[in] src_dir Is the source direction from the face
 */
int GetCenterIndexFromEdge( int src_dir ){
	switch( src_dir ){
		case 6:
			return 15;
			break;
		case 7:
			return 16;
			break;
		case 8:
			return 14;
			break;
		case 9:
			return 13;
			break;
		case 10:
			return 18;
			break;
		case 11:
			return 20;
			break;
		case 12:
			return 10;
			break;
		case 13:
			return 12;
			break;
		case 14:
			return 19;
			break;
		case 15:
			return 17;
			break;
		case 16:
			return 11;
			break;
		case 17:
			return 9;
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting keys from all neighbours by edge when the neighbour is more refined
 *@param[in] direction Is the direcion from where the neighbours are searched
 *@param[in] keys0 Is the neighbour key calculated
 *@param[out] keys1 Is the 2 neighbour key from face
 *@param[in] lvel Is the refinement level on local cell
 */
void GetAllNeighbourKeysFromEdgeLessRefined( int direction , key_type* keys0 , key_type* keys1 , 
																						 int level ){
	switch( direction ){
		case 6:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ];
			keys1[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			break;
		case 7:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ];
			keys1[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			break;
		case 8:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ];
			keys1[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			break;
		case 9:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ];
			keys1[ 2 ] = keys0[ 2 ] + ( 1<<(level-1) );
			break;
		case 10:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys1[ 2 ] = keys0[ 2 ];
			break;
		case 11:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys1[ 2 ] = keys0[ 2 ];
			break;
		case 12:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys1[ 2 ] = keys0[ 2 ];
			break;
		case 13:
			keys1[ 0 ] = keys0[ 0 ];
			keys1[ 1 ] = keys0[ 1 ] + ( 1<<(level-1) );
			keys1[ 2 ] = keys0[ 2 ];
			break;
		case 14:
			keys1[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys1[ 1 ] = keys0[ 1 ];
			keys1[ 2 ] = keys0[ 2 ];
			break;
		case 15:
			keys1[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys1[ 1 ] = keys0[ 1 ];
			keys1[ 2 ] = keys0[ 2 ];
			break;
		case 16:
			keys1[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys1[ 1 ] = keys0[ 1 ];
			keys1[ 2 ] = keys0[ 2 ];
			break;
		case 17:
			keys1[ 0 ] = keys0[ 0 ] + ( 1<<(level-1) );
			keys1[ 1 ] = keys0[ 1 ];
			keys1[ 2 ] = keys0[ 2 ];
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting linear and center node from face
 *@param[in] direction Is the direction face to obtain indexes
 *@param[out] index Is the list of indexes requested
 */
void GetLinearAndCenterIndexFromBoundaryFace( int direction , int* index ){
	switch( direction ){
		case 0:
			index[ 0 ] = 0;
			index[ 1 ] = 3;
			index[ 2 ] = 4;
			index[ 3 ] = 7;
			index[ 4 ] = 25;
			break;
		case 1:
			index[ 0 ] = 1;
			index[ 1 ] = 2;
			index[ 2 ] = 5;
			index[ 3 ] = 6;
			index[ 4 ] = 23;
			break;
		case 2:
			index[ 0 ] = 0;
			index[ 1 ] = 1;
			index[ 2 ] = 4;
			index[ 3 ] = 5;
			index[ 4 ] = 22;
			break;
		case 3:
			index[ 0 ] = 2;
			index[ 1 ] = 3;
			index[ 2 ] = 6;
			index[ 3 ] = 7;
			index[ 4 ] = 24;
			break;
		case 4:
			index[ 0 ] = 0;
			index[ 1 ] = 1;
			index[ 2 ] = 2;
			index[ 3 ] = 3;
			index[ 4 ] = 21;
			break;
		case 5:
			index[ 0 ] = 4;
			index[ 1 ] = 5;
			index[ 2 ] = 6;
			index[ 3 ] = 7;
			index[ 4 ] = 26;
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting linear indexes from face
 *@param[in] direction Is the direction to obtain the indexes
 *@param[out] index Is the list of linear indexes on face
 */
void GetLinearIndexOnFace( int direction , int* index ){
	switch( direction ){
		case 0:
			index[ 0 ] = 0;
			index[ 1 ] = 3;
			index[ 2 ] = 4;
			index[ 3 ] = 7;
			break;
		case 1:
			index[ 0 ] = 1;
			index[ 1 ] = 2;
			index[ 2 ] = 5;
			index[ 3 ] = 6;
			break;
		case 2:
			index[ 0 ] = 0;
			index[ 1 ] = 1;
			index[ 2 ] = 4;
			index[ 3 ] = 5;
			break;
		case 3:
			index[ 0 ] = 2;
			index[ 1 ] = 3;
			index[ 2 ] = 6;
			index[ 3 ] = 7;
			break;
		case 4:
			index[ 0 ] = 0;
			index[ 1 ] = 1;
			index[ 2 ] = 2;
			index[ 3 ] = 3;
			break;
		case 5:
			index[ 0 ] = 4;
			index[ 1 ] = 5;
			index[ 2 ] = 6;
			index[ 3 ] = 7;
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting edge indexes to be tested
 *@param[in] direction is the face direction to be searched
 *@param[out] index Is the list of edge indexes that form the cell face
 */
void GetEdgesIndexesFromFace( int direction , int* index ){
	switch( direction ){
		case 0:
			index[ 0 ] = 6;
			index[ 1 ] = 8;
			index[ 2 ] = 10;
			index[ 3 ] = 12;
			break;
		case 1:
			index[ 0 ] = 7;
			index[ 1 ] = 9;
			index[ 2 ] = 11;
			index[ 3 ] = 13;
			break;
		case 2:
			index[ 0 ] = 6;
			index[ 1 ] = 7;
			index[ 2 ] = 14;
			index[ 3 ] = 16;
			break;
		case 3:
			index[ 0 ] = 8;
			index[ 1 ] = 9;
			index[ 2 ] = 15;
			index[ 3 ] = 17;
			break;
		case 4:
			index[ 0 ] = 10;
			index[ 1 ] = 11;
			index[ 2 ] = 14;
			index[ 3 ] = 15;
			break;
		case 5:
			index[ 0 ] = 12;
			index[ 1 ] = 13;
			index[ 2 ] = 16;
			index[ 3 ] = 17;
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting center index on edge
 *@param[in] direction Is the direction where node index is located
 *@return A int value with the index on local cell
 */
int GetCenterIndexOnEdge( int direction ){
	switch( direction ){
		case 6:
			return 13;
			break;
		case 7:
			return 14;
			break;
		case 8:
			return 16;
			break;
		case 9:
			return 15;
			break;
		case 10:
			return 12;
			break;
		case 11:
			return 10;
			break;
		case 12:
			return 20;
			break;
		case 13:
			return 18;
			break;
		case 14:
			return 9;
			break;
		case 15:
			return 11;
			break;
		case 16:
			return 17;
			break;
		case 17:
			return 19;
			break;
		default:
			assert( 0 );
			break;
	}
}

/**
 *Getting relative position of src_cell over cell when the search direction from src_cell 
 *to cell is direction
 *@param[in] src_cell Is the cell less refined
 *@param[in] cell Is the neighbour cell (more refined) of src_cell
 *@param[in] direction Is the search direction from src_cell
 *@return A int value indicating the relative position of src_cell over cell (0,1,2 or 3)
 */
int GetRelativePositionOverNeighbourLessRefined( OctreeCell* src_cell , OctreeCell* cell , int direction ){
	int position;
	size_t src_index;
	size_t dest_index[ 4 ];
	switch(  direction  ){
		case 0:
			src_index = 0;
			dest_index[ 0 ] = 10;
			dest_index[ 1 ] = 23;
			dest_index[ 2 ] = 14;
			dest_index[ 3 ] = 1;
			break;
		case 1:
			src_index = 1;
			dest_index[ 0 ] = 0;
			dest_index[ 1 ] = 13;
			dest_index[ 2 ] = 25;
			dest_index[ 3 ] = 12;
			break;
		case 2:
			src_index = 0;
			dest_index[ 0 ] = 24;
			dest_index[ 1 ] = 11;
			dest_index[ 2 ] = 3;
			dest_index[ 3 ] = 16;
			break;
		case 3:
			src_index = 3;
			dest_index[ 0 ] = 9;
			dest_index[ 1 ] = 22;
			dest_index[ 2 ] = 13;
			dest_index[ 3 ] = 0;
			break;
		case 4:
			src_index = 0;
			dest_index[ 0 ] = 17;
			dest_index[ 1 ] = 26;
			dest_index[ 2 ] = 20;
			dest_index[ 3 ] = 4;
			break;
		case 5:
			src_index = 7;
			dest_index[ 0 ] = 11;
			dest_index[ 1 ] = 21;
			dest_index[ 2 ] = 12;
			dest_index[ 3 ] = 3;
			break;
		default:
			assert( 0 );
			break;
	}
	bool flag = false;
	key_type src_key[ 3 ];
	src_cell->GetKey( src_index , src_key );
	for(  int i_test = 0  ;  i_test < 4  ;  i_test++  ){
		key_type dest_key[ 3 ];
		cell->GetKey( dest_index[ i_test ] , dest_key );
		if(  AreEquals(  src_key , dest_key  )  ){
			position = i_test;
			flag = true;
			break;
		}
	}
	if(  !flag  ){
		std::cout<<"Error: Not founded local position on direction "<<direction<<std::endl;
		src_cell->PrintInfo();
		std::cout<<std::endl;
		cell->PrintInfo();
		assert( 0 );
	} 
  return position;
}

/**
 *Testing if two keys are equals or not
 */
bool AreEquals(  key_type* key0 , key_type* key1  ){
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		if(   key0[ i_dim ] - key1[ i_dim ]   ){
			return false;
		}
	}
	return true;
}


/**
 *Getting tetrahedra connectivities on local rank cell when neighbour from other rank is 
 *more refined.
 *@param[in] position Is the relative position of neighbour face on local cell face
 *@param[in] src_dir Is the direction from where the position was searched
 *@param[out] tet_index Is the array where the tetrahedra nodes on local cell will be set 
 */
void GetLocalTetIndexWhenNeighbourIsMoreRefined( int position , int src_dir , int* tet_index ){
	int indexes[ 6 ][ 4 ][ 8 ] = { {{ 2,23,15, 8, 2,10,23, 8},{ 6,15,23, 8, 6,23,18, 8},{ 5,18,23, 8, 5,23,14, 8},{ 1,14,23, 8, 1,23,10, 8}}, 
																 {{ 0,12,25, 8, 0,25,13, 8},{ 4,13,25, 8, 4,25,20, 8},{ 7,20,25, 8, 7,25,16, 8},{16,25, 3, 8, 3,25,12, 8}}, 
																 {{ 6,19,24, 8, 6,24,15, 8},{15,24, 2, 8,24,11, 2, 8},{11,24, 3, 8, 3,24,16, 8},{ 7,16,24, 8, 7,24,19, 8}}, 
																 {{ 9,22, 1, 8, 1,22,14, 8},{14,22, 5, 8, 5,22,17, 8},{ 4,17,22, 8, 4,22,13, 8},{13,22, 0, 8, 0,22, 9, 8}}, 
																 {{26,18, 5, 8, 5,17,26, 8},{26, 6,18, 8,26,19, 6, 8},{26, 7,19, 8,26,20, 7, 8},{ 4,20,26, 8, 4,26,17, 8}}, 
																 {{10, 2,21, 8,21, 2,11, 8},{21, 1,10, 8,21, 9, 1, 8},{21,12, 0, 8,21, 0, 9, 8},{12,21, 3, 8,21,11, 3, 8}}
															 };	
	for(  int i_index = 0  ;  i_index < 8  ;  i_index++   ){
		tet_index[ i_index ] = indexes[ src_dir ][ position ][ i_index ];
	}
}

/**
 *Getting local indexes of tetrahedron formed on interfaz when neighbours by face are same 
 *refined, and also neighbours by edge are same or less refined, in this case the edges 
 *direction depending on i_edge for each face source direction are: 
 *face 0: 11- 7- 9-13
 *face 1: 10- 6- 8-12
 *face 2: 15- 9- 8-17
 *face 3: 14- 6- 7-16
 *face 4: 16-13-17-12
 *face 5: 14-11-15-10
 *In this case is created one tetrahedron
 *@param[in] src_dir Is the source direction from other rank
 *@param[in] i_edge Is the index to the edge direction where the tetrahedron will be created
 *@param[out] tet_index Is the array where the tetrahedron indexes will be set
 */
void GetLocalTetIndexOnInterfaceWhenNeighbourIsSameRefinedAndEdgeSame( int src_dir , int i_edge , int* tet_index ){
	int indexes[ 6 ][ 4 ][ 4 ] = { { {  1,23, 2, 8 } , {  5,23, 1, 8 } , {  2,23, 6, 8 } , {  6,23, 5, 8 } } , 
														 		 { {  3,25, 0, 8 } , {  0,25, 4, 8 } , {  7,25, 3, 8 } , {  4,25, 7, 8 } } ,
														 		 { {  2,24, 3, 8 } , {  6,24, 2, 8 } , {  3,24, 7, 8 } , {  7,24, 6, 8 } } ,
														 		 { {  0,22, 1, 8 } , {  4,22, 0, 8 } , {  1,22, 5, 8 } , {  5,22, 4, 8 } } ,
														 		 { {  4,26, 5, 8 } , {  5,26, 6, 8 } , {  6,26, 7, 8 } , {  7,26, 4, 8 } } ,
														 		 { {  1,21, 0, 8 } , {  2,21, 1, 8 } , {  3,21, 2, 8 } , {  0,21, 3, 8 } } };
	for(  int i_index = 0  ;  i_index < 4  ;  i_index++   ){
		tet_index[ i_index ] = indexes[ src_dir ][ i_edge ][ i_index ];
	}
}

/**
 *Getting local indexes of tetrahedra formed on interfaz when neighbours by face are same 
 *refined, and neighbours by edge are more or refined, in this case the edges direction
 *depending on i_edge for each face source direction are: 
 *face 0: 11- 7- 9-13
 *face 1: 10- 6- 8-12
 *face 2: 15- 9- 8-17
 *face 3: 14- 6- 7-16
 *face 4: 16-13-17-12
 *face 5: 14-11-15-10
 *In this case is created one tetrahedron
 *@param[in] src_dir Is the source direction from other rank
 *@param[in] i_edge Is the index to the edge direction where the tetrahedron will be created
 *@param[out] tet_index Is the array where the tetrahedron indexes will be set
 */
void GetLocalTetIndexOnInterfaceWhenNeighbourIsSameRefinedAndEdgeMore( int src_dir , int i_edge , int* tet_index ){
	int indexes[ 6 ][ 4 ][ 8 ] = { { { 1,23,10,8,10,23, 2,8 } , { 5,23,14,8,14,23, 1,8 } , {  6,15,23,8,23,15, 2,8 } , { 6,23,18,8,18,23, 5,8 } } , 
														 		 { { 0,12,25,8,12, 3,25,8 } , { 0,25,13,8,13,25, 4,8 } , {  3,16,25,8,16, 7,25,8 } , { 7,20,25,8,20, 4,25,8 } } ,
														 		 { { 2,24,11,8,11,24, 3,8 } , { 6,24,15,8,15,24, 2,8 } , {  3,24,16,8,16,24, 7,8 } , { 7,24,19,8,19,24, 6,8 } } ,
														 		 { { 9,22, 1,8, 9, 0,22,8 } , { 4,22,13,8,13,22, 0,8 } , {  5,14,22,8, 1,22,14,8 } , { 4,17,22,8,17, 5,22,8 } } ,
														 		 { { 4,26,17,8,17,26, 5,8 } , { 5,26,18,8,18,26, 6,8 } , {  6,26,19,8,19,26, 7,8 } , { 7,26,20,8,20,26, 4,8 } } ,
														 		 { { 1,21, 9,8, 9,21, 0,8 } , { 2,21,10,8,10,21, 1,8 } , { 11,21, 2,8, 3,21,11,8 } , { 0,21,12,8,12,21, 3,8 } } };
	for(  int i_index = 0  ;  i_index < 8  ;  i_index++   ){
		tet_index[ i_index ] = indexes[ src_dir ][ i_edge ][ i_index ];
	}
}

/**
 *Getting indexes of tetrahedra created on intrface cells more refined, in this case the 
 *cell was searched using direction src_dir, and this is the cell i_cell
 *@param[in] i_cell Is the cell index (from 0-3) on interface face
 *@param[in] src_dir Is teh search direction from where cell was obtained
 *@param[out] tet_index Is the array where tetrahedral indexes will be set  
 */
void GetLocalTetIndexOnInterfaceWhenNeighbourIsLessRefined( int i_cell , int src_dir , int* tet_index ){
	int indexes[ 6 ][ 4 ][ 8 ] = { { { 1,5,6,8,1,6,2,8 } , { 1,5,2,8,5,6,2,8 } , { 1,5,6,8,1,6,2,8 } , { 1,5,2,8,5,6,2,8 } } , 
														 		 { { 7,4,0,8,7,0,3,8 } , { 7,4,3,8,4,0,3,8 } , { 7,4,0,8,7,0,3,8 } , { 7,4,3,8,4,0,3,8 } } ,
														 		 { { 7,3,6,8,6,3,2,8 } , { 2,6,7,8,2,7,3,8 } , { 7,3,6,8,6,3,2,8 } , { 2,6,7,8,2,7,3,8 } } ,
														 		 { { 4,5,0,8,0,5,1,8 } , { 0,4,1,8,1,4,5,8 } , { 4,5,0,8,0,5,1,8 } , { 0,4,1,8,1,4,5,8 } } ,
														 		 { { 4,7,6,8,4,6,5,8 } , { 4,7,5,8,7,6,5,8 } , { 4,7,6,8,4,6,5,8 } , { 4,7,5,8,7,6,5,8 } } ,
														 		 { { 3,0,2,8,2,0,1,8 } , { 3,0,1,8,3,1,2,8 } , { 3,0,2,8,2,0,1,8 } , { 3,0,1,8,3,1,2,8 } } };
	for(  int i_index = 0  ;  i_index < 8  ;  i_index++   ){
		tet_index[ i_index ] = indexes[ src_dir ][ i_cell ][ i_index ];
	}
}

/**
 *Getting qudratin indexes on face
 *@param[in] direction Is the search direction where the quadratq nodes are requested
 *@param[out] index Is the array where the nodes willl be setted
 */
void GetQuadraticIndexOnFace( int direction , int* index ){
	switch( direction ){
		case 0:
			index[ 0 ] = 12;
			index[ 1 ] = 13;
			index[ 2 ] = 16;
			index[ 3 ] = 20;
			break;
		case 1:
			index[ 0 ] = 10;
			index[ 1 ] = 14;
			index[ 2 ] = 15;
			index[ 3 ] = 18;
			break;
		case 2:
			index[ 0 ] =  9;
			index[ 1 ] = 13;
			index[ 2 ] = 14;
			index[ 3 ] = 17;
			break;
		case 3:
			index[ 0 ] = 11;
			index[ 1 ] = 15;
			index[ 2 ] = 16;
			index[ 3 ] = 19;
			break;
		case 4:
			index[ 0 ] = 9;
			index[ 1 ] = 10;
			index[ 2 ] = 11;
			index[ 3 ] = 12;
			break;
		case 5:
			index[ 0 ] = 17;
			index[ 1 ] = 18;
			index[ 2 ] = 19;
			index[ 3 ] = 20;
			break;
		default:
			assert( 0 );
			break;		
	}
	
}

/**
 *Getting tetrahedra conectivities on boundary cell when edge is more refined
 *@param[in] direction Is the face where the tetrahedra are created
 *@param[in] i_edge Is the edge where tetrahedra are created
 *@param[out] tet_index Is the array where tetrahedra index will be setted
 */
void GetLocalTetIndexOnBoundaryCellWhenEdgeIsMoreRefined( int direction , int i_edge , int* tet_index ){
	int indexes[ 6 ][ 4 ][ 8 ] = { { {  0,12,25,8,12, 3,25,8 } , {  4,13,25,8,13, 0,25,8 } , {  7,25,16,8,16,25, 3,8 } , { 20, 4,25,8, 7,20,25,8 } } , 
														 		 { {  1,23,10,8,10,23, 2,8 } , {  1,14,23,8,14, 5,23,8 } , {  6,15,23,8,15, 2,23,8 } , {  6,23,18,8,18,23, 5,8 } } ,
														 		 { {  0,22, 9,8, 9,22, 1,8 } , {  4,22,13,8,13,22, 0,8 } , { 14,22, 5,8, 1,22,14,8 } , {  4,17,22,8,17, 5,22,8 } } ,
														 		 { {  2,24,11,8,11,24, 3,8 } , {  2,15,24,8,15, 6,24,8 } , {  3,24,16,8,16,24, 7,8 } , {  7,24,19,8,19,24, 6,8 } } ,
														 		 { {  0, 9,21,8, 9, 1,21,8 } , {  1,10,21,8,10, 2,21,8 } , {  2,11,21,8,11, 3,21,8 } , {  3,12,21,8,12, 0,21,8 } } ,
														 		 { {  4,26,17,8,17,26, 5,8 } , {  5,26,18,8,18,26, 6,8 } , {  6,26,19,8,19,26, 7,8 } , {  7,26,20,8,20,26, 4,8 } } };

	for(  int i_index = 0  ;  i_index < 8  ;  i_index++   ){
		tet_index[ i_index ] = indexes[ direction ][ i_edge ][ i_index ];
	}
}


/**
 *Getting tetrahedron conectivities on boundary cell when edge is same refined
 *@param[in] direction Is the face where the tetrahedra are created
 *@param[in] i_edge Is the edge where tetrahedron are created
 *@param[out] tet_index Is the array where tetrahedron index will be setted
 */
void GetLocalTetIndexOnBoundaryCellWhenEdgeIsSameRefined( int direction , int i_edge , int* tet_index ){
	int indexes[ 6 ][ 4 ][ 4 ] = { { { 0,3,25,8 } , { 4,0,25,8 } , { 3,7,25,8 } , { 7,4,25,8 } } , 
														 		 { { 2,1,23,8 } , { 1,5,23,8 } , { 6,2,23,8 } , { 5,6,23,8 } } ,
														 		 { { 1,0,22,8 } , { 0,4,22,8 } , { 5,1,22,8 } , { 4,5,22,8 } } ,
														 		 { { 3,2,24,8 } , { 2,6,24,8 } , { 7,3,24,8 } , { 6,7,24,8 } } ,
														 		 { { 0,1,21,8 } , { 1,2,21,8 } , { 2,3,21,8 } , { 3,0,21,8 } } ,
														 		 { { 5,4,26,8 } , { 6,5,26,8 } , { 7,6,26,8 } , { 4,7,26,8 } } };

	for(  int i_index = 0  ;  i_index < 4  ;  i_index++   ){
		tet_index[ i_index ] = indexes[ direction ][ i_edge ][ i_index ];
	}
}

/**
 *In this case are created two tetraedra, need to be correct oriented in order to have all 
 *tetrahedra edges conformal between neighbour elements
 *@param[in] position Is the relative position of src_cell on neighbour cell
 *@param[in] direction Is the direction where neighbour is located
 *@param[out] tet_index Is the array to set the tetrahedra indexes
 */
void GetLocalTetIndexWhenNeighbourIsLessRefined( int position , int direction , int* tet_index ){
	int indexes[6][4][8] = {{{7,4,3,8,4,0,3,8},{7,4,0,8,7,0,3,8},{7,4,3,8,4,0,3,8},{7,4,0,8,7,0,3,8}} , 
												  {{1,5,6,8,1,6,2,8},{1,5,2,8,5,6,2,8},{1,5,6,8,1,6,2,8},{1,5,2,8,2,5,6,8}} ,
													{{0,4,5,8,0,5,1,8},{0,4,1,8,1,4,5,8},{0,4,5,8,0,5,1,8},{0,4,1,8,1,4,5,8}} ,
													{{2,6,7,8,2,7,3,8},{7,3,6,8,6,3,2,8},{2,6,7,8,2,7,3,8},{7,3,6,8,6,3,2,8}} ,
													{{0,1,3,8,1,2,3,8},{0,1,2,8,0,2,3,8},{0,1,3,8,1,2,3,8},{0,1,2,8,0,2,3,8}} ,
													{{7,6,4,8,6,5,4,8},{7,6,5,8,7,5,4,8},{7,6,4,8,6,5,4,8},{7,6,5,8,7,5,4,8}} };

	for(  int i_index = 0  ;  i_index < 8  ;  i_index++   ){
		tet_index[ i_index ] = indexes[ direction ][ position ][ i_index ];
	}
}

/**
 *Checkink if source cell is lower than neighbour cell
 *@param[in] src_cell Is the cell that is tested to know is is the lower or not
 *@param[in] other_cell Is the cell that is used to test src_cell
 *@return A bool value indicating is src_cell is lower or not
 */
bool SourceCellIsLower( OctreeCell* src_cell , OctreeCell* other_cell ){
	key_type src_key[ 3 ];	
	key_type other_key[ 3 ];
	src_cell->GetMinKey(     src_key[ 0 ] ,   src_key[ 1 ] ,   src_key[ 2 ] );
	other_cell->GetMinKey( other_key[ 0 ] , other_key[ 1 ] , other_key[ 2 ] );
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		if(  src_key[ i_dim ] < other_key[ i_dim ]  ){
			return true;
		}
	}
	return false;
}

/**
 *Getting tetrahedra indexes when neighbour by face is same refined, and by edge is same 
 *or less. In this case only one tetrahedra is created
 *@param[in] direction Is the direction where the neighbour is located 
 *@param[in] i_edge Is the edge used to create tetrahedra on face located over direction 
 *@param[out] tet_index Is where the tetrahedron indexes will be setted
 */
void GetLocalTetIndexOnCellWhenFaceIsSameRefinedAndEdgeSameOrLess( int direction , int i_edge , int* tet_index ){
	int indexes[6][4][4] = {{{3,0,8,8},{0,4,8,8},{7,3,8,8},{4,7,8,8}} , 
												  {{1,2,8,8},{5,1,8,8},{2,6,8,8},{6,5,8,8}} ,
													{{0,1,8,8},{4,0,8,8},{1,5,8,8},{5,4,8,8}} ,
													{{2,3,8,8},{6,2,8,8},{3,7,8,8},{7,6,8,8}} ,
													{{1,0,8,8},{2,1,8,8},{3,2,8,8},{0,3,8,8}} ,
													{{4,5,8,8},{5,6,8,8},{6,7,8,8},{7,4,8,8}} };

	for(  int i_index = 0  ;  i_index < 4  ;  i_index++   ){
		tet_index[ i_index ] = indexes[ direction ][ i_edge ][ i_index ];
	}
}

/**
 *Getting tetrahedra indexes when neighbour by face is same refined, and by edge is more 
 *refined. In this case only two tetrahedra are created
 *@param[in] direction Is the direction where the neighbour is located 
 *@param[in] i_edge Is the edge used to create tetrahedra on face located over direction 
 *@param[out] tet_index Is where the tetrahedron indexes will be setted
 */
void GetLocalTetIndexOnCellWhenFaceIsSameRefinedAndEdgeMore( int direction , int i_edge , int* tet_index ){
	int indexes[6][4][8] = {{{3,12,8,8,12,0,8,8},{0,13,8,8,13,4,8,8},{7,16,8,8,16,3,8,8},{4,20,8,8,20,7,8,8}} , 
												  {{1,10,8,8,10,2,8,8},{5,14,8,8,14,1,8,8},{2,15,8,8,15,6,8,8},{6,18,8,8,18,5,8,8}} ,
													{{0, 9,8,8, 9,1,8,8},{4,13,8,8,13,0,8,8},{1,14,8,8,14,5,8,8},{5,17,8,8,17,4,8,8}} ,
													{{2,11,8,8,11,3,8,8},{6,15,8,8,15,2,8,8},{3,16,8,8,16,7,8,8},{7,19,8,8,19,6,8,8}} ,
													{{1, 9,8,8, 9,0,8,8},{2,10,8,8,10,1,8,8},{3,11,8,8,11,2,8,8},{0,12,8,8,12,3,8,8}} ,
													{{4,17,8,8,17,5,8,8},{5,18,8,8,18,6,8,8},{6,19,8,8,19,7,8,8},{7,20,8,8,20,4,8,8}} };
	for(  int i_index = 0  ;  i_index < 8  ;  i_index++   ){
		tet_index[ i_index ] = indexes[ direction ][ i_edge ][ i_index ];
	}
}

/**
 *Getting tetrahedra indexes when neighbour by face is more refined, in this case eight 
 *tetrahedra are created on that face.
 *@param[in] direction Indicates the face where the tetrahedra are created
 *@param[out] tet_index Is the array where tetrahedra conectivities will be setted
 */
void GetLocalTetIndexOnFaceWhenNeighbourIsMoreRefined( int direction , int* tet_index ){
	int indexes[6][32] = { { 0,12,25, 8,12, 3,25, 8, 4,13,25, 8,13, 0,25, 8, 3,16,25, 8,16, 7,25, 8, 7,20,25, 8,20, 4,25, 8} ,
												 { 2,10,23, 8,10, 1,23, 8, 1,14,23, 8,14, 5,23, 8, 6,15,23, 8,15, 2,23, 8, 5,18,23, 8,18, 6,23, 8} ,
												 { 1, 9,22, 8, 9, 0,22, 8, 0,13,22, 8,13, 4,22, 8, 5,14,22, 8,14, 1,22, 8, 4,17,22, 8,17, 5,22, 8} ,
												 { 3,11,24, 8,11, 2,24, 8, 2,15,24, 8,15, 6,24, 8, 7,16,24, 8,16, 3,24, 8, 6,19,24, 8,19, 7,24, 8} ,
												 { 0, 9,21, 8, 9, 1,21, 8, 1,10,21, 8,10, 2,21, 8, 2,11,21, 8,11, 3,21, 8, 3,12,21, 8,12, 0,21, 8} ,
												 { 5,17,26, 8,17, 4,26, 8, 6,18,26, 8,18, 5,26, 8, 7,19,26, 8,19, 6,26, 8, 4,20,26, 8,20, 7,26, 8} };
	for(  int i_index = 0  ;  i_index < 32  ;  i_index++   ){
		tet_index[ i_index ] = indexes[ direction ][ i_index ];
	}
}




