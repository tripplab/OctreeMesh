#include "mesher.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                   COMUNICATOR METHODS                                //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Mesher::Mesher(){
	strcpy( input_ , "\0" );
	strcpy( output_ , "\0" );
	r_level_ = 0 ;
	fold_ = -1;
	fold_index_ = -1;
	octree_ = NULL;
	local_root_ = NULL;
	vdb_ = NULL;
	msh_ = NULL;
}

/**
 *Constructor receiving program arguments 
 *This constructor initializes the parallel region in each process, also reads the vdb 
 *file containing the capside information
 *@param[in] argc Is the number inidicating the program arguments
 *@param[in] argv Is the list of arguments for the program
 */
Mesher::Mesher( int argc , char** argv ){
	
	msh_ = NULL;
	local_root_ = NULL;
	input_ = argv[ 1 ];
	vdb_ = new Vdb( input_ , 1 , 1 , atoi( argv[ 2 ] ) , atoi( argv[ 3 ] ) );
	virus_ = argv[ 7 ];
	fold_ = atoi( argv[ 5 ] );
	fold_index_ = atoi( argv[ 6 ] );
	resolution_capside_ = atof( argv[ 4 ] );
	resolution_ = this->InterpolateCapsidResolution(  );

	cone_amplitude_ = atof( argv[ 8 ] );
	cone_amplitude_ = cone_amplitude_ * 3.14159265358979323846264338327950 /180.00;
	prop_variation_ = atof( argv[ 9 ] );
	young_modulus_ = atof( argv[ 10 ] );

	//Checking if folds are correct values
	assert(  ( fold_ == 5 )  ||  ( fold_ == 3 )  ||  ( fold_ == 2 )  );
	assert(  ( fold_index_ >= 0 )  &&  ( fold_index_ < 12 )  );

	//Checking for correct values on virus name and cone amplitude
	assert( !strcmp( virus_ , "1CWP" ) || !strcmp( virus_ , "4G93" ) || 
					!strcmp( virus_ , "3IZG" ) || !strcmp( virus_ , "3J4U" ) );
	assert( ( cone_amplitude_ >= 0.0 ) && 
					( cone_amplitude_ <= 3.14159265358979323846264338327950*0.5 ) );
	//Creating octree
	octree_ = new OctreeDriver();

	local_root_ = octree_->GetRoot();

}

/**
 *Default destructor 
 *Finalizes the parallel region
 */
Mesher::~Mesher(  ){
}

//GETS

/**
 *Getting all leaves on local root
 *@param[out] leaves Is the vector that will contain all leaves on local root
 */
void Mesher::GetAllLocalRootLeaves( OctreeCell_vector& leaves ){
	OctreeCell_vector cells_stack;
	cells_stack.push_back( local_root_ );
	while( !cells_stack.empty() ){
		OctreeCell* cell = cells_stack.back();
		cells_stack.pop_back();
		if( cell->HasChildren() ){
			for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
				cells_stack.push_back( cell->pGetChild( i_child ) );
			}
		}else{
			leaves.push_back( cell );
		}
	}
}

/**
 *Getting cell on local octree using key
 *This method takes a key and returns a pointer to the cell that contains the key
 *@param[in] keys Are the keys to look for the cell
 *@return A OctreeCell pointer to the cell containning keys
 */
OctreeCell* Mesher::GetCellOnLocalOctree( key_type* keys ){

	OctreeCell* cell = local_root_;
	for(  size_t i = 0  ;  i < ROOT_LEVEL  ;  i++  ){
		if(  cell->IsLeaf()  ){
			return cell;
		}
		cell = cell->pGetChild( keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
	}
	return cell;
}

/**
 *Getting cone direction normalized
 *Direction is given by the vector that indicates the fold to be aligned with the Y axis
 *@param[out] dir Is the vector to store the normalized direction vector 
 */
void Mesher::GetConeDirection( double* dir ){
	switch(  fold_  ){
		case 2:
			this->GetConeDirectionFold2( dir );
			break;
		case 3:
			this->GetConeDirectionFold3( dir );
			break;
		case 5:
			this->GetConeDirectionFold5( dir );
			break;
	}
	Normalize( dir );
}

/**
 *Getting cone direction normalized on fold 2
 *@param[out] dir Is the normalized vector indicating the cone direction
 */
void Mesher::GetConeDirectionFold2( double* dir ){
	switch( fold_index_ ){
		case 0:
			dir[ 0 ] = 0.000; dir[ 1 ] = 0.000; dir[ 2 ] = 123.301;
			break;
		case 1:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 2:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 3:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 4:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 5:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 6:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 7:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 8:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 9:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 10:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 11:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
	}
}

/**
 *Getting cone direction normalized on fold 3
 *@param[out] dir Is the normalized vector indicating the cone direction
 */
void Mesher::GetConeDirectionFold3( double* dir ){
	switch( fold_index_ ){
		case 0:
			dir[ 0 ] = 41.100; dir[ 1 ] = 0.000; dir[ 2 ] = 107.601;
			break;
		case 1:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 2:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 3:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 4:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 5:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 6:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 7:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 8:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 9:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 10:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
		case 11:
			std::cout << "Fold " << fold_ << " index " << fold_index_ << " not known " << std::endl;
			assert( 0 );
			break;
	}
}

/**
 *Getting cone direction normalized on fold 5
 *@param[out] dir Is the normalized vector indicating the cone direction
 */
void Mesher::GetConeDirectionFold5( double* dir ){
	switch( fold_index_ ){
		case 0:
			dir[ 0 ] = 000.000; dir[ 1 ] =  77.316; dir[ 2 ] = 125.100;
			break;
		case 1:
			dir[ 0 ] = 000.000; dir[ 1 ] = -77.316; dir[ 2 ] = 125.100;
			break;
		case 2:
			dir[ 0 ] = -77.316; dir[ 1 ] =  125.101; dir[ 2 ] = 000.000;
			break;
		case 3:
			dir[ 0 ] =  77.316; dir[ 1 ] =  125.100; dir[ 2 ] = 000.000;
			break;
		case 4:
			dir[ 0 ] =-125.101; dir[ 1 ] =  000.000; dir[ 2 ] =  77.316;
			break;
		case 5:
			dir[ 0 ] = 125.100; dir[ 1 ] =  000.000; dir[ 2 ] =  77.316;
			break;
		case 6:
			dir[ 0 ] = -77.316; dir[ 1 ] = -125.101; dir[ 2 ] = 000.000;
			break;
		case 7:
			dir[ 0 ] =  77.316; dir[ 1 ] = -125.101; dir[ 2 ] = 000.000;
			break;
		case 8:
			dir[ 0 ] =-125.101; dir[ 1 ] =  000.000; dir[ 2 ] = -77.316;
			break;
		case 9:
			dir[ 0 ] = 000.000; dir[ 1 ] =   77.316; dir[ 2 ] =-125.100;
			break;
		case 10:
			dir[ 0 ] = 000.000; dir[ 1 ] =  -77.316; dir[ 2 ] =-125.100;
			break;
		case 11:
			dir[ 0 ] = 125.100; dir[ 1 ] =  000.000; dir[ 2 ] = -77.316;
			break;
	}
}

/**
 *Getting number of nodes in the final mesh
 *@return A size_t with the number of nodes in the final mesh
 */
size_t Mesher::GetNNodesInFinalMesh(  ){
	return n_final_nodes_;
}

/**
 *Getting number of nodes in the final mesh
 *@return A size_t with the number of nodes in the final mesh
 */
size_t Mesher::GetNElementsInFinalMesh(  ){
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );	
	size_t n_leaves = leaves.size();
	size_t index = 0;
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		OctreeCell* cell = leaves[ i_leaf ];
		DataInsideOctreeCell* data = cell->pGetData();
		if(  data->_GetIsMesh() == 1  ){
			index++;
		}
	}
	return index;
}

//SETS
void Mesher::SetAllOctreeNodes( ){
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );	
	size_t n_leaves = leaves.size();
	int index = 0;
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		OctreeCell* cell = leaves[ i_leaf ];
		DataInsideOctreeCell* data = cell->pGetData();
		if(  data->_GetIsMesh() == 1  ){
			this->SetNodesOnCellAndNeighbours( cell , &index );
		}
	}
	n_final_nodes_ = (size_t)index;
}

/**
 *Setting nodes on cell and its neighbours
 *@param[in] cell Is the cell where the nodes will be set
 *@param[in] index Is the actual index of the nodes in the list of nodes
 */
void Mesher::SetNodesOnCellAndNeighbours( OctreeCell* cell , int* index ){
	double scale = 1.00 / ( 1 << ROOT_LEVEL );
	for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
		key_type keys[ 3 ];
		cell->GetKey( (size_t)i_node , keys );
		double* coord = new double [ 3 ];
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			coord[ i_dim ] = ( static_cast<double>( keys[ i_dim ] ) * scale );
		}
		if(  !this->SetNodeInCellAndNeighbors( i_node , cell , coord , (*index) )  ){
			delete[] coord;
		}else{
			(*index)++;
		}
	}
}

/**
 *Setting a node on a cell and all its neighbors
 *This method set a node on a cell and all its neighbors, in this case exist at most 8 
 *neighbours
 *@param[in] node Is the node index to be set on the local cell
 *@param[in] cell Is the cell where the node have to be set
 *@param[in] coord Is the coordinate to be set 
 *@param[in] interfaz Is a bool value indicating if a node is interfaz beetwen process
 */
bool Mesher::SetNodeInCellAndNeighbors( int node , OctreeCell* cell , double* coord , int index ){

	DataInsideOctreeCell* data = cell->pGetData();
	if(  data->_IsSetNode( node )  ){
		return false; 
	}

	this->SetNodeInCell( node , cell , coord , index );

	key_type point_key[ 3 ];
	cell->GetKey( (size_t)node , point_key );
	for(  size_t i_direction = 0  ;  i_direction < 8  ;  i_direction++  ){
		key_type neighbour_key[ 3 ];
		if(  cell->GetNeighbourKey( (size_t)node , i_direction , neighbour_key )  ){
			OctreeCell* neighbour_cell = this->GetCellOnLocalOctree( neighbour_key );
			if(  !neighbour_cell || ( neighbour_cell == cell )  ){
				continue;
			}
			DataInsideOctreeCell* data_neighbour = neighbour_cell->pGetData();
			if(  data_neighbour->_GetIsMesh() != 1  ){
				continue;
			}
			size_t position = neighbour_cell->GetLocalPosition( point_key );
			this->SetNodeInCell( (int)position , neighbour_cell , coord , index );
		}
	}
	return true;
}

/**
 *Setting a node on a cell
 *This method set a node on a cell
 *@param[in] node Is the local position where the node will be set
 *@param[in] cell Is the cell where the node will be set
 *@param[in] coord Is the pointer to the coord that will be set 
 *@param[in] interfaz Is a bool value indicating if a node is interfaz beetwen process
 */
void Mesher::SetNodeInCell( int node , OctreeCell* cell , double* coord , int index ){

	DataInsideOctreeCell* data = cell->pGetData();
	data->_SetNode( node , coord );
	data->_SetNodeIndex( node , index );
}

/**
 *Setting that node was printed in actual cell and its neighbours
 *@param[in] node Is the local index of node to be set
 *@param[in] cell Is the cell pointer where the information will be set
 */
void Mesher::SetPrintedOnNodeAndNeighbours( int node , OctreeCell* cell ){

	this->SetNodePrintedInCell( node , cell );

	key_type point_key[ 3 ];
	cell->GetKey( (size_t)node , point_key );
	for(  size_t i_direction = 0  ;  i_direction < 8  ;  i_direction++  ){
		key_type neighbour_key[ 3 ];
		if(  cell->GetNeighbourKey( (size_t)node , i_direction , neighbour_key )  ){
			OctreeCell* neighbour_cell = this->GetCellOnLocalOctree( neighbour_key );
			if(  !neighbour_cell || ( neighbour_cell == cell )  ){
				continue;
			}
			DataInsideOctreeCell* data_neighbour = neighbour_cell->pGetData();
			if(  data_neighbour->_GetIsMesh() != 1  ){
				continue;
			}
			size_t position = neighbour_cell->GetLocalPosition( point_key );
			if(  !data_neighbour->_NodeWasPrinted( (int)position )  ){
				this->SetNodePrintedInCell( (int)position , neighbour_cell );
			}
		}
	}
}

/**
 *Setting that fixed node was printed in actual cell and its neighbours
 *@param[in] node Is the local index of node to be set
 *@param[in] cell Is the cell pointer where the information will be set
 */
void Mesher::SetPrintedOnFixedNodeAndNeighbours( int node , OctreeCell* cell ){

	this->SetFixedNodePrintedInCell( node , cell );
	key_type point_key[ 3 ];
	cell->GetKey( (size_t)node , point_key );
	for(  size_t i_direction = 0  ;  i_direction < 8  ;  i_direction++  ){
		key_type neighbour_key[ 3 ];
		if(  cell->GetNeighbourKey( (size_t)node , i_direction , neighbour_key )  ){
			OctreeCell* neighbour_cell = this->GetCellOnLocalOctree( neighbour_key );
			if(  !neighbour_cell || ( neighbour_cell == cell )  ){
				continue;
			}
			DataInsideOctreeCell* data_neighbour = neighbour_cell->pGetData();
			if(  data_neighbour->_GetIsMesh() != 1  ){
				continue;
			}
			size_t position = neighbour_cell->GetLocalPosition( point_key );
			if(  !data_neighbour->_NodeFixedWasPrinted( (int)position )  ){
				this->SetFixedNodePrintedInCell( (int)position , neighbour_cell );
			}
		}
	}
}

/**
 *Setting node printed on cell
 *@param[in] node Is the local node index
 *@param[in] cell is the cell where the node printed will be set
 */
void Mesher::SetNodePrintedInCell( int node , OctreeCell* cell ){
	DataInsideOctreeCell* data = cell->pGetData();
	data->_SetPrinted( node );
}

/**
 *Setting fixed node printed on cell
 *@param[in] node Is the local node index
 *@param[in] cell is the cell where the node printed will be set
 */
void Mesher::SetFixedNodePrintedInCell( int node , OctreeCell* cell ){
	DataInsideOctreeCell* data = cell->pGetData();
	data->_SetPrintedFixed( node );
}

/**
 *This function set all fixed and loaded elements in the mesh
 *The function uses the cone information and calculates whether intersects or not
 */
void Mesher::SetLoadedAndFixedElements(  ){
	double tot_vol = 0.0;//Total volume of capsid mesh
	double prop_vol = 0.0;//proportion of volume to be loaded
	double vertex[ 3 ] = { 0.5 , 0.5 , 0.5 };//cone vertex
	double direction[ 3 ] = { 0.0 , 0.0 , 0.0 };//cone direction
	double cell_vol;
	//Obtaining total volume
	tot_vol =  this->CalculateOctreeMeshVolume( &cell_vol );
	prop_vol = this->CalculateProportionLoaded(  );
	this->GetConeDirection( direction );
	//Seting loaded elements
	this->SetLoadedElements( tot_vol , prop_vol , vertex , direction , cell_vol/tot_vol );
	//Seting fixed elements, the cone direction is inverted because the fixed elements are
  //in the oposite side than the loaded elements
  direction[ 0 ] *= -1.0;		direction[ 1 ] *= -1.0;		direction[ 2 ] *= -1.0;
	this->SetFixedElements( tot_vol , prop_vol , vertex , direction , cell_vol/tot_vol );
}

/**
 *Setting loaded elements based on the intersection between the cone and the elements
 *Loaded elements are set as all elements intersecting the cone until reach a prestablished
 *proportion, the elements in the interior of the capside are removed until reach the 
 *proportion required
 *@param[in] tot_vol Is the total volume of the final mesh
 *@param[in] prop_vol Is the proportion of volume to be loaded
 *@param[in] vertex Is the origin of the cone
 *@param[in] direction Is the direction of the cone
 *@param[in] cell_prop Is the proportion of volume added by each hexahedron
 */
bool Mesher::SetLoadedElements( double tot_vol , double prop_vol , double* vertex , 
											  						 double* direction , double cell_prop ){
	/*First all cells intersecting the cone are marked as loaded and is calculated the 
	  distance between its center and the cone vertex in order to determine whish need to be
	  removed from the selection                                                          */
	std::vector<double> distance;
	std::vector<size_t> leaf_id;
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	size_t n_leaves = leaves.size();
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		OctreeCell* leaf = leaves[ i_leaf ];
		DataInsideOctreeCell* data = leaf->pGetData();
		if(  data->_GetIsMesh(  ) == 1 ){
			//Obtaining cell center
			double min_point[ 3 ];
			double max_point[ 3 ];
			double center[ 3 ];
			leaf->GetMinPointNormalized( min_point ); 
			leaf->GetMaxPointNormalized( max_point );
			center[ 0 ] = 0.5 * ( min_point[ 0 ] + max_point[ 0 ] ); 
			center[ 1 ] = 0.5 * ( min_point[ 1 ] + max_point[ 1 ] ); 
			center[ 2 ] = 0.5 * ( min_point[ 2 ] + max_point[ 2 ] ); 
			//Testing whether cell center intersects cone or not
			if(  ConeIntersectsNode( vertex , direction , cone_amplitude_ , center )  ){
				//Add as loaded cell
				double dist = DistancePointPoint( vertex , center );
				data->_SetLoadedElement( dist );
				distance.push_back( dist );
				leaf_id.push_back( i_leaf );	
			}
		}
	}

	//Copying distance and index vector to a simply vector and sorting based on the distance
	size_t n_positions = distance.size();
	size_t* leaf_index =  new size_t[ n_positions ];
	double* v_dist =  new double[ n_positions ];
	for(  size_t i_pos = 0  ;  i_pos < n_positions  ; i_pos++    ){
		leaf_index[ i_pos ] = leaf_id[ i_pos ];
		v_dist[ i_pos ] = distance[ i_pos ];
	}
	leaf_id.clear();
	distance.clear();
	ShellSortSizeTDouble( leaf_index , v_dist , n_positions );

	//Actualizing loaded proportion aplying the prop_variation_
	prop_vol += prop_variation_;

	//Erasing cells in order to accomplish the total proportion
	double actual_prop = cell_prop * n_positions;
	if(  actual_prop < prop_vol  ){
		assert( 0 );
		return true;
	}	

	for(  size_t i_pos = 0  ;  i_pos < n_positions  ;  i_pos++  ){
		if(  actual_prop > prop_vol  ){
			DataInsideOctreeCell* data = leaves[ leaf_index[ i_pos ] ]->pGetData();
			data->_SetNotLoadedElement( );
			actual_prop -= cell_prop;
		}else{
			break;
		}
	}

	delete[] leaf_index;
	delete[] v_dist;
	return true;
}

/**
 *Setting fixed elements based on the intersection between the cone and the elements
 *Fixed elements are set as all elements intersecting the cone until reach a prestablished
 *proportion, the elements in the interior of the capside are removed until reach the 
 *proportion required
 *@param[in] tot_vol Is the total volume of the final mesh
 *@param[in] prop_vol Is the proportion of volume to be loaded
 *@param[in] vertex Is the origin of the cone
 *@param[in] direction Is the direction of the cone
 *@param[in] cell_prop Is the proportion of volume added by each hexahedron
 */
bool Mesher::SetFixedElements( double tot_vol , double prop_vol , double* vertex , 
													 					double* direction , double cell_prop ){
	/*First all cells intersecting the cone are marked as fixed and is calculated the 
	  distance between its center and the cone vertex in order to determine which need to be
	  removed from the selection                                                          */
	std::vector<double> distance;
	std::vector<size_t> leaf_id;
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	size_t n_leaves = leaves.size();
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		OctreeCell* leaf = leaves[ i_leaf ];
		DataInsideOctreeCell* data = leaf->pGetData();
		if(  data->_GetIsMesh(  ) == 1 ){
			//Obtaining cell center
			double min_point[ 3 ];
			double max_point[ 3 ];
			double center[ 3 ];
			leaf->GetMinPointNormalized( min_point ); 
			leaf->GetMaxPointNormalized( max_point );
			center[ 0 ] = 0.5 * ( min_point[ 0 ] + max_point[ 0 ] ); 
			center[ 1 ] = 0.5 * ( min_point[ 1 ] + max_point[ 1 ] ); 
			center[ 2 ] = 0.5 * ( min_point[ 2 ] + max_point[ 2 ] ); 
			//Testing whether cell center intersects cone or not
			if(  ConeIntersectsNode( vertex , direction , cone_amplitude_ , center )  ){
				//Add as fixed cell
				double dist = DistancePointPoint( vertex , center );
				data->_SetFixedElement( dist );
				distance.push_back( dist );
				leaf_id.push_back( i_leaf );
			}
		}
	}
	//Copying distance and index vector to a simply vector and sorting based on the distance
	size_t n_positions = distance.size();
	size_t* leaf_index =  new size_t[ n_positions ];
	double* v_dist =  new double[ n_positions ];
	for(  size_t i_pos = 0  ;  i_pos < n_positions  ; i_pos++    ){
		leaf_index[ i_pos ] = leaf_id[ i_pos ];
		v_dist[ i_pos ] = distance[ i_pos ];
	}
	leaf_id.clear();
	distance.clear();
	ShellSortSizeTDouble( leaf_index , v_dist , n_positions );

	//Actualizing loaded proportion aplying the prop_variation_
	prop_vol += prop_variation_;

	//Erasing cells in order to accomplish the total proportion
	double actual_prop = cell_prop * n_positions;
	if(  actual_prop < prop_vol  ){
		assert( 0 );
		return true;
	}

	double test_distance = 0.5 * percentage_ * 0.95;

	for(  size_t i_pos = 0  ;  i_pos < n_positions  ;  i_pos++  ){
		if(  v_dist[ i_pos ] < test_distance  ){
			DataInsideOctreeCell* data = leaves[ leaf_index[ i_pos ] ]->pGetData();
			data->_SetNotFixedElement( );
		}else{
			break;
		}
	}
	delete[] leaf_index;
	delete[] v_dist;
	return true;		 
}


//UTILITIES
/**
 *Method to read the information of the vdb file previusly opened
 */
void Mesher::ReadVdbFile(  ){
	vdb_->ReadCompleteFile( fold_  );
}

/**
 *Saving vdb mesh on Gid Spheres format
 */
void Mesher::SaveVdbOnGiDMesh(  ){
	vdb_->SaveGiDSphericalMesh(  );
}

/**
 *Scaling boundary nodes
 *This method calls to the Boundary::ScaleNodes method from Boundary class
 */
void Mesher::ScaleCapside(  ){
	vdb_->ScaleSpheresMesh( percentage_ );
	vdb_->PrintResumenOnWarnings();
}

/**
 *Refining local cells based on the refinement level given as parameter
 *This method subdivide the leaves of the local root until reach the r_level_
 */
void Mesher::RefineLocalRoot(){
	int i_level = 0;
	while( i_level < r_level_ ){
		i_level++;
		//Obtaining all leaves vector
		OctreeCell_vector leaves;
		this->GetAllLocalRootLeaves( leaves );
		size_t n_leaves = leaves.size();

		for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
			OctreeCell* leaf = leaves[ i_leaf ];
			DataInsideOctreeCell* data = leaf->pGetData();
			if(  data  ){
				if( data->_GetIsMesh() == 1 ){
					//Subdividing cell
					leaf->SubdivideCell(  );
					//Assigning empty data to childrens
		      DataInsideOctreeCell** data_children = new DataInsideOctreeCell*[ 8 ];
					for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
						data_children[ i_child ] = new DataInsideOctreeCell;
						OctreeCell* child = leaf->pGetChild( i_child );
						DataInsideOctreeCell** data_child = child->pGetDataPointer();
						(*data_child) = data_children[ i_child ]; 
					}
					//Assigning interior atoms on childs
					size_t n_atoms = data->_GetNInteriorAtoms();
					for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
						OctreeCell* child = leaf->pGetChild( i_child );
						double child_min_coord[ 3 ],child_center[ 3 ], child_radius;;
						double child_max_coord[ 3 ];
						child->GetMinPointNormalized( child_min_coord );
						child->GetMaxPointNormalized( child_max_coord );
						for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
							child_center[ i_pos ] = ( child_min_coord[ i_pos ] + child_max_coord[ i_pos ] ) * 0.5;
						}
						child_radius = ( child_max_coord[ 0 ] - child_min_coord[ 0 ] ) * 0.5;
						for(  size_t i_atom = 0  ;  i_atom < n_atoms  ;  i_atom++  ){
							Atom* atom = data->_GetInteriorAtom( i_atom );
							double coordinates[ 3 ],radius;
							coordinates[ 0 ] = atom->GetCoordinate( 0 );
							coordinates[ 1 ] = atom->GetCoordinate( 1 );
							coordinates[ 2 ] = atom->GetCoordinate( 2 );
							radius = atom->GetRadius( );
							bool flag = PointIsContainedInAABB( coordinates , child_min_coord , child_max_coord );
							if(  flag  ){
								data_children[ i_child ]->_AddInteriorAtom( atom , data->_GetInteriorProteinIndex( i_atom ) , 
																																data->_GetInteriorAminoacidIndex( i_atom ) );
							}else{
								if(  AABBvsAABB( child_center , child_radius , coordinates , radius )  ){
									//possibly exists a intersection between atom and cell, is required to perform
									//the test
									if(  SphereIntersectsCell(  child_min_coord , child_max_coord , coordinates , radius  )  ) {
										data_children[ i_child ]->_AddBoundaryAtom( atom , data->_GetInteriorProteinIndex( i_atom ) , 
																																			data->_GetInteriorAminoacidIndex( i_atom ) );
									}
								}
							}
						}
					}

					//Assigning boundary atoms on childs
					n_atoms = data->_GetNBoundaryAtoms();
					for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
						OctreeCell* child = leaf->pGetChild( i_child );
						double child_min_coord[ 3 ],child_center[ 3 ], child_radius;
						double child_max_coord[ 3 ];
						child->GetMinPointNormalized( child_min_coord );
						child->GetMaxPointNormalized( child_max_coord );
						for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
							child_center[ i_pos ] = ( child_min_coord[ i_pos ] + child_max_coord[ i_pos ] ) * 0.5;
						}
						child_radius = ( child_max_coord[ 0 ] - child_min_coord[ 0 ] ) * 0.5;
						for(  size_t i_atom = 0  ;  i_atom < n_atoms  ;  i_atom++  ){
							Atom* atom = data->_GetBoundaryAtom( i_atom );
							double coordinates[ 3 ],radius;
							radius = atom->GetRadius( );
							for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
								coordinates[ i_pos ] = atom->GetCoordinate( i_pos );
							}
							//performingAABB vs AABB test
							if(  AABBvsAABB( child_center , child_radius , coordinates , radius )  ){
								//possibly exists a intersection between atom and cell, is required to perform
								//the test
								if(  SphereIntersectsCell(  child_min_coord , child_max_coord , coordinates , radius  )  ) {
									data_children[ i_child ]->_AddBoundaryAtom( atom , data->_GetBoundaryProteinIndex( i_atom ) , 
																																		data->_GetBoundaryAminoacidIndex( i_atom ) );
								}
							}
						}
					}
				}
			}
		}
		leaves.clear();
	}
}

/**
 *Assigning the atoms belonging to the local root
 */
void Mesher::AssignAtomsOnLocalRoot( ){
	//Allocating new data structure
	DataInsideOctreeCell* data_ = new DataInsideOctreeCell;
	//testing if atoms belongs to local root or not
	size_t n_proteins = vdb_->GetNProteins();
	for(  size_t i_prot = 0  ;  i_prot < n_proteins  ;  i_prot++  ){
		Protein* protein = vdb_->GetProtein( i_prot );
		size_t n_aminoacids = protein->GetNAminoacids( );
		for(  size_t i_amino = 0  ;  i_amino < n_aminoacids  ;  i_amino++  ){
			Aminoacid* amino = protein->GetAminoacid( i_amino );
			size_t n_atoms = amino->GetNAtoms( );
			for(  size_t i_atom = 0  ;  i_atom < n_atoms  ;  i_atom++  ){
				Atom* atom = amino->GetAtom( i_atom );
				data_->_AddInteriorAtom( atom , i_prot , i_amino );		
			}
		}
	}
	//Assigning allocated data filled with atoms on the local root data pointer
	DataInsideOctreeCell** data_cell = local_root_->pGetDataPointer();
	(*data_cell) = data_;
}

/**
 *This method calculates the refinement level and the percentage of octree that need to be 
 *used based on the resolution required using the expresion 
 *    x = ( DOMAIN )/( ( 2^REFINEMENT )*RESOLUTION )
 */
void Mesher::CalculateRefinementAndPercentage(){
	//Getting the domain in the mesh
	double domain;
	domain = vdb_->CalculateMaximumDomain();

	//Starting to look for the refinement and the percentage
	bool flag = false;
	r_level_ = 0;
	while( !flag ){
		percentage_ = domain / ( ( pow( 2 , r_level_ ) )*( resolution_ ) );
		if(  percentage_ < 1.00  ){
			flag = true;
		}else{
			r_level_++;
		}
	}
}

/**
 *Selecting fold_ to be aligned with Y axis
 *@param[out] coords is the coordinate to be rotated
 */
void Mesher::RotateCoordianteToFold( double* coord ){
	switch( fold_ ){
		case 5:
			this->RotateFold5( coord );
			break;
		case 3:
			this->RotateFold3( coord );
			break;
		case 2:
			this->RotateFold2( coord );
			break;
		default:
			std::cout << " Fold is not valid, need to be 5,3 or 2  and you are using " << fold_ << std::endl;
			assert( 0 );
			break;
	}
}

/**
 *Selecting index from fold 2 to align with Y axis
 *@param[out] coord Is the coordinate to be rotated
 */
void Mesher::RotateFold2( double* coord ){
	switch( fold_index_ ){
		case 0:
			this->RotateFold2Id0( coord );
			break;
		case 1:
			this->RotateFold2Id1( coord );
			break;
		case 2:
			this->RotateFold2Id2( coord );
			break;
		case 3:
			this->RotateFold2Id3( coord );
			break;
		case 4:
			this->RotateFold2Id4( coord );
			break;
		case 5:
			this->RotateFold2Id5( coord );
			break;
		case 6:
			this->RotateFold2Id6( coord );
			break;
		case 7:
			this->RotateFold2Id7( coord );
			break;
		case 8:
			this->RotateFold2Id8( coord );
			break;
		case 9:
			this->RotateFold2Id9( coord );
			break;
		case 10:
			this->RotateFold2Id10( coord );
			break;
		case 11:
			this->RotateFold2Id11( coord );
			break;
		default:
			std::cout << " Fold id not valid " << std::endl;
			assert( 0 );
			break;
	}
}

/**
 *Selecting index from fold 3 to align with Y axis
 *@param[out] coord Is the coordinate to be rotated
 */
void Mesher::RotateFold3( double* coord ){
	switch( fold_index_ ){
		case 0:
			this->RotateFold3Id0( coord );
			break;
		case 1:
			this->RotateFold3Id1( coord );
			break;
		case 2:
			this->RotateFold3Id2( coord );
			break;
		case 3:
			this->RotateFold3Id3( coord );
			break;
		case 4:
			this->RotateFold3Id4( coord );
			break;
		case 5:
			this->RotateFold3Id5( coord );
			break;
		case 6:
			this->RotateFold3Id6( coord );
			break;
		case 7:
			this->RotateFold3Id7( coord );
			break;
		case 8:
			this->RotateFold3Id8( coord );
			break;
		case 9:
			this->RotateFold3Id9( coord );
			break;
		case 10:
			this->RotateFold3Id10( coord );
			break;
		case 11:
			this->RotateFold3Id11( coord );
			break;
		default:
			std::cout << " Fold id not valid " << std::endl;
			assert( 0 );
			break;
	}
}

/**
 *Selecting index from fold 5 to align with Y axis
 *@param[out] coord Is the coordinate to be rotated
 */
void Mesher::RotateFold5( double* coord ){
	switch( fold_index_ ){
		case 0:
			this->RotateFold5Id0( coord );
			break;
		case 1:
			this->RotateFold5Id1( coord );
			break;
		case 2:
			this->RotateFold5Id2( coord );
			break;
		case 3:
			this->RotateFold5Id3( coord );
			break;
		case 4:
			this->RotateFold5Id4( coord );
			break;
		case 5:
			this->RotateFold5Id5( coord );
			break;
		case 6:
			this->RotateFold5Id6( coord );
			break;
		case 7:
			this->RotateFold5Id7( coord );
			break;
		case 8:
			this->RotateFold5Id8( coord );
			break;
		case 9:
			this->RotateFold5Id9( coord );
			break;
		case 10:
			this->RotateFold5Id10( coord );
			break;
		case 11:
			this->RotateFold5Id11( coord );
			break;
		default:
			std::cout << " Fold id not valid " << std::endl;
			assert( 0 );
			break;
	}
}

/**
 *Rotating atom coordinate to put id 0 from fold 2 aligned with Y axis
 *The coordinate is rotated over the X axis -90 degrees (1.570796327 radians)
 */
void Mesher::RotateFold2Id0( double* coord ){
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	double newcoord[ 3 ] , teta;
	teta  = -1.570796327;
	matrix[ 0 ][ 0 ] = 1.0; matrix[ 0 ][ 1 ] = 0.0;         matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = 0.0; matrix[ 1 ][ 1 ] = cos( teta ); matrix[ 1 ][ 2 ] = -sin( teta );
	matrix[ 2 ][ 0 ] = 0.0; matrix[ 2 ][ 1 ] = sin( teta ); matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	coord[ 0 ] = newcoord[ 0 ];	coord[ 1 ] = newcoord[ 1 ];	coord[ 2 ] = newcoord[ 2 ];
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put id 1 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id1( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 2 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id2( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 3 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id3( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 4 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id4( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 5 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id5( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 6 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id6( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 7 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id7( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 8 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id8( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 9 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id9( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 10 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id10( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 11 from fold 2 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold2Id11( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 0 from fold 3 aligned with Y axis
 *The coordinate is rotated over the Y axis -20.905192903 degrees (0.364864447 radians)
 *The coordinate is rotated over the X axis -90 degrees (1.570796327 radians)
 */
void Mesher::RotateFold3Id0( double* coord ){
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	//Y rotation
	double newcoord[ 3 ] , teta;
	teta  = -0.364864447;
	matrix[ 0 ][ 0 ] = cos( teta );  matrix[ 0 ][ 1 ] = 0.0; matrix[ 0 ][ 2 ] = sin( teta ); 
	matrix[ 1 ][ 0 ] = 0.0;          matrix[ 1 ][ 1 ] = 1.0; matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = -sin( teta ); matrix[ 2 ][ 1 ] = 0.0; matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	//X Rotation
	teta  = -1.570796327;
	matrix[ 0 ][ 0 ] = 1.0; matrix[ 0 ][ 1 ] = 0.0; matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = 0.0; matrix[ 1 ][ 1 ] = cos( teta ); matrix[ 1 ][ 2 ] = -sin( teta );
	matrix[ 2 ][ 0 ] = 0.0; matrix[ 2 ][ 1 ] = sin( teta ); matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , newcoord , 3 , coord , 3 );
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put id 1 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id1( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 2 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id2( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 3 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id3( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 4 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id4( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 5 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id5( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 6 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id6( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 7 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id7( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 8 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id8( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 9 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id9( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 10 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id10( double* coord ){
	assert( 0 );//NOT DEFINED YET
}

/**
 *Rotating atom coordinate to put id 11 from fold 3 aligned with Y axis
 *?????????????
 */
void Mesher::RotateFold3Id11( double* coord ){
	assert( 0 );//NOT DEFINED YET
}


/**
 *Rotating atom coordinate to put the fold 0 aligned with Y axis
 *The coordinate is rotated over the X axis 31.71745718-90 degrees (1.017222269 radians)
 */
void Mesher::RotateFold5Id0( double* coord ){
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	double newcoord[ 3 ] , teta;
	teta  = -1.017222269;
	matrix[ 0 ][ 0 ] = 1.0; matrix[ 0 ][ 1 ] = 0.0;         matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = 0.0; matrix[ 1 ][ 1 ] = cos( teta ); matrix[ 1 ][ 2 ] = -sin( teta );
	matrix[ 2 ][ 0 ] = 0.0; matrix[ 2 ][ 1 ] = sin( teta ); matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	coord[ 0 ] = newcoord[ 0 ];	coord[ 1 ] = newcoord[ 1 ];	coord[ 2 ] = newcoord[ 2 ];
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 1 aligned with Y axis
 *The coordinate is rotated over the X axis 121.71745718 degrees (2.124370385 radians)
 */
void Mesher::RotateFold5Id1( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	double newcoord[ 3 ] , teta;
	teta  = 2.124370385;
	matrix[ 0 ][ 0 ] = 1.0; matrix[ 0 ][ 1 ] = 0.0;         matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = 0.0; matrix[ 1 ][ 1 ] = cos( teta ); matrix[ 1 ][ 2 ] = -sin( teta );
	matrix[ 2 ][ 0 ] = 0.0; matrix[ 2 ][ 1 ] = sin( teta ); matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	coord[ 0 ] = newcoord[ 0 ];	coord[ 1 ] = newcoord[ 1 ];	coord[ 2 ] = newcoord[ 2 ];
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 2 aligned with Y axis
 *The coordinate is rotated over the Z axis 31.717252354 degrees (0.553570483 radians)
 */
void Mesher::RotateFold5Id2( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	double newcoord[ 3 ] , teta;
	teta  = 0.553570483;
	matrix[ 0 ][ 0 ] = cos( teta ); matrix[ 0 ][ 1 ] = -sin( teta ); matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = sin( teta ); matrix[ 1 ][ 1 ] = cos( teta );  matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = 0.0;         matrix[ 2 ][ 1 ] = 0.0;          matrix[ 2 ][ 2 ] = 1.0;
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	coord[ 0 ] = newcoord[ 0 ];	coord[ 1 ] = newcoord[ 1 ];	coord[ 2 ] = newcoord[ 2 ];
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 3 aligned with Y axis
 *The coordinate is rotated over the Z axis -31.717252354 degrees (-0.553570483 radians)
 */
void Mesher::RotateFold5Id3( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	double newcoord[ 3 ] , teta;
	teta  = -0.553570483;
	matrix[ 0 ][ 0 ] = cos( teta ); matrix[ 0 ][ 1 ] = -sin( teta ); matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = sin( teta ); matrix[ 1 ][ 1 ] = cos( teta );  matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = 0.0;         matrix[ 2 ][ 1 ] = 0.0;          matrix[ 2 ][ 2 ] = 1.0;
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	coord[ 0 ] = newcoord[ 0 ];	coord[ 1 ] = newcoord[ 1 ];	coord[ 2 ] = newcoord[ 2 ];
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 4 aligned with Y axis
 *The coordinate is rotated over the Y axis -58.282747642 degrees (-1.017225843 radians)
 *The coordinate is rotated over the X axis 90 degrees (1.570796327 radians)
 */
void Mesher::RotateFold5Id4( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	//Y rotation
	double newcoord[ 3 ] , teta;
	teta  = -1.017225843;
	matrix[ 0 ][ 0 ] = cos( teta );  matrix[ 0 ][ 1 ] = 0.0; matrix[ 0 ][ 2 ] = sin( teta ); 
	matrix[ 1 ][ 0 ] = 0.0;          matrix[ 1 ][ 1 ] = 1.0; matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = -sin( teta ); matrix[ 2 ][ 1 ] = 0.0; matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	//X Rotation
	teta  = 1.570796327;
	matrix[ 0 ][ 0 ] = 1.0; matrix[ 0 ][ 1 ] = 0.0; matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = 0.0; matrix[ 1 ][ 1 ] = cos( teta ); matrix[ 1 ][ 2 ] = -sin( teta );
	matrix[ 2 ][ 0 ] = 0.0; matrix[ 2 ][ 1 ] = sin( teta ); matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , newcoord , 3 , coord , 3 );
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 5 aligned with Y axis
 *The coordinate is rotated over the Y axis 58.28254282 degrees (1.017222269 radians)
 *The coordinate is rotated over the X axis 90 degrees (1.570796327 radians)
 */
void Mesher::RotateFold5Id5( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	//Y rotation
	double newcoord[ 3 ] , teta;
	teta  = 1.017222269;
	matrix[ 0 ][ 0 ] = cos( teta );  matrix[ 0 ][ 1 ] = 0.0; matrix[ 0 ][ 2 ] = sin( teta ); 
	matrix[ 1 ][ 0 ] = 0.0;          matrix[ 1 ][ 1 ] = 1.0; matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = -sin( teta ); matrix[ 2 ][ 1 ] = 0.0; matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	//X Rotation
	teta  = 1.570796327;
	matrix[ 0 ][ 0 ] = 1.0; matrix[ 0 ][ 1 ] = 0.0; matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = 0.0; matrix[ 1 ][ 1 ] = cos( teta ); matrix[ 1 ][ 2 ] = -sin( teta );
	matrix[ 2 ][ 0 ] = 0.0; matrix[ 2 ][ 1 ] = sin( teta ); matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , newcoord , 3 , coord , 3 );
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 6 aligned with Y axis
 *The coordinate is rotated over the Z axis 148.282747 degrees (2.588022159 radians)
 */
void Mesher::RotateFold5Id6( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	double newcoord[ 3 ] , teta;
	teta  = 2.588022159;
	matrix[ 0 ][ 0 ] = cos( teta ); matrix[ 0 ][ 1 ] = -sin( teta ); matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = sin( teta ); matrix[ 1 ][ 1 ] = cos( teta );  matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = 0.0;         matrix[ 2 ][ 1 ] = 0.0;          matrix[ 2 ][ 2 ] = 1.0;
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	coord[ 0 ] = newcoord[ 0 ];	coord[ 1 ] = newcoord[ 1 ];	coord[ 2 ] = newcoord[ 2 ];
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 7 aligned with Y axis
 *The coordinate is rotated over the Z axis -148.282747 degrees (-2.588022159 radians)
 */
void Mesher::RotateFold5Id7( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	double newcoord[ 3 ] , teta;
	teta  = -2.588022159;
	matrix[ 0 ][ 0 ] = cos( teta ); matrix[ 0 ][ 1 ] = -sin( teta ); matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = sin( teta ); matrix[ 1 ][ 1 ] = cos( teta );  matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = 0.0;         matrix[ 2 ][ 1 ] = 0.0;          matrix[ 2 ][ 2 ] = 1.0;
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	coord[ 0 ] = newcoord[ 0 ];	coord[ 1 ] = newcoord[ 1 ];	coord[ 2 ] = newcoord[ 2 ];
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 8 aligned with Y axis
 *The coordinate is rotated over the Y axis -31.717252358 degrees (-0.553570483 radians)
 *The coordinate is rotated over the Z axis 90 degrees (1.570796327 radians)
 */
void Mesher::RotateFold5Id8( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	//Y rotation
	double newcoord[ 3 ] , teta;
	teta  = -0.553570483;
	matrix[ 0 ][ 0 ] = cos( teta );  matrix[ 0 ][ 1 ] = 0.0; matrix[ 0 ][ 2 ] = sin( teta ); 
	matrix[ 1 ][ 0 ] = 0.0;          matrix[ 1 ][ 1 ] = 1.0; matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = -sin( teta ); matrix[ 2 ][ 1 ] = 0.0; matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	//Z Rotation
	teta  = 1.570796327;
	matrix[ 0 ][ 0 ] = cos( teta ); matrix[ 0 ][ 1 ] = -sin( teta ); matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = sin( teta ); matrix[ 1 ][ 1 ] = cos( teta );  matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = 0.0;         matrix[ 2 ][ 1 ] = 0.0;          matrix[ 2 ][ 2 ] = 1.0;
	MatrixVectorMultiplication( matrix , 3 , 3 , newcoord , 3 , coord , 3 );
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 9 aligned with Y axis
 *The coordinate is rotated over the X axis -58.282747642 degrees (-1.017225843 radians)
 */
void Mesher::RotateFold5Id9( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	double newcoord[ 3 ] , teta;
	teta  = -1.017225843;
	matrix[ 0 ][ 0 ] = 1.0; matrix[ 0 ][ 1 ] = 0.0;         matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = 0.0; matrix[ 1 ][ 1 ] = cos( teta ); matrix[ 1 ][ 2 ] = -sin( teta );
	matrix[ 2 ][ 0 ] = 0.0; matrix[ 2 ][ 1 ] = sin( teta ); matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	coord[ 0 ] = newcoord[ 0 ];	coord[ 1 ] = newcoord[ 1 ];	coord[ 2 ] = newcoord[ 2 ];
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 10 aligned with Y axis
 *The coordinate is rotated over the X axis -121.71745718 degrees (-2.124370385 radians)
 */
void Mesher::RotateFold5Id10( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	double newcoord[ 3 ] , teta;
	teta  = -2.124370385;
	matrix[ 0 ][ 0 ] = 1.0; matrix[ 0 ][ 1 ] = 0.0;         matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = 0.0; matrix[ 1 ][ 1 ] = cos( teta ); matrix[ 1 ][ 2 ] = -sin( teta );
	matrix[ 2 ][ 0 ] = 0.0; matrix[ 2 ][ 1 ] = sin( teta ); matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	coord[ 0 ] = newcoord[ 0 ];	coord[ 1 ] = newcoord[ 1 ];	coord[ 2 ] = newcoord[ 2 ];
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Rotating atom coordinate to put the fold 11 aligned with Y axis
 *The coordinate is rotated over the Y axis 31.71745718 degrees (0.553574058 radians)
 *The coordinate is rotated over the Z axis -90 degrees (-1.570796327 radians)
 */
void Mesher::RotateFold5Id11( double* coord ){
	assert( 0 );//NEED TO BE CORRECTED
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	//Y rotation
	double newcoord[ 3 ] , teta;
	teta  = 0.553574058;
	matrix[ 0 ][ 0 ] = cos( teta );  matrix[ 0 ][ 1 ] = 0.0; matrix[ 0 ][ 2 ] = sin( teta ); 
	matrix[ 1 ][ 0 ] = 0.0;          matrix[ 1 ][ 1 ] = 1.0; matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = -sin( teta ); matrix[ 2 ][ 1 ] = 0.0; matrix[ 2 ][ 2 ] = cos( teta );
	MatrixVectorMultiplication( matrix , 3 , 3 , coord , 3 , newcoord , 3 );
	//Z Rotation
	teta  = -1.570796327;
	matrix[ 0 ][ 0 ] = cos( teta ); matrix[ 0 ][ 1 ] = -sin( teta ); matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = sin( teta ); matrix[ 1 ][ 1 ] = cos( teta );  matrix[ 1 ][ 2 ] = 0.0;
	matrix[ 2 ][ 0 ] = 0.0;         matrix[ 2 ][ 1 ] = 0.0;          matrix[ 2 ][ 2 ] = 1.0;
	MatrixVectorMultiplication( matrix , 3 , 3 , newcoord , 3 , coord , 3 );
	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;
}

/**
 *Method to calculate total volume of mesh
 *Mesh cells are the cells that in the DataInsideOctreeCell contain the is_mesh_ variable
 * as TRUE
 *@return A double value with the total volume of the mesh
 */
double Mesher::CalculateOctreeMeshVolume(  double* cell_volume  ){
	double tot_vol;
	size_t counter = 0;
	double cell_size;
	bool flag = false;
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	size_t n_leaves = leaves.size();
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		OctreeCell* leaf = leaves[ i_leaf ];
		DataInsideOctreeCell* data = leaf->pGetData();		
		if(  ( data->_GetIsMesh( ) == 1 )  && !flag  ){
			flag = true;
			double min_point[ 3 ];
			double max_point[ 3 ];
			leaf->GetMinPointNormalized( min_point ); 
			leaf->GetMaxPointNormalized( max_point );
			cell_size = max_point[ 0 ] - min_point[ 0 ];
		}
		if(  data->_GetIsMesh( ) == 1  ){
			counter++;
		}
	}
	(*cell_volume) = cell_size * cell_size * cell_size;
	tot_vol = counter * (*cell_volume);
	return tot_vol;
}

/**
 *This method calculates the proportion of volume to be loaded based on the virus, fold
 *and mesh resolucion, the method uses a cuadrati model using cuadratic regresion
 *@return A double value with the proportion required
 */
double Mesher::CalculateProportionLoaded(  ){
	double vol_prop;
	switch( fold_ ){
		case 2:
			if(  !strcmp( virus_ , "1CWP" )  ){
				vol_prop = 0.2565350 - 0.0171050 * resolution_capside_ + 0.004106000000 * ( resolution_capside_ * resolution_capside_ );
			}
			if(  !strcmp( virus_ , "4G93" )  ){
				vol_prop = 0.0780206 + 0.0180684 * resolution_capside_ + 0.000386900000 * ( resolution_capside_ * resolution_capside_ );
			}
			if(  !strcmp( virus_ , "3IZG" )  ){
				vol_prop = 0.0413656 + 0.0204383 * resolution_capside_ - 0.000017227800 * ( resolution_capside_ * resolution_capside_ );
			}
			if(  !strcmp( virus_ , "3J4U" )  ){
				vol_prop = 0.0174092 + 0.0002131 * resolution_capside_ + 0.000784200000 * ( resolution_capside_ * resolution_capside_ );
			}
			break;
		case 3:
			if(  !strcmp( virus_ , "1CWP" )  ){
				vol_prop = 0.4107155 + 0.0246620 * resolution_capside_ + 0.000856200000 * ( resolution_capside_ * resolution_capside_ );
			}
			if(  !strcmp( virus_ , "4G93" )  ){
				vol_prop = 0.2106251 + 0.0229389 * resolution_capside_ - 0.000361400000 * ( resolution_capside_ * resolution_capside_ );
			}
			if(  !strcmp( virus_ , "3IZG" )  ){
				vol_prop = 0.0225870 + 0.0238720 * resolution_capside_ - 0.000171000000 * ( resolution_capside_ * resolution_capside_ );
			}
			if(  !strcmp( virus_ , "3J4U" )  ){
				vol_prop = 0.0076560 - 0.0023530 * resolution_capside_ + 0.000500000000 * ( resolution_capside_ * resolution_capside_ );
			}
			break;
		case 5:
			if(  !strcmp( virus_ , "1CWP" )  ){
				vol_prop = 0.5303590 + 0.0297540 * resolution_capside_ + 0.000409000000 * ( resolution_capside_ * resolution_capside_ );
			}
			if(  !strcmp( virus_ , "4G93" )  ){
				vol_prop = 0.1123588 + 0.0340842 * resolution_capside_ - 0.000005522564 * ( resolution_capside_ * resolution_capside_ );
			}
			if(  !strcmp( virus_ , "3IZG" )  ){
				vol_prop = 0.2403398 + 0.0194864 * resolution_capside_ - 0.000064931300 * ( resolution_capside_ * resolution_capside_ );
			}
			if(  !strcmp( virus_ , "3J4U" )  ){
				vol_prop = 0.4230595 + 0.0228704 * resolution_capside_ + 0.000303800000 * ( resolution_capside_ * resolution_capside_ );
			}
			break;
	}
	vol_prop /= 100.000;
	return vol_prop;
}

/**
 *Calculating density of loaded elements
 *The method obtains the total volume of the loaded elements and uses the masa to obtain
 *The elements density
 *@param[in] masa Is the masa of the elements loaded
 *@return A double value with the elements density 
 */
double Mesher::CalcDensityForLoadedElements( double masa ){
	size_t n_loaded = 0;
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	size_t n_leaves = leaves.size();
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		OctreeCell* leaf = leaves[ i_leaf ];
		DataInsideOctreeCell* data = leaf->pGetData();
		if( data->_GetIsMesh() == 1 ){
			if(  data->_IsLoaded(  )  ){
				n_loaded++;
			}
		}
	}

	return ( masa / ( resolution_ * resolution_ * resolution_ * n_loaded ) );
}

/**
 *This method interpolate the resolution of the CapsidMesh mesher to the equivalent
 *resolution in the octree mesher
 *@return A double value with the real resolutuon on the octree mesher based on the 
 *        CapsidMesh mesher resolution
 */
double Mesher::InterpolateCapsidResolution(  ){
	double real_resolution;
	if(  !strcmp( virus_ , "1CWP" )  ){
		real_resolution = this->InterpolateCapsidResolution1CWP(  );
	}
	if(  !strcmp( virus_ , "4G93" )  ){
		real_resolution = this->InterpolateCapsidResolution4G93(  );
	}
	if(  !strcmp( virus_ , "3IZG" )  ){
		real_resolution = this->InterpolateCapsidResolution3IZG(  );
	}
	if(  !strcmp( virus_ , "3J4U" )  ){
		real_resolution = this->InterpolateCapsidResolution3J4U(  );
	}
	return real_resolution;
}

/**
 *Interpolate resolution for 1CWP
 *@return A double value witht the real resolution
 */
double Mesher::InterpolateCapsidResolution1CWP(  ){
	if(  fabs( resolution_capside_ - 16.00 ) < 1e-5  ){
		return 10.90;
	}
	if(  fabs( resolution_capside_ - 15.00 ) < 1e-5  ){
		return 10.90;
	}
	if(  fabs( resolution_capside_ - 14.00 ) < 1e-5  ){
		return 11.495;
	}
	if(  fabs( resolution_capside_ - 13.00 ) < 1e-5  ){
		return 6.45;
	}
	if(  fabs( resolution_capside_ - 12.00 ) < 1e-5  ){
		return 6.45;
	}
	if(  fabs( resolution_capside_ - 11.00 ) < 1e-5  ){
		return 6.22;
	}
	if(  fabs( resolution_capside_ - 10.00 ) < 1e-5  ){
		return 5.48;
	}
	if(  fabs( resolution_capside_ -  9.00 ) < 1e-5  ){
		return 5.21;
	}
	if(  fabs( resolution_capside_ -  8.00 ) < 1e-5  ){
		return 5.275;
	}
	if(  fabs( resolution_capside_ -  7.00 ) < 1e-5  ){
		return 4.51;
	}
	if(  fabs( resolution_capside_ -  6.00 ) < 1e-5  ){
		return 2.00;
	}
	if(  fabs( resolution_capside_ -  5.00 ) < 1e-5  ){
		return 1.08;
	}
	if(  fabs( resolution_capside_ -  4.00 ) < 1e-5  ){
		return 0.935;
	}
	if(  fabs( resolution_capside_ -  3.00 ) < 1e-5  ){
		return 0.825;
	}
	if(  fabs( resolution_capside_ -  2.00 ) < 1e-5  ){
		return 0.71;
	}
	if(  fabs( resolution_capside_ -  1.00 ) < 1e-5  ){
		return 0.585;
	}
	if(  fabs( resolution_capside_ -  0.95 ) < 1e-5  ){
		return 0.576;
	}
	if(  fabs( resolution_capside_ -  0.90 ) < 1e-5  ){
		return 0.569;
	}
	return 0.5615;
}

/**
 *Interpolate resolution for 4G93
 *@return A double value witht the real resolution
 */
double Mesher::InterpolateCapsidResolution4G93(  ){
	if(  fabs( resolution_capside_ - 16.00 ) < 1e-5  ){
		return 7.32;
	}
	if(  fabs( resolution_capside_ - 15.00 ) < 1e-5  ){
		return 7.32;
	}
	if(  fabs( resolution_capside_ - 14.00 ) < 1e-5  ){
		return 7.76;
	}
	if(  fabs( resolution_capside_ - 13.00 ) < 1e-5  ){
		return 7.47;
	}
	if(  fabs( resolution_capside_ - 12.00 ) < 1e-5  ){
		return 7.2;
	}
	if(  fabs( resolution_capside_ - 11.00 ) < 1e-5  ){
		return 6.01;
	}
	if(  fabs( resolution_capside_ - 10.00 ) < 1e-5  ){
		return 6.425;
	}
	if(  fabs( resolution_capside_ -  9.00 ) < 1e-5  ){
		return 6.35;
	}
	if(  fabs( resolution_capside_ -  8.00 ) < 1e-5  ){
		return 3.965;
	}
	if(  fabs( resolution_capside_ -  7.00 ) < 1e-5  ){
		return 3.71;
	}
	if(  fabs( resolution_capside_ -  6.00 ) < 1e-5  ){
		return 2.352;
	}
	if(  fabs( resolution_capside_ -  5.00 ) < 1e-5  ){
		return 1.27;
	}
	if(  fabs( resolution_capside_ -  4.00 ) < 1e-5  ){
		return 1.14;
	}
	if(  fabs( resolution_capside_ -  3.00 ) < 1e-5  ){
		return 1.025;
	}
	if(  fabs( resolution_capside_ -  2.00 ) < 1e-5  ){
		return 0.87;
	}
	if(  fabs( resolution_capside_ -  1.00 ) < 1e-5  ){
		return 0.7;
	}
	if(  fabs( resolution_capside_ -  0.95 ) < 1e-5  ){
		return 0.696;
	}
	if(  fabs( resolution_capside_ -  0.90 ) < 1e-5  ){
		return 0.687;
	}
	return 0.678;
}

/**
 *Interpolate resolution for 3IZG
 *@return A double value witht the real resolution
 */
double Mesher::InterpolateCapsidResolution3IZG(  ){
	if(  fabs( resolution_capside_ - 16.00 ) < 1e-5  ){
		return 9.5;
	}
	if(  fabs( resolution_capside_ - 15.00 ) < 1e-5  ){
		return 9.5;
	}
	if(  fabs( resolution_capside_ - 14.00 ) < 1e-5  ){
		return 9.5;
	}
	if(  fabs( resolution_capside_ - 13.00 ) < 1e-5  ){
		return 9.11;
	}
	if(  fabs( resolution_capside_ - 12.00 ) < 1e-5  ){
		return 8.60;
	}
	if(  fabs( resolution_capside_ - 11.00 ) < 1e-5  ){
		return 8.60;
	}
	if(  fabs( resolution_capside_ - 10.00 ) < 1e-5  ){
		return 8.429;
	}
	if(  fabs( resolution_capside_ -  9.00 ) < 1e-5  ){
		return 6.00;
	}
	if(  fabs( resolution_capside_ -  8.00 ) < 1e-5  ){
		return 3.99;
	}
	if(  fabs( resolution_capside_ -  7.00 ) < 1e-5  ){
		return 3.3;
	}
	if(  fabs( resolution_capside_ -  6.00 ) < 1e-5  ){
		return 3.35;
	}
	if(  fabs( resolution_capside_ -  5.00 ) < 1e-5  ){
		return 2.03;
	}
	if(  fabs( resolution_capside_ -  4.00 ) < 1e-5  ){
		return 1.80;
	}
	if(  fabs( resolution_capside_ -  3.00 ) < 1e-5  ){
		return 1.58;
	}
	if(  fabs( resolution_capside_ -  2.00 ) < 1e-5  ){
		return 1.335;
	}
	if(  fabs( resolution_capside_ -  1.00 ) < 1e-5  ){
		return 1.05;
	}
	if(  fabs( resolution_capside_ -  0.95 ) < 1e-5  ){
		return 1.039;
	}
	if(  fabs( resolution_capside_ -  0.90 ) < 1e-5  ){
		return 0.599;
	}
	return 0.59;
}

/**
 *Interpolate resolution for 3J4U
 *@return A double value witht the real resolution
 */
double Mesher::InterpolateCapsidResolution3J4U(  ){
	if(  fabs( resolution_capside_ - 16.00 ) < 1e-5  ){
		return 12.00;
	}
	if(  fabs( resolution_capside_ - 15.00 ) < 1e-5  ){
		return 12.00;
	}
	if(  fabs( resolution_capside_ - 14.00 ) < 1e-5  ){
		return 12.00;
	}
	if(  fabs( resolution_capside_ - 13.00 ) < 1e-5  ){
		return 10.3;
	}
	if(  fabs( resolution_capside_ - 12.00 ) < 1e-5  ){
		return 8.15;
	}
	if(  fabs( resolution_capside_ - 11.00 ) < 1e-5  ){
		return 6.89;
	}
	if(  fabs( resolution_capside_ - 10.00 ) < 1e-5  ){
		return 7.12;
	}
	if(  fabs( resolution_capside_ -  9.00 ) < 1e-5  ){
		return 6.85;
	}
	if(  fabs( resolution_capside_ -  8.00 ) < 1e-5  ){
		return 4.465;
	}
	if(  fabs( resolution_capside_ -  7.00 ) < 1e-5  ){
		return 4.465;
	}
	if(  fabs( resolution_capside_ -  6.00 ) < 1e-5  ){
		return 4.105;
	}
	if(  fabs( resolution_capside_ -  5.00 ) < 1e-5  ){
		return 2.35;
	}
	if(  fabs( resolution_capside_ -  4.00 ) < 1e-5  ){
		return 2.15;
	}
	if(  fabs( resolution_capside_ -  3.00 ) < 1e-5  ){
		return 1.85;
	}
	if(  fabs( resolution_capside_ -  2.00 ) < 1e-5  ){
		return 1.587;
	}
	if(  fabs( resolution_capside_ -  1.00 ) < 1e-5  ){
		return 1.24;
	}
	if(  fabs( resolution_capside_ -  0.95 ) < 1e-5  ){
		return 0.720;
	}
	if(  fabs( resolution_capside_ -  0.90 ) < 1e-5  ){
		return 0.719;
	}
	return 0.708;
}




//SAVING ON FILE
/**
 *Saving problem data in order to be solved by using FEMT
 *The method saves three files which are:
 *solver.dat file: Contains the solver information
 *problem.dat: Contains the problem information, such as the files to be saved, materials
 *             and boundary conditions
 *geometry.dat file: Contains the geometry information, number of nodes and elements
 */
void Mesher::PrintDataFilesForFEMT(  ){
	int vir_type = vdb_->GetVirusType(  ); 
  std::string name = std::to_string( fold_ );
	std::string aux = std::to_string( vir_type );
	name = name + "_T" + aux + "_" + virus_ + ".";
	//Saving solver data file
	char solver_name[ 10 ] = "CG";
	this->PrintSolverDataFileForFEMT( solver_name , 1 , 1e-5 , 1000000 , 1 , name );
	this->PrintProblemDataFileForFEMT( name );
	this->PrintGeometryDataFileForFEMT( name );
}

/**
 *Method to save the solver data for the simulation
 *@param[in] Solver Is a char indicating the solver to be used, it can be CG for conjugate
 *           gradiente, other for cholesky
 *@param[in] threads Is the amount of parallel threads to use in the solver
 *@param[in] tol  Is the tolerance to use in CG
 *@param[in] max_steps Maximum steps for CG
 *@param[in] preconditioner Is the preconditioner for CG
 *@param[in] name Is the name of the file to be saved
 */
void Mesher::PrintSolverDataFileForFEMT( char* solver , int threads , double tol , 
																		 					int max_steps , int preconditioner , 
																							std::string name ){
	name = name + "solver.dat";
	ofstream fp( name );
	fp << "; Solver file" << endl << endl;
	fp << "{Solver}" << endl;
	if(  !strcmp( solver , "CG" )  ){
		fp << "1 ; Type (1=Conjugate_gradient, 2=Cholesky_decomposition, 3=Cholesky2_decomposition, 4=Biconjugate_gradient, 5=LU_decomposition)" << endl;
	}else{
		fp << "2 ; Type (1=Conjugate_gradient, 2=Cholesky_decomposition, 3=Cholesky2_decomposition, 4=Biconjugate_gradient, 5=LU_decomposition)" << endl;
	}
	fp << threads << " ; Threads" << endl;
	fp << "; Parameters for iterative solvers" << endl; 
	fp << tol << " ; Tolerance" << endl;
	fp << max_steps << " ; Max_steps" << endl;
	fp << preconditioner << " ; Preconditioner (0=None, 1=Jacobi, 2=Incomplete_Cholesky, 3=Incomplete_Cholesky2, 4=Incomplete_LU, 5=Sparse_Approximate_Inverse)" << endl;
	fp << "0 ; Preconditioner_level" << endl; 
	fp << "0 ; Preconditioner_threshold" << endl << endl;
	fp << "{Substructuring}" << endl;
	fp << "1 ; Substructuring_threads" << endl;
	fp << "0.0001 ; Substructuring_tolerance" << endl;
	fp << "0 ; Multiple_results (0=No, 1=Yes)" << endl;
	fp.close();
}

/**
 *Method to save the problem data file to solve with FEMT
 *@param[in] name Is the name of the file to be saved
 */
void Mesher::PrintProblemDataFileForFEMT( std::string name ){
	//Getting all information to generate the file
	double density = this->CalcDensityForLoadedElements( 0.989 );
	std::vector<int> nodes_fixed;
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	size_t n_leaves = leaves.size();
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		OctreeCell* cell = leaves[ i_leaf ];
		DataInsideOctreeCell* data = cell->pGetData(  );
		if(  data->_GetIsMesh() == 1 ){
			if(  data->_IsFixed(  )  ){
				for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
					if(  !data->_NodeFixedWasPrinted(  i_node  )  ){
						nodes_fixed.push_back( data->_GetIndex( i_node )+1 );
						this->SetPrintedOnFixedNodeAndNeighbours( i_node , cell );
					}
				}
			}
		}
	}
	//Fixed nodes are sorted in order to be printed in order
	sort(nodes_fixed.begin(), nodes_fixed.end());

	//Generating problem data file
	name = name + "problem.dat";
	ofstream fp( name );
	fp << "; Problem file" << endl << endl;

	fp << "{General}" << endl;
	fp << "1 ; Plane_problem (1=PlaneStress, 2=PlaneStrain)" << endl;
	fp << "1 ; Calculate_displacement (0=No, 1=Yes)" << endl;
	fp << "1 ; Calculate_strain (0=No, 1=Yes)" << endl;
	fp << "1 ; Calculate_stress (0=No, 1=Yes)" << endl;
	fp << "1 ; Calculate_Von_Mises (0=No, 1=Yes)" << endl;
	fp << "1 ; Problem_type (1=Stationary, 2=Dynamic)" << endl;
	fp << "1 ; Save_mesh (0=No, 1=Yes)" << endl;
	fp << "0 ; Save_system_of_equations (0=No, 1=Yes)" << endl << endl;

	fp << "{DynamicParameters}" << endl;
	fp << "0.1 ; Time_per_step" << endl;
	fp << "20 ; Steps" << endl;
	fp << "4 ; Result_every_steps" << endl;
	fp << "0 ; Time_scheme_factor (The alpha-mehtod is second-order accurate and unconditionally stable for [-1/3, 0])" << endl;
	fp << "0 ; Rayleigh damping coefficient a" << endl;
	fp << "0 ; Rayleigh damping coefficient b" << endl << endl;

	fp << "{Gravity}" << endl;
	fp << "1 ; Use_mass_forces (0=No, 1=Yes)" << endl;
	fp << "9.80665 ; Gravity" << endl << endl;

	fp << "{UserFunctions}" << endl;
	fp << "; f1, f2, ..., f20 (one per line)" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl;
	fp << "0" << endl << endl;

	fp << "{Materials}" << endl;
	fp << "2 ; Materials count" << endl; 
	fp << "; Poisson_ratio Young_modulus Density Thickness" << endl;
	fp << "0.3 " << young_modulus_ << " 0 1 ; Chain A" << endl;
	fp << "0.3 " << young_modulus_ << " " <<  density << " 1 ; Loaded" << endl << endl;

	size_t n_fixed = nodes_fixed.size();

	fp << "{Displacement}" << endl;
	fp << n_fixed << " ; Displacements count" << endl;
	fp << "; Node User_function_for_x Displacement_x Fixed_x User_function_for_y Displacement_y Fixed_y User_function_for_z Displacement_z Fixed_z" << endl;
	for(  size_t i_node = 0  ;  i_node < n_fixed  ;  i_node++  ){
		fp << nodes_fixed[ i_node ] << " 0 0 1 0 0 1 0 0 1" << endl;
	}
	fp << endl;

	fp << "{NodalForce}" << endl;
	fp << "0 ; Nodal forces count" << endl;
	fp << "; Node User_function_for_x Force_x User_function_for_y Force_y User_function_for_z Force_z" << endl << endl;

	fp << "{NormalForce}" << endl;
	fp << "0 ; Normal forces count" << endl;
	fp.close();
}

/**
 *Method to save the geometry data file to solve with FEMT
 *@param[in] name Is the name of the file to be saved
 */
void Mesher::PrintGeometryDataFileForFEMT( std::string name ){
	//Obtaining the information to save the mesh
	size_t n_nodes = this->GetNNodesInFinalMesh(  );
	size_t n_elements = this->GetNElementsInFinalMesh(  );

	//Generating problem data file
	name = name + "geometry.dat";
	ofstream fp( name );

	fp << "; Geometry file" << endl << endl;
	fp << "{Nodes}" << endl;
	fp << "3 ; Dimension" << endl;
	fp << n_nodes << " ; Nodes count" << endl;
	fp << "; X1 X2 ..." << endl;
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	for(  size_t i = 0  ;  i < leaves.size()  ;  i++  ){
		OctreeCell* cell = leaves[ i ];
		DataInsideOctreeCell* data = cell->pGetData(  );
		if(  data->_GetIsMesh() == 1  ){
			double min_point[ 3 ];
			double max_point[ 3 ];
			cell->GetMinPointNormalized( min_point ); 
			cell->GetMaxPointNormalized( max_point );
			double cell_size = max_point[ 0 ] - min_point[ 0 ];
			for(  size_t j = 0  ;  j < 2  ;  j++  ){
				for(  size_t k = 0  ;  k < 2  ;  k++  ){
					for(  size_t h = 0  ;  h < 2  ;  h++  ){
						int node;
						double coord[ 3 ];
						if(  j == 0  &&  k == 0  &&  h == 0  ){
							node = 0;
							coord[ 0 ] = min_point[ 0 ];
							coord[ 1 ] = min_point[ 1 ];
							coord[ 2 ] = min_point[ 2 ];
						}
						if(  j == 0  &&  k == 0  &&  h == 1  ){
							node = 1;
							coord[ 0 ] = min_point[ 0 ] + cell_size;
							coord[ 1 ] = min_point[ 1 ];
							coord[ 2 ] = min_point[ 2 ];
						}
						if(  j == 0  &&  k == 1  &&  h == 0  ){
							node = 2;
							coord[ 0 ] = min_point[ 0 ] + cell_size;
							coord[ 1 ] = min_point[ 1 ] + cell_size;
							coord[ 2 ] = min_point[ 2 ];
						}
						if(  j == 0  &&  k == 1  &&  h == 1  ){
							node = 3;
							coord[ 0 ] = min_point[ 0 ];
							coord[ 1 ] = min_point[ 1 ] + cell_size;
							coord[ 2 ] = min_point[ 2 ];
						}
						if(  j == 1  &&  k == 0  &&  h == 0  ){
							node = 4;
							coord[ 0 ] = min_point[ 0 ];
							coord[ 1 ] = min_point[ 1 ];
							coord[ 2 ] = min_point[ 2 ] + cell_size;
						}
						if(  j == 1  &&  k == 0  &&  h == 1  ){
							node = 5;
							coord[ 0 ] = min_point[ 0 ] + cell_size;
							coord[ 1 ] = min_point[ 1 ];
							coord[ 2 ] = min_point[ 2 ] + cell_size;
						}
						if(  j == 1  &&  k == 1  &&  h == 0  ){
							node = 6;
							coord[ 0 ] = min_point[ 0 ] + cell_size;
							coord[ 1 ] = min_point[ 1 ] + cell_size;
							coord[ 2 ] = min_point[ 2 ] + cell_size;
						}
						if(  j == 1  &&  k == 1  &&  h == 1  ){
							node = 7;
							coord[ 0 ] = min_point[ 0 ];
							coord[ 1 ] = min_point[ 1 ] + cell_size;
							coord[ 2 ] = min_point[ 2 ] + cell_size;
						}
						
						if(  !data->_NodeWasPrinted(  node  )  ){
							//Rotating coordinate
							double new_coord[ 3 ];
							new_coord[ 0 ] = vdb_->UnscaleCoordinate( 0 , coord[ 0 ] );
							new_coord[ 1 ] = vdb_->UnscaleCoordinate( 1 , coord[ 1 ] );
							new_coord[ 2 ] = vdb_->UnscaleCoordinate( 2 , coord[ 2 ] );
							this->RotateCoordianteToFold( new_coord );
							fp << new_coord[ 0 ] << " " << new_coord[ 1 ] << " " << new_coord[ 2 ] << endl;
							this->SetPrintedOnNodeAndNeighbours( node , cell );
						}
					}
				}
			}
		}
	}
	fp << endl;

	fp << "{Mesh}" << endl;
	fp << "5 ; Element type (2=Triangle, 3=Quadrilateral, 4=Tetrahedra, 5=Hexahedra)" << endl;
	fp << "8 ; Nodes per element" << endl;
	fp << n_elements << " ; Elements count" << endl;
	fp << "; Material Node1 Node2 ..." << endl;
	for(  size_t i = 0  ;  i < leaves.size()  ;  i++  ){
		OctreeCell* cell = leaves[ i ];
		DataInsideOctreeCell* data = cell->pGetData(  );
		if(  data->_GetIsMesh() == 1  ){
			if( data->_IsLoaded() ){
				fp << "2 " << data->_GetIndex( 0 ) + 1 << " " << data->_GetIndex( 1 ) + 1 << " " 
					 << data->_GetIndex( 2 ) + 1 << " " << data->_GetIndex( 3 ) + 1 << " " 
					 << data->_GetIndex( 4 ) + 1 << " " << data->_GetIndex( 5 ) + 1 << " " 
					 << data->_GetIndex( 6 ) + 1 << " " << data->_GetIndex( 7 ) + 1 << endl;
			}else{
				fp << "1 " << data->_GetIndex( 0 ) + 1 << " " << data->_GetIndex( 1 ) + 1 << " " 
					 << data->_GetIndex( 2 ) + 1 << " " << data->_GetIndex( 3 ) + 1 << " " 
					 << data->_GetIndex( 4 ) + 1 << " " << data->_GetIndex( 5 ) + 1 << " " 
					 << data->_GetIndex( 6 ) + 1 << " " << data->_GetIndex( 7 ) + 1 << endl;
			}
		}
	}

	fp.close();
}




