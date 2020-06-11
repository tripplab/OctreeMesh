#include "interpolator.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                   INTERPOLATOR METHODS                                //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Interpolator::Interpolator(){
	strcpy( input_ , "\0" );
	vdb_ = NULL;
	strcpy( results_ , "\0" );
	strcpy( output_ , "\0" );
	r_level_ = 0 ;
	octree_ = NULL;
	local_root_ = NULL;
}

/**
 *Constructor receiving program arguments 
 *This constructor initializes the parallel region in each process, also reads the vdb 
 *file containing the capside information
 *@param[in] argc Is the number inidicating the program arguments
 *@param[in] argv Is the list of arguments for the program
 */
Interpolator::Interpolator( int argc , char** argv ){
	local_root_ = NULL;

	//Parameter one VdB mesh
	input_ = argv[ 1 ];
	vdb_ = new Vdb( input_ );
	this->ReadVdbFile(  );
  this->ScaleCapside(  );

	//Parameter 2 hexahedral mesh
	mesh_ = argv[ 2 ];
	hex_ = new Hexahedral( mesh_ );
  this->ScaleHexahedralMesh(  );

	//Parameter 3 results file
	results_ = argv[ 3 ];
	res_ = new Results( results_ , hex_->GetNNodes() , hex_->GetNElements() );
	
	//Parameter 4 output name
	output_ = argv[ 4 ];

	//parameter 5 refinement level
	r_level_ = atoi( argv[ 5 ] );
	octree_ = new OctreeDriver();
	local_root_ = octree_->GetRoot();
}

/**
 *Default destructor 
 *Finalizes the parallel region
 */
Interpolator::~Interpolator(  ){
}

//GETS

/**
 *Getting all leaves on local root
 *@param[out] leaves Is the vector that will contain all leaves on local root
 */
void Interpolator::GetAllLocalRootLeaves( OctreeCell_vector& leaves ){
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
OctreeCell* Interpolator::GetCellOnLocalOctree( key_type* keys ){

	OctreeCell* cell = local_root_;
	for(  size_t i = 0  ;  i < ROOT_LEVEL  ;  i++  ){
		if(  cell->IsLeaf()  ){
			return cell;
		}
		cell = cell->pGetChild( keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
	}
	return cell;
}

//UTILITIES
/**
 *Method to read the information of the vdb file previusly opened
 */
void Interpolator::ReadVdbFile(  ){
	vdb_->ReadCompleteFile(  );
}


/**
 *Scaling boundary nodes
 *This method calls to the Boundary::ScaleNodes method from Boundary class
 */
void Interpolator::ScaleCapside(  ){
	vdb_->ScaleSpheresMesh();
	vdb_->PrintResumenOnWarnings();
}

/**
 *Refining local cells based on the refinement level given as parameter
 *This method subdivide the leaves of the local root until reach the r_level_
 */
void Interpolator::RefineLocalRoot(){

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

					//ASSIGNING HEXAHEDRONS ON CHILDS
					size_t n_hexahedrons = data->_GetNHexahedrons( );
					for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
						OctreeCell* child = leaf->pGetChild( i_child );
						double child_min_coord[ 3 ],child_max_coord[ 3 ], child_center[ 3 ], child_width[ 3 ];
						child->GetMinPointNormalized( child_min_coord );
						child->GetMaxPointNormalized( child_max_coord );
						for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
							child_center[ i_pos ] = ( child_min_coord[ i_pos ] + child_max_coord[ i_pos ] ) * 0.5;
							child_width[ i_pos ] = ( child_max_coord[ i_pos ] - child_min_coord[ i_pos ] ) * 0.5;
						}
						for(  size_t i_hex = 0  ;  i_hex < n_hexahedrons  ;  i_hex++  ){
							Hexahedron* aux = data->_GetHexahedron( i_hex );
							double min_coord[ 3 ], max_coord[ 3 ], center[ 3 ], width[ 3 ];
							aux->GetMinPoint( min_coord );
							aux->GetMaxPoint( max_coord );
							for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
								center[ i_dim ] = 0.5 * ( min_coord[ i_dim ] + max_coord[ i_dim ] );
								width[ i_dim ] = 0.5 * ( max_coord[ i_dim ]- min_coord[ i_dim ] );
							}
							if(  AABBvsAABBDifferentSize( child_center , child_width , center , width )  ){
								data_children[ i_child ]->_AddHexahedron( aux );
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
void Interpolator::AssignAtomsOnLocalRoot( ){

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
 *Scaling the hexahedral mesh to the domain[0,1]
 */
void Interpolator::ScaleHexahedralMesh(  ){
	hex_->Rotate( );
	hex_->ScaleMesh();
}

/**
 *Translating strain, stress, and vonmises to the nodes in the mesh and savin this results
 *file with the same name that the post.res file but adding the prefix nodes
 *NOTE: The set of values on the gauss points are translated using the next information
 *
 *	GAUSS INDEX    			TRASLATED TO				ELEMENT INDEX
 *			 0  						------------>							 0

 *			 1  						------------>							 4
 *			 2  						------------>							 3
 *			 3  						------------>							 7
 *			 4  						------------>							 1
 *			 5  						------------>							 5
 *			 6  						------------>							 2
 *			 7  						------------>							 6
 *In this case is only required the displacement and the vonmises
 */
void Interpolator::TranslateAllResultsToNodes(  ){

	size_t n_nodes = res_->GetNNodes(  );
	size_t n_elements = res_->GetNElements(  );
	int*	counter; //Counter of number of values added on the mesh
	res_nodes_ = new double*[ 14 ];
	counter = new int[ n_nodes ];
	for(  int i_value = 0  ;  i_value < 14  ;  i_value++  ){
		res_nodes_[ i_value ] =  new double[ n_nodes ];
		for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){//Fill with zeros
			res_nodes_[ i_value ][ i_node ] = 0.0;
			counter[ i_node ] = 0;
		}
	}
	//Setting displacements on row 0
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		res_nodes_[ 0 ][ i_node ] = res_->GetDisplacementMagnitude( i_node );
	}
	//Seting strain,stress and vonmises
	int local_real[ 8 ] = { 0,4,3,7,1,5,2,6 };
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		for(  int i_gauss = 0  ;  i_gauss < 8  ;  i_gauss++  ){
			for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
					res_nodes_[ i_kind+1 ][ hex_->GetIndex( i_elem , local_real[ i_gauss ] ) ] += res_->GetStrain( i_elem , i_kind , i_gauss );
					res_nodes_[ i_kind+7 ][ hex_->GetIndex( i_elem , local_real[ i_gauss ] ) ] += res_->GetStress( i_elem , i_kind , i_gauss );
			}
			res_nodes_[ 13 ][ hex_->GetIndex( i_elem , local_real[ i_gauss ] ) ] += res_->GetVonmises( i_elem , i_gauss );
			counter[ hex_->GetIndex( i_elem , local_real[ i_gauss ] ) ]++;
		}
	}
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		for(  int i_row = 1  ;  i_row < 14  ;  i_row++  ){
			res_nodes_[ i_row ][ i_node ] = res_nodes_[ i_row ][ i_node ] / (double)counter[ i_node ];
		}
	}
	//SCALING TO [0-100] THE DISPLACEMENT
	double min = 1e20;
	double max = -1e20;
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if( res_nodes_[ 0 ][ i_node ] < min ){
			min = res_nodes_[ 0 ][ i_node ];
		}
		if( res_nodes_[ 0 ][ i_node ] > max ){
			max = res_nodes_[ 0 ][ i_node ];
		}
	}
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		res_nodes_[ 0 ][ i_node ] = (( res_nodes_[ 0 ][ i_node ] - min )/( max - min )*99.0)+1.0;
		res_nodes_[ 0 ][ i_node ] = log( res_nodes_[ 0 ][ i_node ] );
	}
	//SCALING TO [0-100] THE VONMISES VALUE
	min = 1e20;
	max = -1e20;
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if( res_nodes_[ 13 ][ i_node ] < min ){
			min = res_nodes_[ 13 ][ i_node ];
		}
		if( res_nodes_[ 13 ][ i_node ] > max ){
			max = res_nodes_[ 13 ][ i_node ];
		}
	}
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		res_nodes_[ 13 ][ i_node ] = (( res_nodes_[ 13 ][ i_node ] - min )/( max - min )*99.00)+1.0;
		res_nodes_[ 13 ][ i_node ] = log( res_nodes_[ 13 ][ i_node ] );
	}
	delete counter;
}

/**
 *Assigning the hexahedral eleemnts on the local roots, in this
 */
void Interpolator::AssignHexahedralsOnLocalRoot( ){
	//Getting local root information
	DataInsideOctreeCell* data = local_root_->pGetData();

	//Assigning information of hexahedral mesh on local_root
	size_t n_elements = hex_->GetNElements();
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		Hexahedron* aux = hex_->GetElement( i_elem );
		data->_AddHexahedron( aux );
	}
}

/**
 *Interpolates the results (displacement magnitude and vonmises) to the center of the atom
 *This method takes the atoms on the vdb file and for the valid atoms looks for the element
 *that contains its center and looks the psi,eta and mu coordinate to interpolate the result
 *of the displacement magnitud and vonmises
 */
void Interpolator::InterpolateResultsToAtomsAndSavePdb(){
	FILE* fin = NULL;
	FILE* fout = NULL;
	
	fin = fopen( input_ , "r" );
	fout = fopen( output_ , "w" );
	//PRinting header on output file
	fprintf( fout , "REMARK 220\n"  );	
	fprintf( fout , "REMARK 220 EXPERIMENTAL DETAILS\n");
	fprintf( fout , "REMARK 220  EXPERIMENT TYPE                : THEORETICAL MODELLING\n");
	fprintf( fout , "REMARK 225\n");
	fprintf( fout , "REMARK 225 CapsidMesh to PDB Model\n");
	fprintf( fout , "REMARK 225\n");
	fprintf( fout , "REMARK 225 Nanoindentation numerical simulation performed with FEM\n");
	fprintf( fout , "REMARK 225 Authors:\n"); 
	fprintf( fout , "REMARK 225 Parameters\n");
	fprintf( fout , "REMARK 225 Mesh resolution 5.0 A, Force 0.5 nN, Poisson 0.3\n");
	fprintf( fout , "REMARK 225 Elasticity 200 MPa, applied on 5-fold\n");
	fprintf( fout , "REMARK 225 Magnitude of indentation is located in the 'occupancy' field\n"); 
	fprintf( fout , "REMARK 225 Magnitude of von Mises is located in the 'tempFactor' field\n");


	while( !feof( fin ) ){
		char params[ 11 ][ 500 ];
		char line[500];
		fgets( line , 500 , fin );
		//if is end of file does not do anyting
		if(  feof( fin ) ){
			break;
		}
		//reading column 1
		int position = 0;
		for(  int i_char = 0  ;  i_char < 4  ;  i_char++  ){
			params[ 0 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 0 ][ 4 ] = 0;
	
		//reading column 2
		position = 6;
		for(  int i_char = 0  ;  i_char < 5  ;  i_char++  ){
			params[ 1 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 1 ][ 5 ] = 0;

		//reading column 3
		position = 12;
		for(  int i_char = 0  ;  i_char < 4  ;  i_char++  ){
			params[ 2 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 2 ][ 4 ] = 0;
	 
		//reading column 4
		position = 17;
		for(  int i_char = 0  ;  i_char < 3  ;  i_char++  ){
			params[ 3 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 3 ][ 3 ] = 0;

		//reading column 5
		position = 21;
		for(  int i_char = 0  ;  i_char < 1  ;  i_char++  ){
			params[ 4 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 4 ][ 1 ] = 0;

		//reading column 6
		position = 22;
		for(  int i_char = 0  ;  i_char < 4  ;  i_char++  ){
			params[ 5 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 5 ][ 4 ] = 0;

		//reading column 7
		position = 30;
		for(  int i_char = 0  ;  i_char < 8  ;  i_char++  ){
			params[ 6 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 6 ][ 8 ] = 0;

		//reading column 8
		position = 38;
		for(  int i_char = 0  ;  i_char < 8  ;  i_char++  ){
			params[ 7 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 7 ][ 8 ] = 0;

		//reading column 9
		position = 46;
		for(  int i_char = 0  ;  i_char < 8  ;  i_char++  ){
			params[ 8 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 8 ][ 8 ] = 0;

		//reading column 10
		position = 54;
		for(  int i_char = 0  ;  i_char < 6  ;  i_char++  ){
			params[ 9 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 9 ][ 6 ] = 0;

		//reading column 11
		position = 60;
		for(  int i_char = 0  ;  i_char < 6  ;  i_char++  ){
			params[ 10 ][ i_char ] = line[ position ];
			position++;
		}
		params[ 10 ][ 6 ] = 0;

		if(  !strcmp( params[ 0 ] , "ATOM" )  ){
			//In this case the information correcponds to an ATOM, but some of the aminoacids are 
		  //not valid so neet to be tested if the aminoacid is valid or not
			int index;
			//testing if is a valid aminoacid or not
			if(  vdb_->AminoacidIsValid( params[ 3 ] , &index )  ){
				//Getting coordinates
				double coords[ 3 ];
				sscanf( params[ 6 ] , "%lf" , &coords[ 0 ] );
				sscanf( params[ 7 ] , "%lf" , &coords[ 1 ] );
				sscanf( params[ 8 ] , "%lf" , &coords[ 2 ] );
				for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
					coords[ i_pos ] = vdb_->ScaleCoord( i_pos , coords[ i_pos ] );
				}
				//Getting Key from coordinates
				size_t key[ 3 ];
				for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
					key[ i_pos ] = static_cast<key_type> ( (1 << 29) * coords[ i_pos ] );
				}
				//Looking for the cell containing the atom
				OctreeCell* cell = local_root_;
				for(  size_t i_level = 0  ;  i_level < ROOT_LEVEL  ;  i_level++  ){
					if(  cell->IsLeaf()  ){
						break;
					}
					cell = cell->pGetChild( key[ 0 ] , key[ 1 ] , key[ 2 ] );
				}
				//looking the hexahedron containing the atom center
				DataInsideOctreeCell* data = cell->pGetData();
				size_t n_hexa = data->_GetNHexahedrons();
				if(  n_hexa == 0  ){
					bool flag_neigh = false;
					key_type neighbour_key[ 3 ];
					for(  int i_direction = 0  ;  i_direction < 18  ;  i_direction++  ){
						if(  cell->GetNeighbourKey( i_direction , neighbour_key )  ){
								OctreeCell* neighbour_cell = this->GetCellOnLocalOctree( neighbour_key );
								DataInsideOctreeCell* neighbour_data = neighbour_cell->pGetData();
								n_hexa = neighbour_data->_GetNHexahedrons();
								if(  n_hexa > 0  ){
									flag_neigh = true;								
									size_t close_index;
									double close_dist = 1e20;
									for(  size_t i_hexa = 0  ;  i_hexa < n_hexa  ;  i_hexa++  ){
										Hexahedron* aux = neighbour_data->_GetHexahedron( i_hexa );
										double dist = aux->CalcDistance( coords );
										if(  dist < close_dist  ){
											close_dist = dist;
											close_index = aux->GetIndex( 0 );
										}
									}
									double displacement = res_nodes_[ 0 ][ close_index ];
									double vonmises = res_nodes_[ 13 ][ close_index ];
									//**********PRINT RESULTS ON FILE*******************************************
									double residuo = 0.0;
									if(  displacement > 999.99  ){
										displacement = 999.99;
									}
									int centenas = displacement/100;
									residuo = displacement-centenas*100;
									int decenas = residuo/10;
									residuo = residuo-decenas*10;
									int unidades = residuo/1;
									residuo = residuo-unidades;
									int decima = residuo/0.1;
									residuo = residuo-0.1*decima;
									int centesima = residuo/0.01;
									residuo = residuo-0.01*centesima;



									char aux_char[2];
									if(  centenas == 0  ){
										line[ 54 ] = ' ';
										if(  decenas == 0  ){
											line[ 55 ] = ' ';										
											if(  unidades == 0  ){
												line[ 56 ] = '0';
												line[ 57 ] = '.';
												sprintf( aux_char , "%d" , decima );
												line[ 58 ] = aux_char[0];
												sprintf( aux_char , "%d" , centesima );
												line[ 59 ] = aux_char[0];											
											}else{
												sprintf( aux_char , "%d" , unidades );
												line[ 56 ] = aux_char[0];
												line[ 57 ] = '.';
												sprintf( aux_char , "%d" , decima );
												line[ 58 ] = aux_char[0];
												sprintf( aux_char , "%d" , centesima );
												line[ 59 ] = aux_char[0];											
											}
										}else{
											sprintf( aux_char , "%d" , decenas );
											line[ 55 ] = aux_char[0];
											sprintf( aux_char , "%d" , unidades );
											line[ 56 ] = aux_char[0];
											line[ 57 ] = '.';
											sprintf( aux_char , "%d" , decima );
											line[ 58 ] = aux_char[0];
											sprintf( aux_char , "%d" , centesima );
											line[ 59 ] = aux_char[0];
										}
									}else{
										sprintf( aux_char , "%d" , centenas );
										line[ 54 ] = aux_char[0];
										sprintf( aux_char , "%d" , decenas );
										line[ 55 ] = aux_char[0];
										sprintf( aux_char , "%d" , unidades );
										line[ 56 ] = aux_char[0];
										line[ 57 ] = '.';
										sprintf( aux_char , "%d" , decima );
										line[ 58 ] = aux_char[0];
										sprintf( aux_char , "%d" , centesima );
										line[ 59 ] = aux_char[0];									
									}
									if(  vonmises > 999.99  ){
										vonmises = 999.99;
									}
									centenas = vonmises/100;
									residuo = vonmises-centenas*100;
									decenas = residuo/10;
									residuo = residuo-decenas*10;
									unidades = residuo/1;
									residuo = residuo-unidades;
									decima = residuo/0.1;
									residuo = residuo-0.1*decima;
									centesima = residuo/0.01;
									residuo = residuo-0.01*centesima;

									if(  centenas == 0  ){
										line[ 60 ] = ' ';
										if(  decenas == 0  ){
											line[ 61 ] = ' ';										
											if(  unidades == 0  ){
												line[ 62 ] = '0';
												line[ 63 ] = '.';
												sprintf( aux_char , "%d" , decima );
												line[ 64 ] = aux_char[0];
												sprintf( aux_char , "%d" , centesima );
												line[ 65 ] = aux_char[0];											
											}else{
												sprintf( aux_char , "%d" , unidades );
												line[ 62 ] = aux_char[0];
												line[ 63 ] = '.';
												sprintf( aux_char , "%d" , decima );
												line[ 64 ] = aux_char[0];
												sprintf( aux_char , "%d" , centesima );
												line[ 65 ] = aux_char[0];											
											}
										}else{
											sprintf( aux_char , "%d" , decenas );
											line[ 61 ] = aux_char[0];
											sprintf( aux_char , "%d" , unidades );
											line[ 62 ] = aux_char[0];
											line[ 63 ] = '.';
											sprintf( aux_char , "%d" , decima );
											line[ 64 ] = aux_char[0];
											sprintf( aux_char , "%d" , centesima );
											line[ 65 ] = aux_char[0];
										}
									}else{
										sprintf( aux_char , "%d" , centenas );
										line[ 60 ] = aux_char[0];
										sprintf( aux_char , "%d" , decenas );
										line[ 61 ] = aux_char[0];
										sprintf( aux_char , "%d" , unidades );
										line[ 62 ] = aux_char[0];
										line[ 63 ] = '.';
										sprintf( aux_char , "%d" , decima );
										line[ 64 ] = aux_char[0];
										sprintf( aux_char , "%d" , centesima );
										line[ 65 ] = aux_char[0];									
									}
									fprintf( fout , "%s" , line );
									//**********PRINT RESULTS ON FILE*******************************************
									break;
								}
						}
					}
					if(  !flag_neigh  ){
						//**********PRINT RESULTS ON FILE*******************************************
						line[ 54 ] = ' ';
						line[ 55 ] = ' ';
						line[ 56 ] = '0';
						line[ 57 ] = '.';
						line[ 58 ] = '0';
						line[ 59 ] = '0';
						line[ 60 ] = ' ';
						line[ 61 ] = ' ';
						line[ 62 ] = '0';
						line[ 63 ] = '.';
						line[ 64 ] = '0';
						line[ 65 ] = '0';
						fprintf( fout , "%s" , line );
						//**********PRINT RESULTS ON FILE*******************************************		
					}
					continue;
				}
				bool flag = false;
				for(  size_t i_hexa = 0  ;  i_hexa < n_hexa  ;  i_hexa++  ){
					Hexahedron* aux = data->_GetHexahedron( i_hexa );
					if(  aux->ContainPoint( coords )  ){
						flag = true;
						double psi,eta,mu,weights[ 8 ];
						double displacement = 0.0;
						double vonmises = 0.0;
						size_t n_index[ 8 ];
						aux->GetNormalizedCoordinates( coords[ 0 ] , coords[ 1 ] , coords[ 2 ] , &psi , &eta , &mu );	
						for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
							weights[ i_node ] = aux->FormFunction( i_node , psi , eta , mu );
							n_index[ i_node ] = aux->GetIndex( i_node );
							displacement += (res_nodes_[ 0 ][ n_index[ i_node ] ]*weights[ i_node ]);
							vonmises     += res_nodes_[ 13 ][ n_index[ i_node ] ]*weights[ i_node ];
						}
						//**********PRINT RESULTS ON FILE*******************************************
						double residuo = 0.0;
						if(  displacement > 999.99  ){
							displacement = 999.99;
						}
						int centenas = displacement/100;
						residuo = displacement-centenas*100;
						int decenas = residuo/10;
						residuo = residuo-decenas*10;
						int unidades = residuo/1;
						residuo = residuo-unidades;
						int decima = residuo/0.1;
						residuo = residuo-0.1*decima;
						int centesima = residuo/0.01;
						residuo = residuo-0.01*centesima;



						char aux_char[2];
						if(  centenas == 0  ){
							line[ 54 ] = ' ';
							if(  decenas == 0  ){
								line[ 55 ] = ' ';										
								if(  unidades == 0  ){
									line[ 56 ] = '0';
									line[ 57 ] = '.';
									sprintf( aux_char , "%d" , decima );
									line[ 58 ] = aux_char[0];
									sprintf( aux_char , "%d" , centesima );
									line[ 59 ] = aux_char[0];											
								}else{
									sprintf( aux_char , "%d" , unidades );
									line[ 56 ] = aux_char[0];
									line[ 57 ] = '.';
									sprintf( aux_char , "%d" , decima );
									line[ 58 ] = aux_char[0];
									sprintf( aux_char , "%d" , centesima );
									line[ 59 ] = aux_char[0];											
								}
							}else{
								sprintf( aux_char , "%d" , decenas );
								line[ 55 ] = aux_char[0];
								sprintf( aux_char , "%d" , unidades );
								line[ 56 ] = aux_char[0];
								line[ 57 ] = '.';
								sprintf( aux_char , "%d" , decima );
								line[ 58 ] = aux_char[0];
								sprintf( aux_char , "%d" , centesima );
								line[ 59 ] = aux_char[0];
							}
						}else{
							sprintf( aux_char , "%d" , centenas );
							line[ 54 ] = aux_char[0];
							sprintf( aux_char , "%d" , decenas );
							line[ 55 ] = aux_char[0];
							sprintf( aux_char , "%d" , unidades );
							line[ 56 ] = aux_char[0];
							line[ 57 ] = '.';
							sprintf( aux_char , "%d" , decima );
							line[ 58 ] = aux_char[0];
							sprintf( aux_char , "%d" , centesima );
							line[ 59 ] = aux_char[0];									
						}

						if(  vonmises > 999.99  ){
							vonmises = 999.99;
						}
						centenas = vonmises/100;
						residuo = vonmises-centenas*100;
						decenas = residuo/10;
						residuo = residuo-decenas*10;
						unidades = residuo/1;
						residuo = residuo-unidades;
						decima = residuo/0.1;
						residuo = residuo-0.1*decima;
						centesima = residuo/0.01;
						residuo = residuo-0.01*centesima;

						if(  centenas == 0  ){
							line[ 60 ] = ' ';
							if(  decenas == 0  ){
								line[ 61 ] = ' ';										
								if(  unidades == 0  ){
									line[ 62 ] = '0';
									line[ 63 ] = '.';
									sprintf( aux_char , "%d" , decima );
									line[ 64 ] = aux_char[0];
									sprintf( aux_char , "%d" , centesima );
									line[ 65 ] = aux_char[0];											
								}else{
									sprintf( aux_char , "%d" , unidades );
									line[ 62 ] = aux_char[0];
									line[ 63 ] = '.';
									sprintf( aux_char , "%d" , decima );
									line[ 64 ] = aux_char[0];
									sprintf( aux_char , "%d" , centesima );
									line[ 65 ] = aux_char[0];											
								}
							}else{
								sprintf( aux_char , "%d" , decenas );
								line[ 61 ] = aux_char[0];
								sprintf( aux_char , "%d" , unidades );
								line[ 62 ] = aux_char[0];
								line[ 63 ] = '.';
								sprintf( aux_char , "%d" , decima );
								line[ 64 ] = aux_char[0];
								sprintf( aux_char , "%d" , centesima );
								line[ 65 ] = aux_char[0];
							}
						}else{
							sprintf( aux_char , "%d" , centenas );
							line[ 60 ] = aux_char[0];
							sprintf( aux_char , "%d" , decenas );
							line[ 61 ] = aux_char[0];
							sprintf( aux_char , "%d" , unidades );
							line[ 62 ] = aux_char[0];
							line[ 63 ] = '.';
							sprintf( aux_char , "%d" , decima );
							line[ 64 ] = aux_char[0];
							sprintf( aux_char , "%d" , centesima );
							line[ 65 ] = aux_char[0];									
						}
						fprintf( fout , "%s" , line );
						//**********PRINT RESULTS ON FILE*******************************************
						break;
					}
				}
				if(  !flag  ){//Center is not inside any hexahedron, looking for the closest 
											//hexahedron and using values on node 0
					size_t close_index;
					double close_dist = 1e20;
					for(  size_t i_hexa = 0  ;  i_hexa < n_hexa  ;  i_hexa++  ){
						Hexahedron* aux = data->_GetHexahedron( i_hexa );
						double dist = aux->CalcDistance( coords );
						if(  dist < close_dist  ){
							close_dist = dist;
							close_index = aux->GetIndex( 0 );
						}
					}
					double displacement = res_nodes_[ 0 ][ close_index ];
					double vonmises = res_nodes_[ 13 ][ close_index ];
					//**********PRINT RESULTS ON FILE*******************************************
					double residuo = 0.0;
					if(  displacement > 999.99  ){
						displacement = 999.99;
					}
					int centenas = displacement/100;
					residuo = displacement-centenas*100;
					int decenas = residuo/10;
					residuo = residuo-decenas*10;
					int unidades = residuo/1;
					residuo = residuo-unidades;
					int decima = residuo/0.1;
					residuo = residuo-0.1*decima;
					int centesima = residuo/0.01;
					residuo = residuo-0.01*centesima;



					char aux_char[2];
					if(  centenas == 0  ){
						line[ 54 ] = ' ';
						if(  decenas == 0  ){
							line[ 55 ] = ' ';										
							if(  unidades == 0  ){
								line[ 56 ] = '0';
								line[ 57 ] = '.';
								sprintf( aux_char , "%d" , decima );
								line[ 58 ] = aux_char[0];
								sprintf( aux_char , "%d" , centesima );
								line[ 59 ] = aux_char[0];											
							}else{
								sprintf( aux_char , "%d" , unidades );
								line[ 56 ] = aux_char[0];
								line[ 57 ] = '.';
								sprintf( aux_char , "%d" , decima );
								line[ 58 ] = aux_char[0];
								sprintf( aux_char , "%d" , centesima );
								line[ 59 ] = aux_char[0];											
							}
						}else{
							sprintf( aux_char , "%d" , decenas );
							line[ 55 ] = aux_char[0];
							sprintf( aux_char , "%d" , unidades );
							line[ 56 ] = aux_char[0];
							line[ 57 ] = '.';
							sprintf( aux_char , "%d" , decima );
							line[ 58 ] = aux_char[0];
							sprintf( aux_char , "%d" , centesima );
							line[ 59 ] = aux_char[0];
						}
					}else{
						sprintf( aux_char , "%d" , centenas );
						line[ 54 ] = aux_char[0];
						sprintf( aux_char , "%d" , decenas );
						line[ 55 ] = aux_char[0];
						sprintf( aux_char , "%d" , unidades );
						line[ 56 ] = aux_char[0];
						line[ 57 ] = '.';
						sprintf( aux_char , "%d" , decima );
						line[ 58 ] = aux_char[0];
						sprintf( aux_char , "%d" , centesima );
						line[ 59 ] = aux_char[0];									
					}

					if(  vonmises > 999.99  ){
						vonmises = 999.99;
					}
					centenas = vonmises/100;
					residuo = vonmises-centenas*100;
					decenas = residuo/10;
					residuo = residuo-decenas*10;
					unidades = residuo/1;
					residuo = residuo-unidades;
					decima = residuo/0.1;
					residuo = residuo-0.1*decima;
					centesima = residuo/0.01;
					residuo = residuo-0.01*centesima;

					if(  centenas == 0  ){
						line[ 60 ] = ' ';
						if(  decenas == 0  ){
							line[ 61 ] = ' ';										
							if(  unidades == 0  ){
								line[ 62 ] = '0';
								line[ 63 ] = '.';
								sprintf( aux_char , "%d" , decima );
								line[ 64 ] = aux_char[0];
								sprintf( aux_char , "%d" , centesima );
								line[ 65 ] = aux_char[0];											
							}else{
								sprintf( aux_char , "%d" , unidades );
								line[ 62 ] = aux_char[0];
								line[ 63 ] = '.';
								sprintf( aux_char , "%d" , decima );
								line[ 64 ] = aux_char[0];
								sprintf( aux_char , "%d" , centesima );
								line[ 65 ] = aux_char[0];											
							}
						}else{
							sprintf( aux_char , "%d" , decenas );
							line[ 61 ] = aux_char[0];
							sprintf( aux_char , "%d" , unidades );
							line[ 62 ] = aux_char[0];
							line[ 63 ] = '.';
							sprintf( aux_char , "%d" , decima );
							line[ 64 ] = aux_char[0];
							sprintf( aux_char , "%d" , centesima );
							line[ 65 ] = aux_char[0];
						}
					}else{
						sprintf( aux_char , "%d" , centenas );
						line[ 60 ] = aux_char[0];
						sprintf( aux_char , "%d" , decenas );
						line[ 61 ] = aux_char[0];
						sprintf( aux_char , "%d" , unidades );
						line[ 62 ] = aux_char[0];
						line[ 63 ] = '.';
						sprintf( aux_char , "%d" , decima );
						line[ 64 ] = aux_char[0];
						sprintf( aux_char , "%d" , centesima );
						line[ 65 ] = aux_char[0];									
					}
					fprintf( fout , "%s" , line );
					//**********PRINT RESULTS ON FILE*******************************************
				}
			}else{//aminoacid not valid, printing line in the output file as was readed
				fprintf( fout , "%s" , line );
			}
		}else{//this line is auxiliary information, is printed equal in the output file
			fprintf( fout , "%s" , line );
		}
	}


	fclose( fin );
	fclose( fout );
}


















