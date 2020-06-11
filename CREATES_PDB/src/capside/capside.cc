#include "capside.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                   ATOM   METHODS                                     //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Atom::Atom(  ){
	strcpy( kind_ , "\0" );
	radius_ = 0.0;
	coords_[ 0 ] = 0.0;
	coords_[ 1 ] = 0.0;
	coords_[ 2 ] = 0.0;
	keys_[ 0 ] = 0;
	keys_[ 1 ] = 0;
	keys_[ 2 ] = 0;
}

/**
 *Constructor receibing the name of the atom and its dimaeter
 *@param[in] kind Is the name of the atom, C, NA, etc.
 *@param[in] radius Is the radius of the atom
 *@param[in] coords Is the spacial position of the atom
 */
Atom::Atom( char* kind , double radius , double* coords ){
	strcpy( kind_ , kind );
	radius_ = radius;
	for(  int i_position = 0  ;  i_position < 3  ;  i_position++  ){
		coords_[ i_position ] = coords[ i_position ]; 
		keys_[ i_position ] = static_cast<size_t> ( ( 1 << 29 ) * coords_[ i_position ]  );
	}
}

/**
 *Default destructor
 */
Atom::~Atom(  ){
	strcpy( kind_ , "\0" );
  radius_ = 0.0;
	coords_[ 0 ] = 0.0;
	coords_[ 1 ] = 0.0;
	coords_[ 2 ] = 0.0;
}

//GETS
/**
 *Getting the kind of atom (its name)
 *@param[out] kind Is the string containing the name of the atom
 *@return A bool value indicating if the atom have the name set 
 */
bool Atom::GetKind( char kind[10] ){
	if(  strcmp( kind_ , "\0" )  ){
		strcpy( kind , kind_ );
		return true;
	}else{
		return false;
	}
	assert( 0 );//this never have to happen
}

/**
 *Getting coordinate from atom
 *@param[in] position Is the index of the axis where the coordinate is required
 *                    1-X axis
 *                    2-Y axis
 *                    3-Z axis
 *return A double value with the coordinate requested
 */
double Atom::GetCoordinate( int position ){
	return coords_[ position ];
}

/**
 *Getting coordinates from atom
 *@param[out] coords Is the array where the coordinates will be stored
 */
void Atom::GetCoordinates( double* coords ){
	for(  int i_position = 0  ;  i_position < 3  ;  i_position++  ){
		coords[ i_position ] =  coords_[ i_position ];
	}
}

/**
 *Getting atom radius
 *@return A double value with the atom radius
 */
double Atom::GetRadius(){
	return radius_;
}

/**
 *Getting the material index for the atom based on:
 *   1-Carbon, 2-Hydrogen , 3-Oxygen , 4-Nitrogen , 5-phosphorus , 6-Sulfur
 *@retutrn a int value witht he material
 */
int Atom::GetMaterialIndex(  ){
	char temp = kind_[ 0 ];
	if(  temp == ' '  ){
		temp = kind_[ 1 ];
		if(  temp == ' '  ){
			temp = kind_[ 2 ];
			if(  temp == ' '  ){
				temp = kind_[ 3 ];
			}
		}
	}
	int index;
	switch( temp ){
		case 'C':
			index = 1;
			break;
		case 'H':
			index = 2;
			break;
		case 'O':
			index = 3;
			break;
		case 'N':
			index = 4;
			break;
		case 'P':
			index = 5;
			break;
		case 'S':
			index = 6;
			break;
		default:
			std::cout << "Error: Atom not valid " << std::endl;
			assert( 0 );
			break;
	}
	return index;
}

/**
 *Getting key from atom
 *@param[in] position Is the dimension of the key requested
 *@return A size_t value with the key requested
 */
size_t Atom::GetKey( int position ){
	return keys_[ position ];
}

//SETS
/**
 *Setting kind of atom
 *@param[in] kind is the name of the atom that will be set
 *@return A bool value indicating if the name was set correctly
 */
bool Atom::SetKind( char* kind ){
	strcpy( kind_ , kind );
	bool flag = false;
	if(  !strcmp( kind_ , kind )  ){
		flag = true;
	}
	return flag;
}		

/**
 *Setting coordinate on atom
 *@param[in] coord  Is the coordinate value to be set
 *@param[in] position Is the index of the axis where the coordinate is required
 *                    1-X axis
 *                    2-Y axis
 *                    3-Z axis
 *@return A bool value insicating if the coordinate was set correctly
 */
bool Atom::SetCoordinate( double coord , int position ){
	coords_[ position ] = coord;
	return ( coords_[ position ] == coord );
}

/**
 *Setting coordinates on atom
 *@param[in] coords Is the array with the atom coordinates
 *@return A bool value indicating if the coordinates where set correctly
 */
bool Atom::SetCoordinates( double* coords ){
	for(  int i_position = 0  ;  i_position < 3  ;  i_position++  ){
		coords_[ i_position ] = coords[ i_position ];
		if(  coords_[ i_position ] != coords[ i_position ]  ){
			return false;
		}
	}
	return true;
}

/**
 *Setting radius on atom
 *@param[in] radius Is the radius to be set
 */
void Atom::SetRadius( double radius ){
	radius_ =  radius;
}

//UTILITIES
//CLEANERS
//DEBUG
/**
 *Printing atom parameters
 */
void Atom::PrintParameters(){
	char kind[ 10 ];
	if( this->GetKind( kind ) ){
	std::cout << "   correcto kind            : " << kind   << std::endl;
	}else{
	std::cout << "   incorrecto kind            : " << kind   << std::endl;
	}


}

/**
 *Printing atom information
 */
void Atom::PrintInfo(){
	this->PrintParameters();
	std::cout << "   Radius    	 : " << this->GetRadius() << std::endl;
	std::cout << "   Coordinates : " << this->GetCoordinate( 0 )  << " "; 
	std::cout << this->GetCoordinate( 1 ) << " " << this->GetCoordinate( 2 ) << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////////
//                              AMINOACID   METHODS                                     //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Aminoacid::Aminoacid(  ){
	kind_ = 20;
	atoms_.clear(  );
}

/**
 *Constructor receiving the kind of aminoacid
 *@param[in] kind  Is the kind of aminoacid that will be added
 *NOTE: The kind of aminoacid need to be a valid aminoacid (Se valids on capside.h file).
 */
Aminoacid::Aminoacid( int kind ){
	assert(  ( kind >= 0 )  &&  ( kind < 20 )  );
	kind_ = kind;
}

/**
 *Default destructor
 */
Aminoacid::~Aminoacid(  ){
	kind_ = 20;
	atoms_.clear(  );
}

//GETS
/**
 *Getting number of atoms that compose the aminoacid
 *@return A size_t value indicating the number of atoms contained in the aminoacid 
 */
size_t Aminoacid::GetNAtoms(  ){
	return atoms_.size();
}

/**
 *Getting kind of aminoacid
 *@return A int value indicating the kind of aminoacid processes (must be a number 
 *        between 0-19 )
 */
int Aminoacid::GetKind(  ){
	assert(  ( kind_ >= 0 )  &&  ( kind_ < 20 )  );
	return kind_;
}

/**
 *Getting diameter for a kind of atom
 *@param[in] kind Is a string containing the kind of atom
 *NOTE: The name of the atom can have more than a letter, in this case is only considered
 *      the first letter, each letter is associated to a number.
 *    1-C-Carbon
 *    2-H-Hydrogen
 *    3-O-Oxygen
 *    4-N-Nitrogen
 *    5-P-Phosphorus
 *    6-S-Sulfur
 *@return A double value indicating the diameter of the atom.
 */
double Aminoacid::GetAtomRadius( char* kind ){
	char temp = kind[ 0 ];
	if(  temp == ' '  ){
		temp = kind[ 1 ];
		if(  temp == ' '  ){
			temp = kind[ 2 ];
			if(  temp == ' '  ){
				temp = kind[ 3 ];
			}
		}
	}
	double real_radius;
	switch( temp ){
		case 'C':
			real_radius = radius_[ 0 ];
			break;
		case 'H':
			real_radius = radius_[ 1 ];
			break;
		case 'O':
			real_radius = radius_[ 2 ];
			break;
		case 'N':
			real_radius = radius_[ 3 ];
			break;
		case 'P':
			real_radius = radius_[ 4 ];
			break;
		case 'S':
			real_radius = radius_[ 5 ];
			break;
		default:
			std::cout << "Error: Atom not valid --" << kind << "--" << std::endl;		
			assert( 0 );
			break;
	}
	return real_radius;
}

/**
 *Getting coordinate from atom
 *@param[in] index Is the index of the atom on the array of atoms
 *@param[in] position Is the index of the axis where the coordinate is required
 *                    1-X axis
 *                    2-Y axis
 *                    3-Z axis
 *@return A double value with the coordinate requested
 */
double Aminoacid::GetAtomCoordinate( size_t index , int position ){
	//validating that the atom index exists
	assert( ( index >= 0 )  &&  ( index < atoms_.size() ) );
	return atoms_[ index ]->GetCoordinate( position );
}
	
/**
 *Getting coordinates from atom
 *@param[in] index Is the index of the atom on the array of atoms
 *@param[out] coords Is the array where the coordinates will be stored 
 */
void Aminoacid::GetAtomCoordinates( size_t index , double* coords ){
	//validating that the atom index exists
	assert( ( index >= 0 )  &&  ( index < atoms_.size() ) );
	atoms_[ index ]->GetCoordinates( coords );
}

/**
 *Getting atom from aminoacid
 *@param[in] index Is the index of the atom on the aminoacid
 *@return A pointer to the atom requested
 */
Atom* Aminoacid::GetAtom( size_t index ){
	assert( ( index >= 0 )  &&  ( index < atoms_.size() ) );
	return atoms_[ index ];
}

//SETS
/**
 *Setting the kind of aminoacid in the class
 *@param[in] kind Is the kind of aminoacid to be set
 *@return A bool value indicating if the value was set correctly
 */
bool Aminoacid::SetKind( int kind ){
	assert(  ( kind >= 0 )  &&  ( kind < 20 )  );
	kind_ = kind;
	return ( kind_ == kind );
}

/**
 *Setting coordinate on atom
 *@param[in] index Is the index of the atom where coordinate will be set
 *@param[in] coord Is the value of the coordinate to be set
 *@param[in] position Is the index of the axis where the coordinate is required
 *                    1-X axis
 *                    2-Y axis
 *                    3-Z axis
 *@return A bool value indicating if coordinate was set correctly
 */
bool Aminoacid::SetAtomCoordinate( size_t index , double coord , int position ){
	//validating that the atom index exists
	assert( ( index >= 0 )  &&  ( index < atoms_.size() ) );
	atoms_[ index ]->SetCoordinate( coord , position );
	return ( atoms_[ index ]->GetCoordinate( position ) == coord ); 
}

/**
 *Setting coordinates on atom
 *@param[in] index Is the index on the list of atoms
 *@param[in] coords Is the array with the atom coordinates
 *@return A bool value indicating if coordinates were added correctly 
 */
bool Aminoacid::SetAtomCoordinates( size_t index , double* coords ){
	//validating that the atom index exists
	assert( ( index >= 0 )  &&  ( index < atoms_.size() ) );
	atoms_[ index ]->SetCoordinates( coords );
	for(  int i_position = 0  ;  i_position < 3  ;  i_position++  ){
		if(  atoms_[ index ]->GetCoordinate( i_position ) != coords[ i_position ]  ){
			return false;
		}
	}
	return true;
}

//UTILITIES
/**
 *Adding an atom to the list of atoms forming the aminoacid
 *@param[in] kind Is the name of theatom to be added
 *@param[in] diameter Is the dimater eof the atom to be added
 *@param[in] coords Is the array withe spacial coordinates of the atom
 *@return A bool value indicating if the atom was added or not
 */
bool Aminoacid::AddAtom( char* kind , double diameter , double* coords ){
	size_t num_atoms = atoms_.size();
	Atom* new_atom = new Atom( kind , diameter , coords );
	atoms_.push_back( new_atom );
	return ( atoms_.size() > num_atoms );
}

/**
 *Testing if an aminoacid is valid
 *@param[in] kind Is a string containing the name of the aminoacid in order to validate
 *                If is one of the 20 valids
 *@param[out] position Is the index of the kind of aminoacid
 *@return A bool value indicating if the aminoacid is valid or not  
 */
bool Aminoacid::IsValid( char* kind , int*  position ){
	bool flag = false;
	(*position) = -1;
	for(  int i_kind = 0  ;  i_kind < 20  ;  i_kind++  ){
		if(  !strcmp( kind , valids_[ i_kind ] )  ){
			flag = true;
			(*position) = i_kind;
			break;
		}
	}
	return flag;
}

/**
 *Testing if the class have the same aminoacid than kind
 *@param[in] index Is the aminoacid to be tested
 *@return A bool value indicating if aminoacids are equals
 */
bool Aminoacid::AreEqual( int index ){
	return ( kind_ == index );
}

//CLEANERS
/**
 *Erasing all atoms from list of atoms
 *@return A bool value indicating if the list of atoms was cleaned
 */
bool Aminoacid::EmptyAtoms(  ){
	atoms_.clear();
	return ( atoms_.size() == 0 );
}

//DEBUG
/**
 *Printing parameters of aminoacid
 */
void Aminoacid::PrintParameters(){
	std::cout << "  kind            : " << this->GetKind()   << std::endl;
	std::cout << "  Number of atoms : " << this->GetNAtoms() << std::endl;
}

/**
 *Printing info contained on aminoacid
 */
void Aminoacid::PrintInfo(){
	this->PrintParameters();
	std::cout << "   ATOMS: " << this->GetNAtoms() << std::endl;
	for(  size_t i_atom = 0  ;  i_atom < this->GetNAtoms()  ;  i_atom++  ){
		std::cout << "   Atom id: " << i_atom << std::endl;
		Atom* aux = atoms_[ i_atom ];
		aux->PrintInfo();
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                PROTEIN   METHODS                                     //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Protein::Protein(  ){
	kind_ = 0;
	aminoacids_.clear(  );
}

/**
 *Constructor receiving the name of protein (this is a char value)
 *@param[in] kind Is the name of the protein to be added 
 */
Protein::Protein( char kind ){
	kind_ = kind;
}

/**
 *Default destructor
 */
Protein::~Protein(  ){
	kind_ = 0;
	aminoacids_.clear(  );
}

//GETS
/**
 *Getting the amount of aminoacids forming the protein
 *@return A size_t value indicating the amount of aminoacids forming the protein
 */
size_t Protein::GetNAminoacids(  ){
	return aminoacids_.size(  );
}

/**
 *Getting the name of the protein
 *@return A char value indicating the name of the protein
 */
char Protein::GetKind(  ){
	return kind_;
}

/**
 *Getting an aminoacid from a protein
 *@param[in] index Is the index on the vector containing the aminoacids
 *@return A pointer to Aminoacid class with the aminoacid requested
 */
Aminoacid* Protein::GetAminoacid( size_t index ){
	assert( index < aminoacids_.size() );
	return aminoacids_[ index ];
}

//SETS
/**
 *Setting kind of protein
 *@param[in] kind Is the name of the protein to be added
 *@return A bool value indicating if the kind was correctly set
 */
bool Protein::SetKind( char kind ){
	kind_ = kind;
	return ( kind_ == kind );
}

//UTILITIES
/**
 *Testing if an aminoacid is valid
 *@param[in] kind Is a string containing the name of the aminoacid in order to validate
 *                If is one of the 20 valids
 *@param[out] position Is the index of the kind of aminoacid
 *@return A bool value indicating if the aminoacid is valid or not  
 */
bool Protein::IsValid( char* kind , int*  position ){
	bool flag = false;
	(*position) = -1;
	for(  int i_kind = 0  ;  i_kind < 20  ;  i_kind++  ){
		if(  !strcmp( kind , valids_[ i_kind ] )  ){
			flag = true;
			(*position) = i_kind;
			break;
		}
	}
	return flag;
}

/**
 *Adding a new aminoacid
 *@param[in] kind Is the string containing the kind of aminoacid to be added
 *@return A bool value indicating if the aminoacid was added correctly
 */
bool Protein::AddAminoacid( char* kind ){
	int position;
	bool flag = false;
	if(  this->IsValid( kind , &position )  ){
		Aminoacid* new_aminoacid = new Aminoacid( position );
		aminoacids_.push_back( new_aminoacid );
		flag = true;
	}else{
		//In this case the aminoacid is not valid and the information was printed in the 
		//warnings file before, so this neve have to happen
		assert( 0 );
	}
	return flag;
}

//CLEANERS
/**
 *Cleaning the aminoacid list
 */
bool Protein::EmptyAminoacids(  ){
	aminoacids_.clear(  );
	return ( aminoacids_.size() == 0 );
}

//DEBUG
/**
 *Printing protein parameters
 */
void Protein::PrintParameters(){
	std::cout << " kind                 : " << this->GetKind()        << std::endl;
	std::cout << " Number of aminoacids : " << this->GetNAminoacids() << std::endl;
}

/**
 *Printing info contained on protein
 */
void Protein::PrintInfo(){
	this->PrintParameters();
	std::cout << "  AMINOACIDS: " << this->GetNAminoacids() << std::endl;
	for(  size_t i_amino = 0  ;  i_amino < this->GetNAminoacids()  ;  i_amino++  ){
		std::cout << "  Aminoacid id: " << i_amino << std::endl;
		Aminoacid* aux = aminoacids_[ i_amino ];
		aux->PrintInfo();
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                    VDB   METHODS                                     //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Vdb::Vdb(  ){
	strcpy( name_ , "\0" );
	ptr_ = NULL;
	mesh_ = NULL;
	proteins_.clear();
}

/**
 *Constructor receiving all parameters required
 *This constructor open the file required and set all the variables in the class
 *@param[in] name Is the name of the file to be opened
 */
Vdb::Vdb( char* name ){
	name_ = name;
	ptr_ = fopen( name_ , "r" );
	mesh_ = new Spheres;
	proteins_.clear( );
}

/**
 *Default destructor
 */
Vdb::~Vdb(  ){

	proteins_.clear(  );
	if(  ptr_  ){
		fclose( ptr_ );
		ptr_ = NULL;
	}
	strcpy( name_ , "\0" );
	delete mesh_;
}

//GETS
/**
 *Getting number of proteins in the capside
 *@return A size_t value with the amount of proteins in the capside
 */
size_t Vdb::GetNProteins(  ){
	return proteins_.size(  );
}

/**
 *Getting a protein pointer
 *@param[in] index Is the index of protein requested
 *@return A pointer to class protein
 */
Protein* Vdb::GetProtein( size_t index ){
	assert( ( index >= 0 ) && ( index < this->GetNProteins( ) ) );
	return proteins_[ index ];
}

FILE* Vdb::GetPointer(){
	return ptr_;
}

char* Vdb::GetName(){
	return name_;
}
//SETS
/**
 *Setting an atom on a capside
 *Need to be tested if the line readed is valid or not valid, if it is valid is added to 
 *the data added before or if data does not exist, then the classes are instanced
 *@param[in] params Is the set of parameters readed from file 
 */
void Vdb::SetAtomOnCapside( char params[11][500] ){

	if(  !strcmp( params[ 0 ] , "ATOM" )  ){
		//In this case the information correcponds to an ATOM, but some of the aminoacids are 
    //not valid so neet to be tested if the aminoacid is valid or not
		int position;
		//testing if is a valid aminoacid or not
		if(  this->AminoacidIsValid( params[ 3 ] , &position )  ){
			this->AddAtom( params , position );

		}else{//aminoacid not valid, printing message on warnings file
			char aux[1000];
			strcpy( aux , params[ 0 ] );
			strcat( aux , " " );
			for(  int i_par = 1  ;  i_par < 11  ;  i_par++  ){
				strcat( aux , params[ i_par ] );
				strcat( aux , " " );
			}
			strcat( aux , "\n" );
		}
	}else{//this line is auxiliary information
		char aux[1000];
		strcpy( aux , params[ 0 ] );
		strcat( aux , " " );
		for(  int i_par = 1  ;  i_par < 11  ;  i_par++  ){
			strcat( aux , params[ i_par ] );
			strcat( aux , " " );
		}
		strcat( aux , "\n" );
	}
}

//UTILITIES
/**
 *Adding a valid atom to the capside information
 *@param[in] params Is the line readed on file vdb
 */
void Vdb::AddAtom( char params[ 11 ][ 500 ] , int aminoacid_index ){
	if(  !this->GetNProteins(  )  ){//Not proteins added before
		char kind = params[ 4 ][ 0 ];
		this->AddProtein( kind );
		Protein* aux_prot = proteins_[ 0 ];
		aux_prot->AddAminoacid( params[ 3 ] );
		Aminoacid* aux_amino = aux_prot->GetAminoacid( 0 );
		double radius;
		double coords[ 3 ];
		radius = aux_amino->GetAtomRadius( params[ 2 ] );
		sscanf( params[ 6 ] , "%lf" , &coords[ 0 ] );
		sscanf( params[ 7 ] , "%lf" , &coords[ 1 ] );
		sscanf( params[ 8 ] , "%lf" , &coords[ 2 ] );
		aux_amino->AddAtom( params[ 2 ] , radius , coords );
	}else{//Tesst if atom belongs to the last protein added, if not add new protein 
		Protein* aux_prot = proteins_.back();
		char new_kind = params[ 4 ][ 0 ];
		char actual_kind = aux_prot->GetKind(  );
		if(  new_kind == actual_kind  ){//Both belongs to the same protein test if belongs to the same aminoacid
			Aminoacid* aux_amino = aux_prot->GetAminoacid( aux_prot->GetNAminoacids( ) - 1 );
			if(  aux_amino->AreEqual( aminoacid_index )  ){//the atom belongs to the last aminoacid
				double radius;
				double coords[ 3 ];
				radius = aux_amino->GetAtomRadius( params[ 2 ] );
				sscanf( params[ 6 ] , "%lf" , &coords[ 0 ] );
				sscanf( params[ 7 ] , "%lf" , &coords[ 1 ] );
				sscanf( params[ 8 ] , "%lf" , &coords[ 2 ] );
				aux_amino->AddAtom( params[ 2 ] , radius , coords );
			}else{//the atom belongs to a new aminoacid
				aux_prot->AddAminoacid( params[ 3 ] );
				Aminoacid* new_amino = aux_prot->GetAminoacid( aux_prot->GetNAminoacids( ) - 1 );
				double radius;
				double coords[ 3 ];
				radius = new_amino->GetAtomRadius( params[ 2 ] );
				sscanf( params[ 6 ] , "%lf" , &coords[ 0 ] );
				sscanf( params[ 7 ] , "%lf" , &coords[ 1 ] );
				sscanf( params[ 8 ] , "%lf" , &coords[ 2 ] );
				new_amino->AddAtom( params[ 2 ] , radius , coords );
			}
		}else{//this atom belongs to a new protein
			this->AddProtein( new_kind );
			Protein* new_prot = proteins_.back();
			new_prot->AddAminoacid( params[ 3 ] );
			Aminoacid* aux_amino = new_prot->GetAminoacid( 0 );
			double radius;
			double coords[ 3 ];
			radius = aux_amino->GetAtomRadius( params[ 2 ] );
			sscanf( params[ 6 ] , "%lf" , &coords[ 0 ] );
			sscanf( params[ 7 ] , "%lf" , &coords[ 1 ] );
			sscanf( params[ 8 ] , "%lf" , &coords[ 2 ] );
			aux_amino->AddAtom( params[ 2 ] , radius , coords );
		}
	}
}

/**
 *Adding a new protein
 *@param[in] kind Is the name of the new protein
 *@return a bool value indicating if the protein was added or not
 */
bool Vdb::AddProtein( char kind ){
	size_t num_proteins = proteins_.size(  );
	Protein* new_protein = new Protein( kind );
	proteins_.push_back( new_protein );
	return ( proteins_.size() > num_proteins );
}

/**
 *Reading the complete vdb file
 */
void Vdb::ReadCompleteFile(  ){
	FILE* fp = this->GetPointer();
	while( !feof( fp ) ){
		int read;
		read = this->ReadLine( );
		if( read == EOF ){
			break;
		}
	}
	this->FillSphericalMeshInformation(  );

}

/**
 *Reading a line from file
 *This method adds the required information on the data structure
 */
int Vdb::ReadLine( ){
	char params[ 11 ][ 500 ];
	FILE* fp = this->GetPointer();
	char line[500];
	fgets( line , 500 , fp );
	//if is end of file does not do anyting
	if(  feof( fp ) ){
		return feof( fp );
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

	this->SetAtomOnCapside( params );

	return 1;
}

/**
 *Testing if an aminoacid is valid
 *@param[in] kind Is a string containing the name of the aminoacid in order to validate
 *                If is one of the 20 valids
 *@param[out] position Is the index of the kind of aminoacid
 *@return A bool value indicating if the aminoacid is valid or not  
 */
bool Vdb::AminoacidIsValid( char* kind , int*  position ){
	bool flag = false;
	(*position) = -1;
	for(  int i_kind = 0  ;  i_kind < 20  ;  i_kind++  ){
		if(  !strcmp( kind , valids_[ i_kind ] )  ){
			flag = true;
			(*position) = i_kind;
			break;
		}
	}
	return flag;
}

/**
 *Filling the spherical mesh information based on the atoms information of the capside
 */
void Vdb::FillSphericalMeshInformation(  ){
	size_t n_proteins = this->GetNProteins(  );
	for(  size_t i_prot = 0  ;  i_prot < n_proteins  ;  i_prot++  ){
		Protein* aux_prot = proteins_[ i_prot ];
		size_t n_aminoacids = aux_prot->GetNAminoacids(  );
		for(  size_t i_amino = 0  ;  i_amino < n_aminoacids  ;  i_amino++  ){
			Aminoacid* aux_amino = aux_prot->GetAminoacid( i_amino );
			size_t n_atoms = aux_amino->GetNAtoms(  );
			for(  size_t i_atom = 0  ;  i_atom < n_atoms  ;  i_atom++  ){
				Atom* aux_atom = aux_amino->GetAtom( i_atom );
				double coords[ 3 ], radius;
				int material;
				aux_atom->GetCoordinates( coords );
				radius = aux_atom->GetRadius(  );
				material = aux_atom->GetMaterialIndex(  );
				Sphere *aux = new Sphere( coords , radius , material );
				mesh_->AddElement( aux );
			}
		}
	}
}

/**
 *Saving the mesh of spheres representing the atoms
 */
void Vdb::SaveGiDSphericalMesh(  ){
	mesh_->SaveMeshOnGiDFile( (char*)"capside_mesh.msh" );
}

/**
 *Scaling the spherical mesh
 */
void Vdb::ScaleSpheresMesh(  ){
	//First the spheres mesh is rotated in order to have the fold 5 pointing to the Y axis 

	//mesh_->Rotate(  );

  mesh_->FillScaler( );
  mesh_->ScaleMesh( );
	//Assigning scaled mesh to atoms information
	size_t index = 0;
	size_t n_proteins = this->GetNProteins(  );
	for(  size_t i_prot = 0  ;  i_prot < n_proteins  ;  i_prot++  ){
		Protein* aux_prot = proteins_[ i_prot ];
		size_t n_aminoacids = aux_prot->GetNAminoacids(  );
		for(  size_t i_amino = 0  ;  i_amino < n_aminoacids  ;  i_amino++  ){
			Aminoacid* aux_amino = aux_prot->GetAminoacid( i_amino );
			size_t n_atoms = aux_amino->GetNAtoms(  );
			for(  size_t i_atom = 0  ;  i_atom < n_atoms  ;  i_atom++  ){
				Atom* aux_atom = aux_amino->GetAtom( i_atom );
				aux_atom->SetCoordinate( mesh_->GetCoord( index , 0 ) , 0 );
				aux_atom->SetCoordinate( mesh_->GetCoord( index , 1 ) , 1 );
				aux_atom->SetCoordinate( mesh_->GetCoord( index , 2 ) , 2 );
 				aux_atom->SetRadius(  mesh_->GetRadius( index  )  );
				index++;
			}
		}
	}
}

/**
 *Llevando las coordenadas del octree al dominio original de la capside
 *@param[in] i_pos Is the position X,Y,Z of coordinate
 *@param[in] coord  Is the coordinate to be scaled
 *@return A double value with the coordinate in the original domain
 */
double Vdb::UnscaleCoordinate( int i_pos , double coord ){
	double value;
	value = mesh_->UnscaleCoordinate( i_pos , coord );
	return value;
}

/**
 *Printing amount of proteins aminoacids and atoms on capside
 */
void Vdb::PrintResumenOnWarnings( ){

	size_t t_aminoacids = 0 , t_atoms = 0;
	size_t n_proteins = this->GetNProteins(  );
	size_t* n_aminoacids;
	n_aminoacids = new size_t [ n_proteins ];
	for(  size_t i_prot = 0 ;  i_prot < n_proteins  ;  i_prot++  ){
		n_aminoacids[ i_prot ] = 0;
	}
	size_t** n_atoms;
	n_atoms = new size_t* [ n_proteins ];


	for(  size_t i_prot = 0  ;  i_prot < n_proteins  ;  i_prot++  ){
		Protein* aux_prot = proteins_[ i_prot ];
		n_aminoacids[ i_prot ] = aux_prot->GetNAminoacids(  );
		t_aminoacids += n_aminoacids[ i_prot ];
		n_atoms[ i_prot ] = new size_t [ n_aminoacids[ i_prot ] ];
		for(  size_t i_amino = 0  ;  i_amino < n_aminoacids[ i_prot ]  ;  i_amino++  ){
			Aminoacid* aux_amino = aux_prot->GetAminoacid( i_amino );
			n_atoms[ i_prot ][ i_amino ] = aux_amino->GetNAtoms(  );
			t_atoms += n_atoms[ i_prot ][ i_amino ];
		}
	}
}

/**
 *Scaling coordinate to the normalized space
 *@param[in] i_pos Is the index of the dimention that will be scaled 0-X,1-Y,2-Z
 *@param[in] coord Is the coordinate to be scaled
 *@return A double value with the coordinate scaled in the domain 0-1
 */
double Vdb::ScaleCoord( int i_pos , double coord ){
	double coordinate = mesh_->ScaleCoord( i_pos , coord );
	return coordinate;
}

//CLEANERS
/**
 *Cleaning proteins of file
 *@return A bool value indicating if the vector is clean or not
 */
bool Vdb::EmptyProteins(  ){
	proteins_.clear(  );
	return ( proteins_.size(  ) == 0 );
}
		
//DEBUG
/**
 *Printing parameters of file
 */
void Vdb::PrintParameters(){
	std::cout << "Name of input file   : " << this->GetName()      << std::endl;
	std::cout << "Pointer to file      : " << this->GetPointer()   << std::endl;
	std::cout << "Number of proteins   : " << proteins_.size()     << std::endl;
}

/**
 *Printing all information from 
 */
void Vdb::PrintInfo(){
	std::cout << "Printing information contained on Vdb file " << std::endl;
	this->PrintParameters();
	std::cout << "PROTEINS: " << this->GetNProteins() << std::endl;
	for(  size_t i_prot = 0  ;  i_prot < this->GetNProteins()  ;  i_prot++  ){
		std::cout << " Protein id: " << i_prot << std::endl;
		Protein* aux = proteins_[ i_prot ];
		aux->PrintInfo();
	}
}










