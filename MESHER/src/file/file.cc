#include "file.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                     FILE METHODS                                     //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
File::File(){
	name_ = new char[100];
	strcpy( name_ , "\0" );
	ptr_ = NULL;
	direction_ = 0;
	type_ = 0;
}

/**
 *Constructor receiving all parameters required
 *This constructor open the file required and set all the variables in the class
 *@param[in] name Is the name of the file to be opened
 *@param[in] direction Is the direction of the file 1=input, 2=output
 *@param[in] type Is the kind of file 1=ASCII, 2=binary
 */
File::File( char* name , int direction , int type ){
	assert( ( direction > 0 ) && ( direction < 3 )  );
	assert( ( type > 0 ) && ( type < 3 )  );
	name_ = name;
  direction_ = direction;
	type_ = type;
	if(  direction == 1  ){//Is to read a file
		if(  type == 1  ){//Is a plane text file
			ptr_ = fopen( name , "r" );
		}else{//Is a binary file
			ptr_ = fopen( name , "rb" );
		}
	}else{//Is to write a file
		if(  type == 1  ){//Is a plane text file
			ptr_ = fopen( name , "w" );
		}else{//Is a binary file
			ptr_ = fopen( name , "wb" );
		}
	}
}

/**
 *Default destructor
 */
File::~File(){
	if(  ptr_  ){
		fclose( ptr_ );
	}
	direction_ = 0;
	type_ = 0;
	ptr_ = NULL;
}
	
//GETS
/**
 *Getting pointer to file openend
 *@return A FILE pointer
 */
FILE* File::GetPointer( ){
	return ptr_;
}

/**
 *Getting name of the file opened
 *@return A char string with the name of the file opened
 */
char* File::GetName( ){
	return name_;
}

/**
 *Getting diretion of the file, 0-No direction, 1-input, 2-output
 *@return A int value indicating the direction
 */
int File::GetDirection(  ){
	return direction_;
}

/**
 *Getting type of file opened, 0-Not specified, 1-ASCII, 2-binary
 *@return A int value inticating the kind of file.
 */
int File::GetType(  ){
	return type_;
}

//SETS
/**
 *Setting name of the file to be opened
 *@param[in] name Is the name of the file that will be opened
 *@return A bool value indicating if name was set
 */
bool File::SetName( char* name ){
	name_ = name;
	return ( name_ == name );
}

/**
 *Setting file direction
 *@param[in] direction Is a int value that indicates if file is input or output
 *@return A bool value indicating if direction was setted correctly
 */
bool File::SetDirection( int direction ){
	direction_ = direction;
	return ( direction_ == direction );
}

/**
 *Setting type of file
 *@param[in] type Is a int value indicating if type of file is ASCII or binary
 *@return A bool value indicating if kind of file was set correctly
 */
bool File::SetType( int type ){
	type_ = type;
	return ( type_ == type );
}

/**
 *Setting pointer to file on null pointer
 *@return A bool value indicating if pointer was set correctly
 */
bool File::SetPointerToNull(  ){
	ptr_ = NULL;
	return ( !ptr_ );
}

/**
 *Setting pointer to file opened
 *@param[in] ptr Is the pointer to the file opened
 *@return A bool value indicating if pointer was set correctly
 */
bool File::SetPointer( FILE* ptr ){
	ptr_ = ptr;
	return ( ptr_ == ptr );
}

//UTILITIES
/**
 *Making a copy of the file
 *NOTE: In order to avoid wrong writing or reading the name of the copied file need to be 
 *      different, if have the same name, is added "copy" at the end of the name. Also 
 *      the file to be copied need to be closed.
 *@param[in] copy_name Is the name of the file copied
 *@return a bool value indicating if the copy h=was performed succesfully
 */
bool File::MakeCopy( char* copy_name ){
	//Assuring different name

	//assuring file closed
	
	assert( 0 );//Not implemented yet
	return true;
}

/**
 *Changing the file name 
 *@param[in] new_name Is the new name of the file.
 *@return The new name of the file
 */
char* File::ChangeName( char* new_name ){
	
	assert( 0 );//Not implemented yet
	return name_;
}

/**
 *Closing file
 *@return A bool value indicating if the file was closed correctly
 */
bool File::CloseFile(  ){
	if( ptr_ ){
		fclose( ptr_ );
		ptr_ = NULL;
	}
	return ptr_;
}

//CLEANERS
/*
 *Erasing the contents of a file.
 *@return A bool value indicating if the file was erased correctly.
 */
bool File::CleanFile(  ){

	assert( 0 ); //Not implemented yet
	return true;
}

//DEBUG
/**
 *Searching a string in file
 *@param[in] str Is the string to be searched
 *@param[out] row Indicates the row where the string begins
 *@param[out] col Indicates the column where the strinf begins
 *@return A bool value indicating if string was founded or not.
 */
bool File::SearchString( char* str , int* row, int* col ){
	
	assert( 0 );//Not implemented yet
	return true;
}








