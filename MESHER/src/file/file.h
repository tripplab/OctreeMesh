#pragma once

//C system files
#include <math.h>
#include <assert.h>
#include <string.h>

//C++ system files
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>

//Other includes


/**
 *File class
 *This class is used to manage files in the execution
 *@param name_ Is the name of the file
 *@param ptr_ Is the pointer to the file
 *@param direction_ Is the direction to work with the file:
 * 0-Not indicated
 * 1-Inout
 * 2-Output
 *@param type_ Is the kind of file to be managed
 * 0-Not indicated
 * 1-ASCII file
 * 2-Binary file
 */
class File{

	char* name_;
	FILE* ptr_;
	int direction_;
	int type_;
		
	public:
		//CONSTRUCTOR AND DESTRUCTOR
		File();
		File( char* name , int direction , int type );
		~File();

		//GETS
		FILE* GetPointer( );
		char* GetName( );
		int GetDirection(  );
		int GetType(  );

		//SETS
		bool SetName( char* name );
		bool SetDirection( int direction );
		bool SetType( int type );
		bool SetPointerToNull(  );
		bool SetPointer( FILE* ptr );

		//UTILITIES
		bool MakeCopy( char* copy_name );
		char* ChangeName( char* new_name );
		bool CloseFile(  );

		//CLEANERS
		bool CleanFile(  );

		//DEBUG
		bool SearchString( char* str , int* row, int* col );

};
