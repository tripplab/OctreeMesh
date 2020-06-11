//C system files
#include <math.h>
#include <assert.h>
#include <string.h>

//C++ system files
#include <vector>
#include <iostream>
#include <cstdio>
/**
 *This program takes a vdb file and proccess its information in order to eliminate the 
 *rows that contain irrelevant information to generate FEM meshes
 *Receives two arguments
 *@param[in] name of vdb file to be cleaned
 *@param[in] name of file to save the relevant information for generating the mesh
 */

void 	CheckParameters( int argc );


int main( int argc , char** argv ){
	
	CheckParameters( argc );

	FILE *fin, *fout;
	fin = fopen( argv[ 1 ] , "r" );
	fout = fopen( argv[ 2 ] , "w" );

	char line[128];
	
	while( !feof( fin ) ){
		fgets( line , 128 , fin );
		if(  line[ 0 ] == 'A'  ){
			if(  line[ 1 ] == 'T'  ){
				if(  line[ 2 ] == 'O'  ){
					if(  line[ 3 ] == 'M'  ){
						fprintf( fout , "%s" , line );
					}
				}
			}			
		}
		

	}
	
	fclose( fin );
	fclose( fout );
  return 0;
}

/**
 *Checking if number of arguments received is correct or not
 */
void CheckParameters( int argc ){
	if(  argc != 3  ){
		std::cout << "\n\n";
		std::cout << "**************************************************************************************" << std::endl;
		std::cout << "*********************************CLEAN PDB APPLICATION********************************" << std::endl;
		std::cout << "Received wrong number of arguments " << std::endl; 
		std::cout << "Arguments must be: " << std::endl;
		std::cout << "  --[1][string][input ] Name of vdb file containing the atoms information" << std::endl;
		std::cout << "  --[2][string][output] Name of vdb file to save the relevant information for meshing the biomeloceule" << std::endl;
		std::cout << "*********************************CLEAN PDB APPLICATION********************************" << std::endl;
		std::cout << "**************************************************************************************" << std::endl;
		std::cout << "\n\n";
		assert( 0 );
	}
}


