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

	char params[ 11 ][ 500 ];
	
	FILE *fin, *fout;
	fin = fopen( argv[ 1 ] , "r" );
	fout = fopen( argv[ 2 ] , "w" );
	 
	printf(" - Will extract aminoacid ATOMs from %s and save then in %s \n", argv[1],argv[2]);

	char line[128];
	
	char aminoacids[]="ALA /*Alanine*/,  ARG /*Arginine*/,  ASN /*Asparagine*/,  ASP /*Aspartic acid*/,  CYS /*Cysteine*/,  GLU /*Glutamic acid*/,  GLN /*Glutamine*/,  GLY /*Glycine*/,  HIS /*Histidine*/,  ILE /*Isoleucine*/,  LEU /*Leucine*/,  LYS /*Lysine*/,  MET /*Methionine*/,  PHE /*Phenylalanine*/,  PRO /*Proline*/,  SER /*Serine*/,  THR /*Threonine*/,  TRP /*Tryptophan*/,  TYR /*Tyrosine*/,  VAL /*Valine*/";

	while( !feof( fin ) ){
		fgets( line , 128 , fin );

	//reading column 1
	//reading column 1 Record name ATOM 1-6
	int position = 0;
	for(  int i_char = 0  ;  i_char < 6  ;  i_char++  ){
		params[ 0 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 0 ][ 6 ] = 0;
	
	//reading column 2
	//reading column 2 Atom serial number 7-11
	position = 6;
	for(  int i_char = 0  ;  i_char < 5  ;  i_char++  ){
		  params[ 1 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 1 ][ 5 ] = 0;
	int AtomSerialNumber = atoi(params[ 1 ]);

	//reading column 3
	//reading column 3 Atom name 13-16
	position = 12;
	for(  int i_char = 0  ;  i_char < 4  ;  i_char++  ){
		params[ 2 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 2 ][ 4 ] = 0;
 
	//reading column 4
	//reading column 4 Resiude name 18-20
	position = 17;
	for(  int i_char = 0  ;  i_char < 3  ;  i_char++  ){
		params[ 3 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 3 ][ 3 ] = 0;

	//reading column 5
	//reading column 5 Chain identifier 22
	position = 21;
	for(  int i_char = 0  ;  i_char < 1  ;  i_char++  ){
		params[ 4 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 4 ][ 1 ] = 0;

	//reading column 6
	//reading column 6 Residue sequence number 23-26
	position = 22;
	for(  int i_char = 0  ;  i_char < 4  ;  i_char++  ){
		params[ 5 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 5 ][ 4 ] = 0;

	//reading column 7
	//reading column 7 Orthogonal coordinates for X in Angstroms 31-38
	position = 30;
	for(  int i_char = 0  ;  i_char < 8  ;  i_char++  ){
		params[ 6 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 6 ][ 8 ] = 0;
	double coordX=strtod(params[ 6 ],NULL);

	//reading column 8
	//reading column 8 Orthogonal coordinates for Y in Angstroms 39-46
	position = 38;
	for(  int i_char = 0  ;  i_char < 8  ;  i_char++  ){
		params[ 7 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 7 ][ 8 ] = 0;
	double coordY=strtod(params[ 7 ],NULL);

	//reading column 9
	//reading column 9 Orthogonal coordinates for Z in Angstroms 47-54
	position = 46;
	for(  int i_char = 0  ;  i_char < 8  ;  i_char++  ){
		params[ 8 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 8 ][ 8 ] = 0;
	double coordZ=strtod(params[ 8 ],NULL);

	//reading column 10
	//reading column 10 Occupancy 55-60
	position = 54;
	for(  int i_char = 0  ;  i_char < 6  ;  i_char++  ){
		params[ 9 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 9 ][ 6 ] = 0;

	//reading column 11
	//reading column 11 Temperature factor 61-66
	position = 60;
	for(  int i_char = 0  ;  i_char < 6  ;  i_char++  ){
		params[ 10 ][ i_char ] = line[ position ];
		position++;
	}
	params[ 10 ][ 6 ] = 0;

	if(strstr(params[0],"ATOM")){ // extract only RECORD ATOM 
	  if(strstr(aminoacids,params[3])){ // extract only aminoacids
	     //printf("%4s%5d %4s %3s %s%s    %8.3f%8.3f%8.3f%6s%6s \n", params[ 0 ], AtomSerialNumber, params[ 2 ], params[ 3 ], params[ 4 ], params[ 5 ], coordX, coordY, coordZ, params[ 9 ], params[ 10 ]);
             fprintf( fout , "%s" , line );
	  }
	}


	} // end of while( !feof( fin ) )

	/*
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
	*/
	
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


