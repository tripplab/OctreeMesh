//C system files
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//C++ system files
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

using namespace std;

/**
 *This program takes a geometry.dat file and rotates all its coordinates over X axis 90
 *degrees 
 *	This program receives 2 arguments
 *		--Name of input file
 *		--Name of output file
 */

void CalculateVector( double* A , double* B , double* vec ){
	vec[ 0 ] = B[ 0 ] - A[ 0 ]; 
	vec[ 1 ] = B[ 1 ] - A[ 1 ]; 
	vec[ 2 ] = B[ 2 ] - A[ 2 ]; 
}

double Det3( double* A , double* B , double* C ){
	double det = 0.0;
	det -=    A[ 0 ]   * ( ( B[ 1 ]*C[ 2 ] ) - ( B[ 2 ] * C[ 1 ] ) );	
	det -= ( -A[ 1 ] ) * ( ( B[ 0 ]*C[ 2 ] ) - ( B[ 2 ] * C[ 0 ] ) );
	det -=    A[ 2 ]   * ( ( B[ 0 ]*C[ 1 ] ) - ( B[ 1 ] * C[ 0 ] ) );
	return det;
}

double VoluemenHexaedro( double** nodes , uint* element ){
	double volumen;
	uint T_indexes[ 24 ][ 4 ] = { { 0,1,0,0 } , { 1,2,0,0 } , { 2,3,0,0 } , { 3,0,0,0 } ,
													 		  { 0,3,1,0 } , { 3,7,1,0 } , { 7,4,1,0 } , { 4,0,1,0 } ,
													 		  { 0,4,2,0 } , { 4,5,2,0 } , { 5,1,2,0 } , { 1,0,2,0 } ,
													 		  { 1,5,3,0 } , { 5,6,3,0 } , { 6,2,3,0 } , { 2,1,3,0 } ,
													 		  { 2,6,4,0 } , { 6,7,4,0 } , { 7,3,4,0 } , { 3,2,4,0 } ,
													 		  { 6,5,5,0 } , { 5,4,5,0 } , { 4,7,5,0 } , { 7,6,5,0 } };
	uint F_indexes[ 6 ][ 4 ] = { { 0,1,2,3 } , { 0,3,4,7 } , { 0,1,5,4 } , { 1,2,5,6 } , 
															 { 2,3,6,7 } , { 4,5,6,7 } };
	double FC[ 6 ][ 3 ];
	double C[ 3 ];
	//calculando el centro de cada cara
	for(  uint i_face = 0  ;  i_face < 6  ;  i_face++  ){
		FC[ i_face ][ 0 ] = 0.0;
		FC[ i_face ][ 1 ] = 0.0;
		FC[ i_face ][ 2 ] = 0.0;
		for(  uint i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			for(  uint i_node = 0  ;  i_node < 4  ;  i_node++  ){
				FC[ i_face ][ i_dim ] += nodes[ element[ F_indexes[ i_face ][ i_node ] ] ][ i_dim ];
			}
		}
		FC[ i_face ][ 0 ] /= 4.0;
		FC[ i_face ][ 1 ] /= 4.0;
		FC[ i_face ][ 2 ] /= 4.0;
	}
	//calculando el centro del hexaedro
	C[ 0 ] = 0.0;
	C[ 1 ] = 0.0;
	C[ 2 ] = 0.0;
	for(  uint i_node = 0  ;  i_node < 8  ;  i_node++  ){
		for(  uint i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			C[ i_dim ] += nodes[ element[ i_node ] ][ i_dim ];
		}
	}
	C[ 0 ] /= 8.0;
	C[ 1 ] /= 8.0;
	C[ 2 ] /= 8.0;
	
	volumen = 0.0;
	for(  uint i_tet = 0  ;  i_tet < 24  ;  i_tet++  ){
		double A[ 3 ] , B[ 3 ] , CC[ 3 ], D[ 3 ], DA[ 3 ] , DB[ 3 ] , DC[ 3 ];
		for(  uint i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
			A[ i_dim ]  = nodes[ element[ T_indexes[ i_tet ][ 0 ] ] ][ i_dim ];
			B[ i_dim ]  = nodes[ element[ T_indexes[ i_tet ][ 1 ] ] ][ i_dim ];
			CC[ i_dim ] = FC[ 						T_indexes[ i_tet ][ 2 ] ][ i_dim ];
			D[ i_dim ] = C[ i_dim ];
		}
		CalculateVector( D , A , DA );
		CalculateVector( D , B , DB );
		CalculateVector( D , CC , DC );
		volumen += ( 0.1666666 * Det3( DA , DB , DC ) );
	}
	return volumen;
}

double CalculateVolumen( uint n_elements , double** nodes , uint** elements ){
	double volumen = 0.0;
	for(  uint i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		volumen += VoluemenHexaedro( nodes , elements[ i_elem ] );
	}
	return volumen;
}

/**
 *Performing matrix by vector operation
 *@param[in] mat Is the matrix to be multiplied
 *@param[in] rows Is the number of rows in the matrix
 *@param[in] cols Is the number of cols in the matrix
 *@param[in] vec Is the vector by the one the matrix will be mutiplied
 *@param[in] dim Is the length of the vector
 *@param[out] res is the vector where the matrix vector multiplication will be stored
 *@param[in] res_dim Is the number of elements of the results vector
 */
void MatrixVectorMultiplication( double** mat , int rows , int cols , double* vec , 
																 int dim , double* res , int res_dim ){
	//Validating dimension
	assert( rows == res_dim );
	assert( cols == dim );
	for(  int i_row = 0  ;  i_row < rows  ;  i_row++  ){
		res[ i_row ] = 0.0;
		for(  int i_col = 0  ;  i_col < cols  ;  i_col++  ){
			res[ i_row ] += ( mat[ i_row ][ i_col ] * vec[ i_col ] );
		}
	}
}


void RotateFold2Id0( double* coord ){

}
/**
 *Rotating coordinate to put id 0 from fold 2 aligned with Y axis
 *The coordinate is rotated over the X axis 90 degrees (1.570796327 radians)
 */
void RotateNodes( size_t n_nodes , double** nodes ){
	double** matrix;
	matrix = new double*[ 3 ];
	matrix[ 0 ] = new double[ 3 ];
	matrix[ 1 ] = new double[ 3 ];
	matrix[ 2 ] = new double[ 3 ];
	double newcoord[ 3 ] , teta;
	teta  = 1.570796327;
	matrix[ 0 ][ 0 ] = 1.0; matrix[ 0 ][ 1 ] = 0.0;         matrix[ 0 ][ 2 ] = 0.0; 
	matrix[ 1 ][ 0 ] = 0.0; matrix[ 1 ][ 1 ] = cos( teta ); matrix[ 1 ][ 2 ] = -sin( teta );
	matrix[ 2 ][ 0 ] = 0.0; matrix[ 2 ][ 1 ] = sin( teta ); matrix[ 2 ][ 2 ] = cos( teta );

	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		MatrixVectorMultiplication( matrix , 3 , 3 , nodes[ i_node ] , 3 , newcoord , 3 );	
		nodes[ i_node ][ 0 ] = newcoord[ 0 ];	nodes[ i_node ][ 1 ] = newcoord[ 1 ];	nodes[ i_node ][ 2 ] = newcoord[ 2 ];
	}

	delete matrix[ 0 ];
	delete matrix[ 1 ];
	delete matrix[ 2 ];
	delete matrix;	

}



/**
 *Esta funcion lee el numero de nodos y elementos del archivo de entrada
 *@param[in] name Es el nombre del archivo que ser\'a leido
 *@param[out] n_nodes Es la referencia a la variable donde se almacena la cantidad de nodos
 										  que tiene la malla
 *@param[ou] n_elements Es la referencia a la variable donde se almacena la cantidad de 
 *											elementos que tiene la malla
 */
void ReadNNodesAndNElements(  char* name , size_t* n_nodes , size_t* n_elements  ){

	FILE* ptr = NULL;

	ptr = fopen( name , "r" );
	if(  !ptr  ){
		std::cout << "File could not be opened---PROGRAM FINISHED" << name << std::endl;
		assert( 0 );
	}
	int nodos, elementos;
	char buffer[ 1000 ];
	for(  size_t i_trash = 0  ;  i_trash < 8  ;  i_trash++  ){
		fscanf( ptr , "%s" , buffer );
	}	
	nodos = atoi( buffer );
	(*n_nodes) = nodos;
	for(  size_t i_trash = 0  ;  i_trash < 7  ;  i_trash++  ){
		fscanf( ptr , "%s" , buffer );
	}
	for(  size_t i_trash = 0  ;  i_trash < (*n_nodes)  ;  i_trash++  ){
		for(  size_t j_trash = 0  ;  j_trash < 3  ;  j_trash++  ){
			fscanf( ptr , "%s" , buffer );
		}
	}
	for(  size_t i_trash = 0  ;  i_trash < 15  ;  i_trash++  ){
		fscanf( ptr , "%s" , buffer );
	}
	elementos = atoi( buffer );
	(*n_elements) = elementos;

	fclose( ptr );

	/*ifstream file;
	file.open( name );
	if(  !file.is_open()  ){
		std::cout << "File could not be opened---PROGRAM FINISHED" << name << std::endl;
		assert( 0 );
	}
	std::string word;
	for(  size_t i_trash = 0  ;  i_trash < 8  ;  i_trash++  ){
		file >> word;	
	}
	(*n_nodes) = stoi( word );
	for(  size_t i_trash = 0  ;  i_trash < 7  ;  i_trash++  ){
		file >> word;	
	}
	for(  size_t i_trash = 0  ;  i_trash < (*n_nodes)  ;  i_trash++  ){
		for(  size_t j_trash = 0  ;  j_trash < 3  ;  j_trash++  ){
			file >> word;	
		}
	}
	for(  size_t i_trash = 0  ;  i_trash < 15  ;  i_trash++  ){
		file >> word;	
	}

	(*n_elements) = stoi( word );

	file.close();*/
}

/**
 *Este metodo lee las coordenadas de los nodos y los elementos de la malla
 *@param[in] name Es el nombre del archivo donde estan guardados los datos de la malla
 *@param[in] n_nodes Es el numero de nodos que tiene la malla
 *@param[in] n_elements Es el numero de elementos que pertenecen a la malla
 *@param[out] nodes Es el arreglo donde se almacenar\'an las coordenadas de los nodos
 *@param[out] elements Es el arreglo donde se almacenar\'an los elementos de la malla
 *@param[out] material Es el arreglo donde se almacenan los materiales de los elemenotos
 */
void ReadNodesAndElements( char* name , size_t n_nodes , size_t n_elements , 
													 double** nodes , size_t** elements , int* material ){
	FILE* ptr = NULL;	
	ptr = fopen( name , "r" ); 
	if(  !ptr  ){
		std::cout << "File could not be opened---PROGRAM FINISHED" << name << std::endl;
		assert( 0 );
	}
	char buffer[1000];
	for(  size_t i_trash = 0  ;  i_trash < 15  ;  i_trash++  ){
		fscanf(ptr,"%s",buffer);	
	}

	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		double coord[ 3 ];
		for(  size_t i_pos = 0  ;  i_pos < 3  ;  i_pos++){
			fscanf(ptr,"%s",buffer);
			coord[ i_pos ] = atof( buffer );
		}	
		nodes[ i_node ][ 0 ] = coord[ 0 ];
		nodes[ i_node ][ 1 ] = coord[ 1 ];
		nodes[ i_node ][ 2 ] = coord[ 2 ];
	}

	for(  size_t i_trash = 0  ;  i_trash < 23  ;  i_trash++  ){
		fscanf(ptr,"%s",buffer);	
	}
	
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		size_t elem[ 8 ];
		fscanf(ptr,"%s",buffer);	
		material[ i_elem ] = atoi( buffer );
		for(  size_t i_pos = 0  ;  i_pos < 8  ;  i_pos++){
			fscanf(ptr,"%s",buffer);	
			elem[ i_pos ] = atoi( buffer );
		}	
		elements[ i_elem ][ 0 ] = elem[ 0 ];
		elements[ i_elem ][ 1 ] = elem[ 1 ];
		elements[ i_elem ][ 2 ] = elem[ 2 ];
		elements[ i_elem ][ 3 ] = elem[ 3 ];
		elements[ i_elem ][ 4 ] = elem[ 4 ];
		elements[ i_elem ][ 5 ] = elem[ 5 ];
		elements[ i_elem ][ 6 ] = elem[ 6 ];
		elements[ i_elem ][ 7 ] = elem[ 7 ];
	}

	fclose( ptr );

	/*ifstream file;
	file.open( name );
	if(  !file.is_open()  ){
		std::cout << "File could not be opened---PROGRAM FINISHED" << name << std::endl;
		assert( 0 );
	}
	std::string word;
	for(  size_t i_trash = 0  ;  i_trash < 15  ;  i_trash++  ){
		file >> word;	
	}

	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		double coord[ 3 ];
		for(  size_t i_pos = 0  ;  i_pos < 3  ;  i_pos++){
			file >> word;
			coord[ i_pos ] = stod( word );
		}	
		nodes[ i_node ][ 0 ] = coord[ 0 ];
		nodes[ i_node ][ 1 ] = coord[ 1 ];
		nodes[ i_node ][ 2 ] = coord[ 2 ];
	}

	for(  size_t i_trash = 0  ;  i_trash < 23  ;  i_trash++  ){
		file >> word;	
	}

	for(  size_t i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		size_t elem[ 8 ];
		file >> word;
		material[ i_elem ] = stoi( word );
		for(  size_t i_pos = 0  ;  i_pos < 8  ;  i_pos++){
			file >> word;
			elem[ i_pos ] = stoi( word );
		}	
		elements[ i_elem ][ 0 ] = elem[ 0 ];
		elements[ i_elem ][ 1 ] = elem[ 1 ];
		elements[ i_elem ][ 2 ] = elem[ 2 ];
		elements[ i_elem ][ 3 ] = elem[ 3 ];
		elements[ i_elem ][ 4 ] = elem[ 4 ];
		elements[ i_elem ][ 5 ] = elem[ 5 ];
		elements[ i_elem ][ 6 ] = elem[ 6 ];
		elements[ i_elem ][ 7 ] = elem[ 7 ];
	}

	file.close();*/
}

/**
 *Funcion para guardar la malla de coordenadas rotadas
 *@param[in] name Es el nombre del archivo donde se guarda la malla rotada
 *@param[in] n_nodes Es el numero de nodos de la malla
 *@param[in] n_elements Es el numero de elementos de la malla
 *@param[in] nodes Es el arreglo que contiene la informacion de los nodos
 *@param[in] elements Es el arreglo que contiene las conectividades de la malla
 *@param[in] material Es el arreglo que contiene los materiales de los elementos
 */
void SaveRotatedMesh( char* name , size_t n_nodes , size_t n_elements , double** nodes , 
											size_t** elements , int* material ){
	FILE* fp = NULL;
	fp = fopen( name , "w" );
	fprintf( fp , "; Geometry file\n\n" );	
	fprintf( fp , "{Nodes}\n" );
	fprintf( fp , "3 ; Dimension\n" );
	fprintf( fp , "%ld ; Nodes count\n" , n_nodes );
	fprintf( fp , "; X1 X2 ...\n" );
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		fprintf( fp , "%lf %lf %lf\n" , nodes[ i_node ][ 0 ] , nodes[ i_node ][ 1 ] , nodes[ i_node ][ 2 ]  );
	}
	fprintf( fp , "\n" );
	fprintf( fp , "{Mesh}\n" );
	fprintf( fp , "5 ; Element type (2=Triangle, 3=Quadrilateral, 4=Tetrahedra, 5=Hexahedra)\n" );
	fprintf( fp , "8 ; Nodes per element\n" );
	fprintf( fp , "%ld ; Elements count\n" , n_elements );
	fprintf( fp , "; Material Node1 Node2 ...\n" );
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		fprintf( fp , "%d %ld %ld %ld %ld %ld %ld %ld %ld\n" , material[ i_elem ] , elements[ i_elem ][ 0 ] , elements[ i_elem ][ 1 ] , elements[ i_elem ][ 2 ] , elements[ i_elem ][ 3 ] , elements[ i_elem ][ 4 ] , elements[ i_elem ][ 5 ] , elements[ i_elem ][ 6 ] , elements[ i_elem ][ 7 ]  );
	}

	fclose( fp );
}

/**
 *Checking if number of arguments is correct or not
 *@param[in] argc Number of arguments in program
 */
void CheckParameters( int argc ){
	if(  argc != 3  ){
		std::cout << "\n\n";
		std::cout << "**************************************************************************************" << std::endl;
		std::cout << "********************************MESH ROTATE APPLICATION*******************************" << std::endl;
		std::cout << "Received wrong number of arguments " << std::endl; 
		std::cout << "Arguments must be: " << std::endl;
		std::cout << "  --[1][string][input ] Name of file containing geometry informaiton in FEM format" << std::endl;
		std::cout << "  --[2][string][output] Name of file containing geometry informaiton in FEM format" << std::endl;
		std::cout << "********************************MESH ROTATE APPLICATION*******************************" << std::endl;
		std::cout << "**************************************************************************************" << std::endl;
		std::cout << "\n\n";
		assert( 0 );
	}
}

int main( int argc , char** argv ){

	CheckParameters( argc );

	size_t n_nodes , n_elements;
	double **nodes;
	size_t **elements;
	int* material;

	ReadNNodesAndNElements(  argv[ 1 ] , &n_nodes , &n_elements  );
	nodes = new double*[ n_nodes ];
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++ ){
		nodes[ i_node ] = new double[ 3 ];
	}
	elements = new size_t*[ n_elements ];
	material = new int[ n_elements ];
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ; i_elem++  ){
		elements[ i_elem ] = new size_t[ 8 ];
	}

	ReadNodesAndElements( argv[ 1 ] , n_nodes , n_elements , nodes , elements , material );
	RotateNodes( n_nodes , nodes );
	SaveRotatedMesh( argv[ 2 ] , n_nodes , n_elements , nodes , elements , material );


	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		delete[] nodes[ i_node ];
	}
	delete[] nodes;
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ; i_elem++  ){
		delete[] elements[ i_elem ];
	}
	delete[] elements;
	delete[] material;
  return 0;
}
