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
#include "../file/file.h"
#include "../boundary/boundary.h"
#include "../mesh/mesh.h"

/**
 *Atom class
 *This class contains the information of one atom in the capside
 *@param coords_ Are the spacial coordinates of the atom
 *@param keys_ Is the binary code for the atom on the octree
 *@param kind_ Is the kind of atom (C,N,CB,etc)
 *@param radius_ Is the radius of the atom used
 */
class Atom{
	
	double coords_[ 3 ];
	size_t keys_[ 3 ];
	char kind_[ 10 ];
	double radius_;

	public:

		//CONSTRUCTOR AND DESTRUCTOR
		Atom(  );
		Atom( char* kind , double radius , double* coords );
		~Atom(  );

		//GETS
		bool GetKind( char kind[10] );
		double GetCoordinate( int position );
		void GetCoordinates( double* coords );
		double GetRadius();
		int GetMaterialIndex(  );
		size_t GetKey( int position );

		//SETS
		bool SetKind( char* kind );
		bool SetCoordinate( double coord , int position );
		bool SetCoordinates( double* coords );
		void SetRadius( double radius );

		//UTILITIES

		//CLEANERS 

		//DEBUG
		void PrintParameters();
		void PrintInfo();


};

/**
 *Class aminoacid
 *This class contains the information referred to one aminoacid in the protein
 *@param kind_ Is the kind of the aminoacid, this must be a value between 0-19 indicating 
 *             a index in the variable valids_.
 *@param atoms_ Is the set of atoms that belongs to this aminoacid
 *@param valids_ Is the variable containig the list of valid names for aminoacids
 *@param radius_ Is the list of radius for the atoms
 */
class Aminoacid{

	int kind_;
	std::vector<Atom*> atoms_;
	char* valids_[ 20 ] = {  {"ALA"}/*Alanine*/ ,       {"ARG"}/*Arginine*/ ,      {"ASN"}/*Asparagine*/ , 
															 {"ASP"}/*Aspartic acid*/ , {"CYS"}/*Cysteine*/ ,      {"GLU"}/*Glutamic acid*/ , 
															 {"GLN"}/*Glutamine*/ ,     {"GLY"}/*Glycine*/ ,       {"HIS"}/*Histidine*/ , 
															 {"ILE"}/*Isoleucine*/ ,    {"LEU"}/*Leucine*/ ,       {"LYS"}/*Lysine*/ , 
															 {"MET"}/*Methionine*/ ,    {"PHE"}/*Phenylalanine*/ , {"PRO"}/*Proline*/ , 
															 {"SER"}/*Serine*/ ,        {"THR"}/*Threonine*/ ,     {"TRP"}/*Tryptophan*/ , 
															 {"TYR"}/*Tyrosine*/ ,      {"VAL"}/*Valine*/ };
	double radius_[ 6 ] = { 1.78/*Carbon*/ , 1.6/*Hydrogen*/ , 1.6/*Oxygen*/ , 1.8/*Nitrogen*/ ,1.83/*Phosphorus*/ , 1.2/*Sulfur*/ };

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Aminoacid(  );
		Aminoacid( int kind );
		~Aminoacid(  );

		//GETS
		size_t GetNAtoms(  );
		int GetKind(  );
		double GetAtomRadius( char* kind );	
		double GetAtomCoordinate( size_t index , int position );	
		void GetAtomCoordinates( size_t index , double* coords );
		Atom* GetAtom( size_t index );

		//SETS
		bool SetKind( int kind );
		bool SetAtomCoordinate( size_t index , double coord , int position );
		bool SetAtomCoordinates( size_t index , double* coords );

		//UTILITIES
		bool AddAtom( char* kind , double diameter , double* coords );
		bool IsValid( char* kind , int* position );
		bool AreEqual( int index );

		//CLEANERS
		bool EmptyAtoms(  );

		//DEBUG
		void PrintParameters();
		void PrintInfo();
};

/**
 *Class protein
 *This class contains the name of the protein and the set of aminoacids that are part of 
 *the protein
 *@param kind_ Is the kind of protein, this is a letter from A-Z.
 *@param aminoacids_ Is the vector of aminoacids composing the protein
 */
class Protein{

	char kind_;
	std::vector<Aminoacid*> aminoacids_;
	char* valids_[ 20 ] = {  {"ALA"}/*Alanine*/ ,       {"ARG"}/*Arginine*/ ,      {"ASN"}/*Asparagine*/ , 
															 {"ASP"}/*Aspartic acid*/ , {"CYS"}/*Cysteine*/ ,      {"GLU"}/*Glutamic acid*/ , 
															 {"GLN"}/*Glutamine*/ ,     {"GLY"}/*Glycine*/ ,       {"HIS"}/*Histidine*/ , 
															 {"ILE"}/*Isoleucine*/ ,    {"LEU"}/*Leucine*/ ,       {"LYS"}/*Lysine*/ , 
															 {"MET"}/*Methionine*/ ,    {"PHE"}/*Phenylalanine*/ , {"PRO"}/*Proline*/ , 
															 {"SER"}/*Serine*/ ,        {"THR"}/*Threonine*/ ,     {"TRP"}/*Tryptophan*/ , 
															 {"TYR"}/*Tyrosine*/ ,      {"VAL"}/*Valine*/ };

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Protein(  );
		Protein( char kind );
		~Protein(  );

		//GETS
		size_t GetNAminoacids(  );
		char GetKind(  );	
		Aminoacid* GetAminoacid( size_t index );

		//SETS
		bool SetKind( char kind );

		//UTILITIES
		bool IsValid( char* kind , int* posiiton );
		bool AddAminoacid( char* kind );

		//CLEANERS
		bool EmptyAminoacids(  );

		//DEBUG
		void PrintParameters();
		void PrintInfo();
};

/**
 *Vdb Is a class to manage the reading of a vdb file containing the capside information
 *@param name_ Is the name of the file that contains the input file
 *@param ptr_ Is the pointer to the file that contains the 
 *@param mesh_ Is the mesh of the atoms
 *@param proteins_ Is the list of proteins forming the capside
 *@param valids_ Is the list of valid aminoacids
 */
class Vdb{
	char* name_;
	FILE* ptr_;
	Spheres* mesh_;
	std::vector<Protein*> proteins_;
	char* valids_[ 20 ] = {  {"ALA"}/*Alanine*/ ,       {"ARG"}/*Arginine*/ ,      {"ASN"}/*Asparagine*/ , 
															 {"ASP"}/*Aspartic acid*/ , {"CYS"}/*Cysteine*/ ,      {"GLU"}/*Glutamic acid*/ , 
															 {"GLN"}/*Glutamine*/ ,     {"GLY"}/*Glycine*/ ,       {"HIS"}/*Histidine*/ , 
															 {"ILE"}/*Isoleucine*/ ,    {"LEU"}/*Leucine*/ ,       {"LYS"}/*Lysine*/ , 
															 {"MET"}/*Methionine*/ ,    {"PHE"}/*Phenylalanine*/ , {"PRO"}/*Proline*/ , 
															 {"SER"}/*Serine*/ ,        {"THR"}/*Threonine*/ ,     {"TRP"}/*Tryptophan*/ , 
															 {"TYR"}/*Tyrosine*/ ,      {"VAL"}/*Valine*/ };
	public:

		//CONSTRUCTOR AND DESTRUCTOR
		Vdb(  );
		Vdb( char* name );
		~Vdb(  );

		//GETS
		size_t GetNProteins(  );
		Protein* GetProtein( size_t index );
		FILE* GetPointer();
		char* GetName(); 

		//SETS
		void SetAtomOnCapside( char params[11][500] );

		//UTILITIES
		void AddAtom( char params[ 11 ][ 500 ] , int aminocid_index );
		bool AddProtein( char kind );
		void ReadCompleteFile(  );
		int ReadLine( );
		bool AminoacidIsValid( char* kind , int* posiiton );
		void FillSphericalMeshInformation(  );
		void SaveGiDSphericalMesh(  );
    void ScaleSpheresMesh(  ); 
		double UnscaleCoordinate( int i_pos , double coord );
		void PrintResumenOnWarnings( );
		double ScaleCoord( int i_pos , double coord );

		//CLEANERS
		bool EmptyProteins(  );		

		//DEBUG
		void PrintParameters();
		void PrintInfo();
		
};




