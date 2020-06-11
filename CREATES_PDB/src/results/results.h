#pragma once

//C system files
#include <math.h>
#include <assert.h>
#include <string.h>

//C++ system files
#include <vector>
#include <iostream>
#include <fstream>


//Pther includes
#include "../geometry/geometry.h"
/**
 *This class contains the displacements for a node on a mesh
 *@param flag_ Is the flag that indicates if the displacements were asigned
 *@param magnitude_ Is the magnitude of the displacement vector
 *@param vector_ Is the list of displacements in each direction
 */
class Displacement{
		bool flag_;
		double magnitude_;
		double vector_[ 3 ];

	public:

		//CONSTRUCTOR AND DESTRUCTOR
		Displacement();
		Displacement( double* values );
		Displacement( double X, double Y, double Z );
		~Displacement();

		//GETS
		double GetDisplacement( int position );
		double GetNormalizedDisplacement( int position );
		double GetDisplacementMagnitude( );

		//SETS
		void SetDisplacement( int position , double value );
		void SetDisplacement( double X, double Y, double Z );
		void SetDisplacement( double* values );

		//UTILITIES
		void CalculateMagnitude(  );

		//MODIFIERS
 		void NormalizeVector(  );
		
		//DEBUG
		void PrintInfo(  );
		bool IsNormalized(  );
		bool IsSet();

};

/**
 *This class contains the strains evaluated on the gauss points from an element
 *@param n_gauss_point_ Number of gauss points in the element 
 *@param set_flag_ flag indicating if the strain was set
 *@param average_flag_ Flag indicating if the average for the satrain was calculated
 *@param averages_ Is the average of each kind of result
 *@param values_ Is a matrix where the rows contain the values in the following order
 **** 		Sxx		Syy		Szz		Sxy		Syz		Sxz				****
 */
class Strain{
    int n_gauss_points_;
		bool set_flag_[ 6 ];
		bool average_flag_[ 6 ];
		double averages_[ 6 ];
		double* values_[ 6 ];

	public:

		//CONSTRUCTOR AND DESTRUCTOR
		Strain(  );
		Strain( double** values , int n_gauss_points );
		Strain( double* Sxx , double*  Syy , double* Szz , double*  Sxy , double* Syz , double* Sxz , int n_gauss_points );
		~Strain(  );

		//GETS
		void GetValue( int kind_requested ,double* res_vector );
		double GetValue( int kind_requested ,int gauss_index );
		double GetAverage( int kind_requested );
		void GetAverage( double* averages );
    int GetNGaussPoints(  );

		//SETS
		void SetValue( int kind , double* values );
		void SetValue( int kind , int gauss_index , double value );
		void SetValue( double** values , int n_gauss_points );
    void SetNGaussPoints( int n_gauss_points );

		//UTILITIES
		void CalculateAverages( );
		void ActualizeAverage( int kind );
		void ActualizeAverage(  );

		//MODIFIERS

		//CLEANERS
		void ClearMemory(  );

		//DEBUG
		void PrintInfo(  );
		bool IsSet( int kind );
		bool AverageIsCalculated( int kind  );
};


/**
 *This class contains the stress evaluated on the gauss points from an element
 *@param n_gauss_point_ Number of gauss points in the element 
 *@param set_flag_ flag indicating if the stress was set
 *@param average_flag_ Flag indicating if the average for the stress was calculated
 *@param averages_ Is the average of each kind of result
 *@param values_ Is a matrix where the columns contain the values in the following order
 **** 		Sxx		Syy		Szz		Sxy		Syz		Sxz				****
 */
class Stress{
    int n_gauss_points_;
		bool set_flag_[ 6 ];
		bool average_flag_[ 6 ];
		double averages_[ 6 ];
		double* values_[ 6 ];

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Stress(  );
		Stress( double** values , int n_gauss_points  );
		Stress( double* Sxx , double*  Syy , double* Szz , double*  Sxy , double* Syz , double* Sxz , int n_gauss_points  );
		~Stress(  );

		//GETS
		void GetValue( int kind_requested ,double* res_vector );
		double GetValue( int kind_requested ,int gauss_index );
		double GetAverage( int kind_requested );
		void GetAverage( double* averages );
    int GetNGaussPoints(  );

		//SETS
		void SetValue( int kind , double* values );
		void SetValue( int kind , int gauss_index , double value );
		void SetValue( double** values , int n_gauss_points );
    void SetNGaussPoints( int n_gauss_points );

		//UTILITIES
		void CalculateAverages( );
		void ActualizeAverage( int kind );
		void ActualizeAverage(  );

		//MODIFIERS

		//CLEANERS
		void ClearMemory(  );

		//DEBUG
		void PrintInfo(  );
		bool IsSet( int kind );
		bool AverageIsCalculated( int kind );

};

/**
 *This class contains the vonmisses evaluated on the gauss points from an element
 *@param n_gauss_point_ Number of gauss points in the element 
 *@param set_flag_ flag indicating if the vonmises was set
 *@param average_flag_ Flag indicating if the average for the vonmises was calculated
 *@param average_ Is the average of the vonmises
 *@param values_ Is the vector of the vonmises in each gauss point from the element
 */
class Vonmises{
    int n_gauss_points_;
		bool set_flag_;
		bool average_flag_;		
		double average_;
		double* values_;

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Vonmises(  );
		Vonmises( double* values , int n_gauss_points );
		~Vonmises(  );

		//GETS
		double GetValue( int gauss_index );
		void GetValue( double* values );
		double GetAverage( );
    int GetNGaussPoints(  );

		//SETS
		void SetValue( int gauss_index , double values );
		void SetValue( double* values , int n_gauss_points );
    void SetNGaussPoints( int n_gauss_points );

		//UTILITIES
		void CalculateAverage( );
		void ActualizeAverage(  );

		//MODIFIERS

		//CLEANERS
		void ClearMemory(  );

		//DEBUG
		void PrintInfo(  );
		bool IsSet(  );
		bool AverageIsCalculated(  );

};

/**
 *This class contains the solution from a capside mesh, it contains the solution of 
 *displacements, strain, stress and vonmises.
 *@param name_ Is the name of the file containing the results
 *@param ptr_ Is the pointer to the file that contains the results
 *@param n_nodes_ Is the number of nodes of the mesh
 *@param n_elements_ Is the number of elements in the mesh
 *@param displacements_ Is the list of displacements on the mesh nodes
 *@param strains_ Is the list of strains on the gaus points of each element
 *@param stress_ Is the list of stresses on the gauss points of each element
 *@param vonmises_ Is the list of vonmisses on the gauss points
 */
class Results{
		char* name_;
		FILE* ptr_;
		size_t n_nodes_;
		size_t n_elements_;
		std::vector<Displacement*> displacements_;
		std::vector<Strain*> strains_;
		std::vector<Stress*> stress_;
		std::vector<Vonmises*> vonmises_;

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Results(  );
		Results( char* name , size_t n_nodes , size_t n_elements );
		~Results(  );

		//GETS
		size_t GetNNodes(  );
		size_t GetNElements(  );
		double GetDisplacementMagnitude( size_t index );
		double GetDisplacement( size_t index , int position );
		void GetDisplacement( size_t index , double* values );
		Displacement* GetDisplacement( size_t index );
		void GetStrain( size_t index , int kind_requested ,double* res_vector );
		double GetStrain( size_t index , int kind_requested ,int gauss_index );
		double GetAverageStrain( size_t index , int kind_requested );
		Strain* GetStrain( size_t index );
		void GetStress( size_t index , int kind_requested ,double* res_vector );
		double GetStress( size_t index , int kind_requested ,int gauss_index );
		double GetAverageStress( size_t index , int kind_requested );
		Stress* GetStress( size_t index );
		double GetVonmises( size_t index , int gauss_index );
		void GetVonmises( size_t index , double* values );
		double GetAverageVonmises( size_t index );
		Vonmises* GetVonmises( size_t index );
		char* GetName();

		//SETS
		void SetDisplacement( size_t index , int position , double value );
		void SetDisplacement( size_t index , double X, double Y, double Z );
		void SetDisplacement( size_t index , double* values );
		void SetStrain( size_t index , int kind , double* values ,int n_gauss_points );
		void SetStrain( size_t index , int kind , int gauss_index , double value );
		void SetStrain( size_t index , double** values , int n_gauss_points );
		void SetStress( size_t index , int kind , double* values , int n_gauss_points);
		void SetStress( size_t index , int kind , int gauss_index , double value );
		void SetStress( size_t index , double** values , int n_gauss_points );
		void SetVonmises( size_t index , int gauss_index , double values );
		void SetVonmises( size_t index , double* values , int n_gauss_points );

		//UTILITIES
		void AddDisplacement( Displacement* disp );
		void AddStrain( Strain* strain );
		void AddStress( Stress* stress );
		void AddVonmises( Vonmises* von );


		//MODIFIERS

		//DEBUG
		void PrintInfo();	
};






