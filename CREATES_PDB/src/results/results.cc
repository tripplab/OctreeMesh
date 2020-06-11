#include "results.h"

//////////////////////////////////////////////////////////////////////////////////////////
/////														DISPLACEMENT METHODS																		//
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Displacement::Displacement(){
  flag_ = false;
  magnitude_ = -1.0;
  for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
    vector_[ i_dim ] = 0.0;
}

/**
 *Constructor receiving the displacements
 *This method set the displacements in each dimension and calculates the magnitude of the 
 *vector in order to have the value if the normalized vector is required
 *@param[in] values Is the vector containing the displacements
 */
Displacement::Displacement( double* values ){
  flag_ = true;
  for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
    vector_[ i_dim ] = values[ i_dim ];  
  magnitude_ = Magnitude( vector_ );
}

/**
 *Constructor receiving the displacements
 *This method set the displacement vector and calculates iis magnitude in order to have 
 *this value when this is requested
 *@param[in] X Is the displacement in the X direction
 *@param[in] Y Is the displacement in the Y direction
 *@param[in] Z is the displacement in the Z direction
 */ 
Displacement::Displacement( double X, double Y, double Z ){
  flag_ = true;
  vector_[ 0 ] = X;
  vector_[ 1 ] = Y;
  vector_[ 2 ] = Z;
  magnitude_ = Magnitude( vector_ );
}

/**
 *Default destructor
 */
Displacement::~Displacement(){
  flag_ = false;
  magnitude_ = -1.0;
  for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
    vector_[ i_dim ] = 0.0;  
  
}

//GETS
/**
 *Getting displacement over an axis
 *@param[in] position Is the index that indicates the direction of the displacement requested
 *@return A double value with the displacement requested 
 */
double Displacement::GetDisplacement( int position ){
  assert( ( position >=0 ) && ( position < 3 )  );
  return vector_[ position ];
}

/**
 *Gettijng normalized displacement over an axis
 *@param[in] position Is the index indicating the direction of the dispalcement requested
 *@return A double value with the normalized displacement
 */
double Displacement::GetNormalizedDisplacement( int position ){
  assert( ( position >=0 ) && ( position < 3 )  );
  return ( vector_[ position ] / magnitude_ );
}

/**
 *Getting the magnitude of the displacement
 *@return A double value with the magnitude of the displacement
 */
double Displacement::GetDisplacementMagnitude( ){
  return magnitude_;
}

//SETS
/**
 *Setting displacement over a direction
 *@param[in] position Index that indicates the direction of the displacement to be set
  *param[in] value Is the displacement value to be set
 */
void Displacement::SetDisplacement( int position , double value ){
  assert( ( position >=0 ) && ( position < 3 )  );
  vector_[ position ] = value; 
}

/**
 *Setting displacements over all directions
 *@param[in] X Is the displacement over the X direction
 *@param[in] Y Is the displacement over the Y direction
 *@param[in] Z Is the displacement over the Z direction
 */
void Displacement::SetDisplacement( double X, double Y, double Z ){
  vector_[ 0 ] = X;
  vector_[ 1 ] = Y;
  vector_[ 2 ] = Z;
}

/**
 *Setting displacement vector
 *@param[in] values Is the vector that contains the displacements over all directions
 */
void Displacement::SetDisplacement( double* values ){
  for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
    vector_[ i_dim ] = values[ i_dim ];
  }
}

//UTILITIES
/**
 *Calculate the magnitude of the displacement
 */
void Displacement::CalculateMagnitude(  ){
  magnitude_ = Magnitude( vector_ );
}

//MODIFIERS
/**
 *Normalizing the displacement vector
 */
void Displacement::NormalizeVector(  ){
  for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  )
    vector_[ i_dim ] /= magnitude_; 
}

//DEBUG
/**
 *Printing displacements information
 */
void Displacement::PrintInfo(  ){
  std::cout << "Disp: " << vector_[ 0 ] << " " << vector_[ 1 ] << " " << vector_[ 2 ] << 
  "|| Disp ||: " << magnitude_ << std::endl; 
}

/**
 *Testing if vector is normalized
 */
bool Displacement::IsNormalized(  ){
  return ( fabs( Magnitude( vector_ ) - 1.0 ) < 1e-10 );
}

/**
 *Testing if the diaplcement vector is set
 */
bool Displacement::IsSet(){
  return flag_;
}


//////////////////////////////////////////////////////////////////////////////////////////
/////																	STRAIN METHODS																		//
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */ 
Strain::Strain(  ){
  n_gauss_points_ = 0;
  for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
    set_flag_[     i_kind ] = false;
    average_flag_[ i_kind ] = false;
    averages_[ i_kind ] = 0.0;
    values_[ i_kind ] = NULL;
  }
}

/**
 *Constructor receiving the matrix of strain and the number of gauss points
 *@param[in] values is the matrix that contains the strain value
 *@param[in] n_gauss_points Value indicating the number of gauss point on the element
 */ 
Strain::Strain( double** values , int n_gauss_points ){
  n_gauss_points_ = n_gauss_points;
  for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
    set_flag_[ i_kind ] = true;
    averages_[ i_kind ] = 0.0;
    values_[ i_kind ] = new double[ n_gauss_points_ ];
    for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
      values_[ i_kind ][ i_point ] = values[ i_kind ][ i_point ];
      averages_[ i_kind ] += values_[ i_kind ][ i_point ];
    }
    averages_[ i_kind ] /= (double)n_gauss_points_;
    average_flag_[ i_kind ] = true;
  }
}

/**
 *Constructor receiving the strains as vectors
 *@param[in] Sxx Is the strain value in the xx component
 *@param[in] Syy Is the strain value in the yy component
 *@param[in] Szz Is the strain value in the zz component
 *@param[in] Sxy Is the strain value in the xy component
 *@param[in] Syz Is the strain value in the yz component
 *@param[in] Sxz Is the strain value in the xz component
 *@param[in] n_gauss_points Is the number of gauss points on the element
 */
Strain::Strain( double* Sxx , double*  Syy , double* Szz , double*  Sxy , double* Syz , double* Sxz , int n_gauss_points ){
  n_gauss_points_ = n_gauss_points;
  for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
    set_flag_[ i_kind ] = true;
		averages_[ i_kind ] = 0.0;
		values_[ i_kind ] = new double[ n_gauss_points_ ];
  }
  for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		//Sxx
		values_[ 0 ][ i_point ] = Sxx[ i_point ];	
		averages_[ 0 ] += values_[ 0 ][ i_point ];
		//Syy
		values_[ 1 ][ i_point ] = Syy[ i_point ];	
		averages_[ 1 ] += values_[ 1 ][ i_point ];
		//Szz
		values_[ 2 ][ i_point ] = Szz[ i_point ];	
		averages_[ 2 ] += values_[ 2 ][ i_point ];
		//Sxy
		values_[ 3 ][ i_point ] = Sxy[ i_point ];	
		averages_[ 3 ] += values_[ 3 ][ i_point ];
		//Syz
		values_[ 4 ][ i_point ] = Syz[ i_point ];	
		averages_[ 4 ] += values_[ 4 ][ i_point ];
		//Sxz
		values_[ 5 ][ i_point ] = Sxz[ i_point ];	
		averages_[ 5 ] += values_[ 5 ][ i_point ];
  }
  for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		average_flag_[ i_kind ] = true;
		averages_[ i_kind ] /= (double)n_gauss_points_;
	}
}

/**
 *Default destructor
 */
Strain::~Strain(  ){
  n_gauss_points_ = 0;
  for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
    set_flag_[     i_kind ] = false;
    average_flag_[ i_kind ] = false;
    averages_[ i_kind ] = 0.0;
		delete[] values_[ i_kind ];
    values_[ i_kind ] = NULL;
  }
}

//GETS
/**
 *Getting a kind of strain on all gauss points
 *@param[in] kind_requested Index of the strain requested as showed in next lines:
	0-Sxx						1-Syy						2-Szz						3-Sxy						4-Syz						5-Sxz
 *@param[out] res_vector Is the array where the strain will be stored
 */
void Strain::GetValue( int kind_requested ,double* res_vector ){
	assert( ( kind_requested >= 0 ) && ( kind_requested < 6 ) );
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		res_vector[ i_point ] = values_[ kind_requested ][ i_point ];
	}		
}

/**
 *Getting strain value on a gauss point
 *@param[in] kind_requested Index of the strain requested as showed in next lines:
	0-Sxx						1-Syy						2-Szz						3-Sxy						4-Syz						5-Sxz
 *@param[in] gauss_index Is the index of the gaus index value requested
 *@return A double value with the strain requested
 */
double Strain::GetValue( int kind_requested ,int gauss_index ){
	assert( ( kind_requested >= 0 ) && ( kind_requested < 6 ) );
	assert( ( gauss_index >=0 ) && ( gauss_index <n_gauss_points_ ) );
	return values_[ kind_requested ][ gauss_index ];
}

/**
 *Getting average frome one kind of strain
 *@param[in] kind_requested Index of the strain requested as showed in next lines:
	0-Sxx						1-Syy						2-Szz						3-Sxy						4-Syz						5-Sxz
 *@return A double value with the average of the strain
 */
double Strain::GetAverage( int kind_requested ){
	assert( ( kind_requested >= 0 ) && ( kind_requested < 6 ) );
	return averages_[ kind_requested ];
}

/**
 *Getting all averages of strain
 *@param[out] averages Is the array were the averages will be stored
 */
void Strain::GetAverage( double* averages ){
	for( int i_kind = 0  ;  i_kind < 6  ;  i_kind++ ){
		averages[ i_kind ] = averages_[ i_kind ];	
	}
}

/**
 *Getting number of gauss points on the element
 *@return A int value with the number of gauyss points
 */
int Strain::GetNGaussPoints(  ){
	return n_gauss_points_;
}

//SETS
/**
 *Setting a kind of strain in the matrix
 *@param[in] kind Is the kind of strain to be set
	0-Sxx						1-Syy						2-Szz						3-Sxy						4-Syz						5-Sxz
 *@param[in] values Is the array containing the values to be set
 */
void Strain::SetValue( int kind , double* values ){
	assert( ( kind >= 0 ) && ( kind < 6 ) );
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		values_[ kind ][ i_point ] =  values[ i_point ];
	}
}

/**
 *Setting a value of strain in the matrix
 *@param[in] kind Is the kind of strain to be set
	0-Sxx						1-Syy						2-Szz						3-Sxy						4-Syz						5-Sxz
 *@param[in] gauss_index Is the index of the gauss point where the value will be set
 *@param[in] value Is the value to be set
 */
void Strain::SetValue( int kind , int gauss_index , double value ){
	assert( ( kind >= 0 ) && ( kind < 6 ) );
	assert( ( gauss_index >= 0 ) && ( gauss_index < n_gauss_points_ ) );
	values_[ kind ][ gauss_index ] = value;
}

/**
 *Setting matrix of values of strain
 *@param[in] values Is the array containing the values of strain
 *@param[in] n_gauss_points
 */
void Strain::SetValue( double** values , int n_gauss_points ){
	n_gauss_points_ = n_gauss_points;
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
			values_[ i_kind ][ i_point ] = values[ i_kind ][ i_point ];
		}		
	}
}

/**
 *Setting number of gauss points
 *@param[in] n_gauss_points Is the number of gauss points to be set
 */
void Strain::SetNGaussPoints( int n_gauss_points ){
	n_gauss_points_ = n_gauss_points;
}

//UTILITIES
/**
 *Calculate averages of strain by kind
 */
void Strain::CalculateAverages( ){
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		averages_[ i_kind ] = 0.0;
		for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
			averages_[ i_kind ] += values_[ i_kind ][ i_point ];
		}
		averages_[ i_kind ] /= (double)n_gauss_points_;
		average_flag_[ i_kind ] =  true;
	}
}

/**
 *Actualize average from a variable
 *@param[in] kind Is the kind of strain average to actualize
 */
void Strain::ActualizeAverage( int kind ){
	assert( ( kind >= 0 ) && ( kind < 6 ) );
	averages_[ kind ] = 0.0;
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		averages_[ kind ] += values_[ kind ][ i_point ];
	}
	averages_[ kind ] /= (double)n_gauss_points_;
}

/**
 *Actualizing the average from all kind of strains in the class
 */
void Strain::ActualizeAverage(  ){
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		averages_[ i_kind ] = 0.0;
		for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
			averages_[ i_kind ] += values_[ i_kind ][ i_point ];
		}
		averages_[ i_kind ] /= (double)n_gauss_points_;
	}
}

//MODIFIERS

//CLEANERS
/**
 *Cleaning the memory allocated before in case that exists
 */
void Strain::ClearMemory(  ){
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		if(  values_[ i_kind ]  ){
			delete[] values_[ i_kind ];
		}
	} 
}

//DEBUG
/**
 *Printing information contained in the class
 */
void Strain::PrintInfo(  ){
	std::cout << "Number of gauss points: " << n_gauss_points_ << std::endl;
	std::cout << "Flags:      set      average: " << std::endl;
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		std::cout << set_flag_[ i_kind ] << "        " << average_flag_[ i_kind ] << std::endl;	
	}
	std::cout << "Values:    Sxx    Syy    Szz    Sxy    Syz    Sxz  :  Average" << std::endl;
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
			std::cout << values_[ i_kind ][ i_point ] << "    ";
		}
		std::cout << averages_[ i_kind ] << std::endl;
	}
}

/**
 *Testing if values for a kind of strain is set
 *@param[in] kind Is the kind of strain average to actualize
 *@return A bool value indicating if values were set or not
 */
bool Strain::IsSet( int kind ){
	assert( ( kind >= 0 ) && ( kind < 6 ) );
	return set_flag_[ kind ];
}

/**
 *Testing if average is calculated on a kind of strain
 *@param[in] kind Is the kind of strain average to actualize
 *@return A bool value indicating if values were set or not
 */
bool Strain::AverageIsCalculated( int kind  ){
	assert( ( kind >= 0 ) && ( kind < 6 ) );
	return average_flag_[ kind ];
}


//////////////////////////////////////////////////////////////////////////////////////////
/////																	STRESS METHODS																		//
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Stress::Stress(  ){
  n_gauss_points_ = 0;
  for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
    set_flag_[     i_kind ] = false;
    average_flag_[ i_kind ] = false;
    averages_[ i_kind ] = 0.0;
    values_[ i_kind ] = NULL;
  }
}

/**
 *Constructor receiving the stress matrix and the number of gauss points
 *@param[in] values Is the matrix of the stresses
 *@param[in] n_gauss_points Is the number of gauss points used in the result
 */
Stress::Stress( double** values , int n_gauss_points  ){
  n_gauss_points_ = n_gauss_points;
  for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
    set_flag_[ i_kind ] = true;
    averages_[ i_kind ] = 0.0;
    values_[ i_kind ] = new double[ n_gauss_points_ ];
    for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
      values_[ i_kind ][ i_point ] = values[ i_kind ][ i_point ];
      averages_[ i_kind ] += values_[ i_kind ][ i_point ];
    }
    averages_[ i_kind ] /= (double)n_gauss_points_;
    average_flag_[ i_kind ] = true;
  }
}

/**
 *Construtor receiving arrays of the stresses
 *@param[in] Sxx Is the strain value in the xx component
 *@param[in] Syy Is the strain value in the yy component
 *@param[in] Szz Is the strain value in the zz component
 *@param[in] Sxy Is the strain value in the xy component
 *@param[in] Syz Is the strain value in the yz component
 *@param[in] Sxz Is the strain value in the xz component
 *@param[in] n_gauss_points Is the number of gauss points on the element
 */
Stress::Stress( double* Sxx , double*  Syy , double* Szz , double*  Sxy , double* Syz , double* Sxz , int n_gauss_points  ){
  n_gauss_points_ = n_gauss_points;
  for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
    set_flag_[ i_kind ] = true;
		averages_[ i_kind ] = 0.0;
		values_[ i_kind ] = new double[ n_gauss_points_ ];
  }
  for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		//Sxx
		values_[ 0 ][ i_point ] = Sxx[ i_point ];	
		averages_[ 0 ] += values_[ 0 ][ i_point ];
		//Syy
		values_[ 1 ][ i_point ] = Syy[ i_point ];	
		averages_[ 1 ] += values_[ 1 ][ i_point ];
		//Szz
		values_[ 2 ][ i_point ] = Szz[ i_point ];	
		averages_[ 2 ] += values_[ 2 ][ i_point ];
		//Sxy
		values_[ 3 ][ i_point ] = Sxy[ i_point ];	
		averages_[ 3 ] += values_[ 3 ][ i_point ];
		//Syz
		values_[ 4 ][ i_point ] = Syz[ i_point ];	
		averages_[ 4 ] += values_[ 4 ][ i_point ];
		//Sxz
		values_[ 5 ][ i_point ] = Sxz[ i_point ];	
		averages_[ 5 ] += values_[ 5 ][ i_point ];
  }
  for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		average_flag_[ i_kind ] = true;
		averages_[ i_kind ] /= (double)n_gauss_points_;
	}
}

/**
 *Default destructor
 */
Stress::~Stress(  ){
  n_gauss_points_ = 0;
  for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
    set_flag_[     i_kind ] = false;
    average_flag_[ i_kind ] = false;
    averages_[ i_kind ] = 0.0;
		delete[] values_[ i_kind ];
    values_[ i_kind ] = NULL;
  }
}

//GETS
/**
 *Getting values of a kind of stress
 *@param[in] kind_requested Index of the strain requested as showed in next lines:
	0-Sxx						1-Syy						2-Szz						3-Sxy						4-Syz						5-Sxz
 *@param[out] res_vector Is the array where the requested values will be stored 
 */
void Stress::GetValue( int kind_requested ,double* res_vector ){
	assert( ( kind_requested >=0 ) && ( kind_requested < 6 ) );
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		res_vector[ i_point ] =  values_[ kind_requested ][ i_point ];
	}
}

/**
 *Getting value from stress 
 *@param[in] kind_requested Is the kind of stress requested
 *@param[in] gauss_index Is the gauss index of the requested value
 *@return A double value with the value requested
 */
double Stress::GetValue( int kind_requested ,int gauss_index ){
	assert( ( kind_requested >=0 ) && ( kind_requested < 6 ) );
	assert( ( gauss_index >=0 ) && ( gauss_index < n_gauss_points_ ) );
	return values_[ kind_requested ][ gauss_index ];
}

/**
 *Getting average from a kind of stress
 *@param[in] kind_requested Is the index of the average stress requested
 *@return A double value with the average
 */
double Stress::GetAverage( int kind_requested ){
	assert( ( kind_requested >=0 ) && ( kind_requested < 6 ) );
	return averages_[ kind_requested ];
}

/**
 *Getting averages of stresses
 *@param[out] averages Is the array where the averages will be stored
 */
void Stress::GetAverage( double* averages ){
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		averages[ i_kind ] = averages_[ i_kind ];
	}
}

/**
 *Getting number of gauss points 
 *@return A int value with the number of gauss points in the element
 */
int Stress::GetNGaussPoints(  ){
	return n_gauss_points_;
}

//SETS
/**
 *Setting values of stresses
 *@param[in] kind Is the kind of stress to be set
 *@param[in] values Is the array containing the values to be set
 */
void Stress::SetValue( int kind , double* values ){
	assert( ( kind >=0 ) && ( kind < 6 ) );
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		values_[ kind ][ i_point ] = values[ i_point ];
	}
}

/**
 *Setting value of stress on a gauss point
 *@param[in] kind Is the kind of stress to be set
 *@param[in] gauss_index Is the gaus point index where value wil  
 *@param[in] value Is the value that will be set in the class
 */
void Stress::SetValue( int kind , int gauss_index , double value ){
	assert( ( kind >= 0 ) && ( kind < 6 ) );
	assert( ( gauss_index >= 0 ) && ( gauss_index < n_gauss_points_ ) );
	values_[ kind ][ gauss_index ] = value;
}

/**
 *Setting matrix of stresses in the class
 *@param[in] values Is the matrix of values to be set in the class
 *@param[in] n_gauss_points Is the number of gauss points in the class
 */
void Stress::SetValue( double** values , int n_gauss_points ){
	n_gauss_points_ = n_gauss_points;
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
			values_[ i_kind ][ i_point ] = values[ i_kind ][ i_point ];
		}		
	}
}

/**
 *Setting number of gauss points in the class
 *@param[in] n_gauss_points Is the number of gauss points to be set in the class
 */
void Stress::SetNGaussPoints( int n_gauss_points ){
	n_gauss_points_ = n_gauss_points;
}

//UTILITIES
/**
 *Calculating averages of stress in each kind
 */
void Stress::CalculateAverages( ){
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		averages_[ i_kind ] = 0.0;
		for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
			averages_[ i_kind ] += values_[ i_kind ][ i_point ];
		}
		averages_[ i_kind ] /= (double)n_gauss_points_;
		average_flag_[ i_kind ] =  true;
	}
}

/**
 *Actualizing average of a kind of stress
 *@param[in] kind Is the kind of average stress to be actualized
 */
void Stress::ActualizeAverage( int kind ){
	assert( ( kind >= 0 ) && ( kind < 6 ) );
	averages_[ kind ] = 0.0;
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		averages_[ kind ] += values_[ kind ][ i_point ];
	}
	averages_[ kind ] /= (double)n_gauss_points_;
}

/**
 *Actualizing the average from all kind of strains in the class
 */
void Stress::ActualizeAverage(  ){
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		averages_[ i_kind ] = 0.0;
		for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
			averages_[ i_kind ] += values_[ i_kind ][ i_point ];
		}
		averages_[ i_kind ] /= (double)n_gauss_points_;
	}
}


//MODIFIERS

//CLEANERS
/**
 *Cleaning the requested memory
 */
void Stress::ClearMemory(  ){
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		if(  values_[ i_kind ]  ){
			delete[] values_[ i_kind ];
		}
	} 
}

//DEBUG
/**
 *Printing information contained in class stress
 */
void Stress::PrintInfo(  ){
	std::cout << "Number of gauss points: " << n_gauss_points_ << std::endl;
	std::cout << "Flags:      set      average: " << std::endl;
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		std::cout << set_flag_[ i_kind ] << "        " << average_flag_[ i_kind ] << std::endl;	
	}
	std::cout << "Values:    Sxx    Syy    Szz    Sxy    Syz    Sxz  :  Average" << std::endl;
	for(  int i_kind = 0  ;  i_kind < 6  ;  i_kind++  ){
		for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
			std::cout << values_[ i_kind ][ i_point ] << "    ";
		}
		std::cout << averages_[ i_kind ] << std::endl;
	}
}

/**
 *Testing if values on class were set or not
 *@param[in] kind is the kind of stress to be tested 
 *@return A bool value indicating if the values were set or not
 */
bool Stress::IsSet( int kind ){
	assert( ( kind >= 0 ) && ( kind < 6 ) );
	return set_flag_[ kind ];
}

/**
 *Testing if average was calculated or not
 *@param[in] kind is the kind of stress to be tested 
 *@return A bool value indicating if the average stress was calculated or not
 */
bool Stress::AverageIsCalculated( int kind ){
	assert( ( kind >= 0 ) && ( kind < 6 ) );
	return average_flag_[ kind ];
}


//////////////////////////////////////////////////////////////////////////////////////////
/////																VONMISES METHODS																		//
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Vonmises::Vonmises(  ){
	n_gauss_points_ = 0;
	set_flag_ = false;
	average_flag_ = false;		
	average_ = 0.0;
	values_ = NULL;
}

/**
 *Constructor receiving the vector of vonmises and the number of gauss points
 *@param[in] values Is the vector containning the vonmises to be set
 *@param[in] n_gauss_points Is the number of gauss points in the element
 */
Vonmises::Vonmises( double* values , int n_gauss_points ){
	n_gauss_points_ = n_gauss_points;
	set_flag_ = true;
	average_flag_ = true;		
	values_ = new double[ n_gauss_points_ ];
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		values_[ i_point ] = values[ i_point ];
		average_ += values_[ i_point ];
	}
	average_ /= (double)n_gauss_points_;
}

/**
 *Default destructor
 */
Vonmises::~Vonmises(  ){
	n_gauss_points_ = 0;
	set_flag_ = false;
	average_flag_ = false;		
	average_ = 0.0;
	delete[] values_;
	values_ = NULL;
}

//GETS
/**
 *Getting value of vonmisses
 *@param[in] gauss_index Is the index of the gauss point index requested
 *@return A double value with the vonmises value
 */
double Vonmises::GetValue( int gauss_index ){
	assert( ( gauss_index >= 0 ) && ( gauss_index < n_gauss_points_ ) );
	return values_[ gauss_index ];
}

/**
 *Getting array of vonmises
 *@param[out] values Is the array were values will be saved
 */
void Vonmises::GetValue( double* values ){
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		values[ i_point ] = values_[ i_point ];
	}
}

/**
 *Getting average of vonmises
 *@return A double value with the average of the vonmises for the element
 */
double Vonmises::GetAverage( ){
	return average_;
}

/**
 *Getting number of gauss points
 *@return A int value with the number of gauss points
 */
int Vonmises::GetNGaussPoints(  ){
	return n_gauss_points_;
}

//SETS
/**
 *Setting vonmises on a gaus point
 *@param[in] gauss_index Is the index where the value will be set
 *@param[in] values Is the value to be set
 */
void Vonmises::SetValue( int gauss_index , double values ){
	assert( ( gauss_index >= 0 ) && ( gauss_index < n_gauss_points_ ) );
	values_[ gauss_index ] = values;
}

/**
 *Setting values of vonmises on the class
 *@param[in] values Is the array containing the values to be set
 *@param[in] n_gauss_points Is the number of gauss points in the element
 */
void Vonmises::SetValue( double* values , int n_gauss_points ){
	if(  !values_  ){
		values_ = new double[ n_gauss_points ];
	}
	n_gauss_points_ =  n_gauss_points;
	set_flag_ = true;
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		values_[ i_point ] =  values[ i_point ];
	}
}

/**
 *Setting number of gauss points
 *@param[in] n_gauss_points Is the number of gauss points to be set in the class
 */
void Vonmises::SetNGaussPoints( int n_gauss_points ){
	n_gauss_points_ = n_gauss_points;
}

//UTILITIES
/**
 *Calculating average vonmises on the class
 */
void Vonmises::CalculateAverage( ){
	average_ = 0.0;
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		average_ += values_[ i_point ]; 	
	}
	average_ /= (double)n_gauss_points_;
	average_flag_ = true;
}

/**
 *Actualize average of vonmises
 */
void Vonmises::ActualizeAverage(  ){
	average_ = 0.0;
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		average_ += values_[ i_point ]; 	
	}
	average_ /= (double)n_gauss_points_;
	average_flag_ = true;
}

//MODIFIERS

//CLEANERS
/**
 *Cleaning the memory requested
 */
void Vonmises::ClearMemory(  ){
	if(  values_  ){
		delete[] values_;
	}
}

//DEBUG
/**
 *Printing information contained in the class
 */
void Vonmises::PrintInfo(  ){
	std::cout << "Number of gauss points: " << n_gauss_points_ << std::endl;
	std::cout << "Flags:      set      average: " << std::endl;
	std::cout << set_flag_ << "        " << average_flag_ << std::endl;	
	std::cout << "Values: " ;
	for(  int i_point = 0  ;  i_point < n_gauss_points_  ;  i_point++  ){
		std::cout << values_[ i_point ] << "  ";
	}
	std::cout << average_ << std::endl;
}

/**
 *Testing if values are set or not
 *@return A bool value indicating if vonmises were set or not
 */
bool Vonmises::IsSet(  ){
	return set_flag_;
}

/**
 *Testing if average vonmises was calculated
 *@return A bool value testing if the average was calculated or not
 */
bool Vonmises::AverageIsCalculated(  ){
	return average_flag_;
}


//////////////////////////////////////////////////////////////////////////////////////////
/////																	RESULTS METHODS																		//
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Results::Results(  ){
	strcpy( name_ , "\0" );
	ptr_ = NULL;
	n_nodes_ = 0.0;
	n_elements_ = 0.0;
	displacements_.clear();
	strains_.clear();
	stress_.clear();
	vonmises_.clear();
}

/**
 *Constructor using the name of the filecontaining the results
 *@param[in] name Is the name of the file containing the results
 *@param[in] n_nodes Is the number of nodes in the mesh readed
 *@param[in] n_elements Is the number of elements in the mesh
 */
Results::Results( char* name , size_t n_nodes , size_t n_elements ){
	n_nodes_ = n_nodes;
	n_elements_ = n_elements;

	name_ = name;
	ptr_ = fopen( name_ , "r" );
	if( !ptr_ ){
		std::cout << "Results file not opened " << std::endl;
		assert( 0 );
	}
	while( !feof( ptr_ ) ){
		char aux[ 100 ];
		fscanf( ptr_ , "%s" , aux );
		if(  !strcmp( aux , "Values" )  ){
			//Reading Displacements
			for(  size_t i_node = 0  ;  i_node < this->GetNNodes()  ;  i_node++  ){
				double values[ 3 ];
				fscanf( ptr_ , "%s" , aux );
				fscanf( ptr_ , "%s" , aux );
				values[ 0 ] = atof( aux );
				fscanf( ptr_ , "%s" , aux );
				values[ 1 ] = atof( aux );
				fscanf( ptr_ , "%s" , aux );
				values[ 2 ] = atof( aux );
				Displacement* disp = new Displacement( values );
				this->AddDisplacement( disp );
			}
			for(  int i_trash = 0  ;  i_trash < 10  ;  i_trash++  ){
				fscanf( ptr_ , "%s" , aux );
			}


			//Reading Strains
			double** values;
			values = new double*[ 6 ];
			for(  int i_row = 0  ;  i_row < 6  ;  i_row++  ){
				values[ i_row ] = new double[ 8 ];
			}
			for(  size_t i_elem = 0  ;  i_elem < this->GetNElements()  ;  i_elem++  ){
				fscanf( ptr_ , "%s" , aux );
				for(  int i_row = 0  ;  i_row < 8  ;  i_row++  ){
					for(  int i_col = 0  ;  i_col < 6  ;  i_col++  ){
						fscanf( ptr_ , "%s" , aux );
						values[ i_col ][ i_row ] = atof( aux );
					}
				}
				Strain* strain =  new Strain( values , 8 );
				this->AddStrain( strain );
			}
			for(  int i_trash = 0  ;  i_trash < 10  ;  i_trash++  ){
				fscanf( ptr_ , "%s" , aux );
			}

			//Reading Stress
			for(  size_t i_elem = 0  ;  i_elem < this->GetNElements()  ;  i_elem++  ){
				fscanf( ptr_ , "%s" , aux );
				for(  int i_row = 0  ;  i_row < 8  ;  i_row++  ){
					for(  int i_col = 0  ;  i_col < 6  ;  i_col++  ){
						fscanf( ptr_ , "%s" , aux );
						values[ i_col ][ i_row ] = atof( aux );
					}
				}
				Stress* stress =  new Stress( values , 8 );
				this->AddStress( stress );
			}
			for(  int i_trash = 0  ;  i_trash < 11  ;  i_trash++  ){
				fscanf( ptr_ , "%s" , aux );
			}
			for(  int i_row = 0  ;  i_row < 6  ;  i_row++  ){
				delete[] values[ i_row ];
			}

			//Reading Vonmisses
			values[ 0 ] = new double[ 8 ];
			for(  size_t i_elem = 0  ;  i_elem < this->GetNElements()  ;  i_elem++  ){
				fscanf( ptr_ , "%s" , aux );
				for(  int i_row = 0  ;  i_row < 8  ;  i_row++  ){
					fscanf( ptr_ , "%s" , aux );
					values[ 0 ][ i_row ] = atof( aux );
				}
				Vonmises* von =  new Vonmises( values[ 0 ] , 8 );
				this->AddVonmises( von );
			}
			delete[] values[ 0 ];
			delete[] values;
			break;
		}
	}
	fclose( ptr_ );
}

/**
 *Default destructor
 */
Results::~Results(  ){
	strcpy( name_ , "\0" );
	fclose( ptr_ );
	ptr_ = NULL;
	for(  size_t i_elem = 0  ;  i_elem < n_elements_  ;  i_elem++  ){
		strains_[ i_elem ]->ClearMemory();	
		stress_[ i_elem ]->ClearMemory();	
		vonmises_[ i_elem ]->ClearMemory();	
	}
	n_nodes_ = 0;
	n_elements_ = 0;
}

//GETS
/**
 *Getting number of nodes in the mesh
 *@return A size_t value with the number of nodes contained in the mesh
 */
size_t Results::GetNNodes(  ){
	return n_nodes_;
}

/**
 *Getting number of elements in the mesh
 *@return A size_t value with the amount of elements in the mesh
 */
size_t Results::GetNElements(  ){
	return n_elements_;
}

/**
 *Getting the magnitude of the displacement on a node in the mesh
 *@param[in] index Is the index of the node in the mesh
 *@return A double value with the magnitude of the displacement
 */
double Results::GetDisplacementMagnitude( size_t index ){
	assert( ( index >= 0 ) && ( index < n_nodes_ ) );
	return displacements_[ index ]->GetDisplacementMagnitude();
}

/**
 *Getting diaplcement from a node in one direction
 *@param[in] index Is the index in which the displacement is requested
 *@param[in] position Is direction where the displacement is requested
 *@return A double value with the displacement requested
 */
double Results::GetDisplacement( size_t index , int position ){
	assert( ( index >= 0 ) && ( index < n_nodes_ ) );
	assert( ( position >= 0 ) && ( position < 3 ) );
	return displacements_[ index ]->GetDisplacement( position );
}


void Results::GetDisplacement( size_t index , double* values ){
	assert( 0 );//NOT IMPLEMENTED
}
Displacement* Results::GetDisplacement( size_t index ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::GetStrain( size_t index , int kind_requested ,double* res_vector ){
	assert( 0 );//NOT IMPLEMENTED
}

/**
 *Getting strain from an element
 *@param[in] index Is the index of the element that contains the strain value
 *@param[in] kind_requested Is the kind of strain requested.
 *@param[in] gauss_index Is the index of the gauss point where the strain is requested
 *@return A double value witht he strain requested
 */
double Results::GetStrain( size_t index , int kind_requested ,int gauss_index ){
	return strains_[ index ]->GetValue( kind_requested , gauss_index );
}
double Results::GetAverageStrain( size_t index , int kind_requested ){
	assert( 0 );//NOT IMPLEMENTED
}
Strain* Results::GetStrain( size_t index ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::GetStress( size_t index , int kind_requested ,double* res_vector ){
	assert( 0 );//NOT IMPLEMENTED
}

/**
 *Getting stress from an element
 *@param[in] index Is the index of the element that contains the stress value
 *@param[in] kind_requested Is the kind of stress requested.
 *@param[in] gauss_index Is the index of the gauss point where the stress is requested
 *@return A double value witht he stress requested
 */
double Results::GetStress( size_t index , int kind_requested ,int gauss_index ){
	return stress_[ index ]->GetValue( kind_requested , gauss_index );
}
double Results::GetAverageStress( size_t index , int kind_requested ){
	assert( 0 );//NOT IMPLEMENTED
}
Stress* Results::GetStress( size_t index ){
	assert( 0 );//NOT IMPLEMENTED
}

/**
 *Getting vonmises from a gauss point on the element
 *@param[in] index Is the index of the element that contains the conmises value
 *@param[in] gauss_index Is the index of the gauss point that contains the vonmises
 *@return A double value with the conmises requested
 */
double Results::GetVonmises( size_t index , int gauss_index ){
	return vonmises_[ index ]->GetValue( gauss_index );
}

void Results::GetVonmises( size_t index , double* values ){
	assert( 0 );//NOT IMPLEMENTED
}
double Results::GetAverageVonmises( size_t index ){
	assert( 0 );//NOT IMPLEMENTED
}
Vonmises* Results::GetVonmises( size_t index ){
	assert( 0 );//NOT IMPLEMENTED
}

/**
 *Getting name of the file
 *@return A char value witht he name of the input results file
 */
char* Results::GetName(){
	return name_;
}

//SETS
void Results::SetDisplacement( size_t index , int position , double value ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::SetDisplacement( size_t index , double X, double Y, double Z ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::SetDisplacement( size_t index , double* values ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::SetStrain( size_t index , int kind , double* values ,int n_gauss_points ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::SetStrain( size_t index , int kind , int gauss_index , double value ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::SetStrain( size_t index , double** values , int n_gauss_points ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::SetStress( size_t index , int kind , double* values , int n_gauss_points){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::SetStress( size_t index , int kind , int gauss_index , double value ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::SetStress( size_t index , double** values , int n_gauss_points ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::SetVonmises( size_t index , int gauss_index , double values ){
	assert( 0 );//NOT IMPLEMENTED
}
void Results::SetVonmises( size_t index , double* values , int n_gauss_points ){
	assert( 0 );//NOT IMPLEMENTED
}

//UTILITIES
/**
 *Adding a new displacement
 *@param[in] disp Is th pointer to the displacement to be added in the displacements vector
 */
void Results::AddDisplacement( Displacement* disp ){
	displacements_.push_back( disp );
}

/**
 *Adding a new strain to the list of straoins
 *@param[in] strain Is the pointer to the strain added
 */
void Results::AddStrain( Strain* strain ){
	strains_.push_back( strain );
}

/**
 *Adding a new stress to the list of stress
 *@param[in] stress Is the pointer to the stress added
 */
void Results::AddStress( Stress* stress ){
	stress_.push_back( stress );
}

/**
 *Adding a new vonmises in the list of vonmises
 *@param[in] von Is the pointer to the vonmises to be added
 */
void Results::AddVonmises( Vonmises* von ){
	vonmises_.push_back( von );
}
		
//MODIFIERS


//DEBUG
/**
 *Printing information contained in class results
 */
void Results::PrintInfo(){
	std::cout << "Name: " << name_ << std::endl;
	std::cout << "Nodos: " << n_nodes_ << std::endl;
	std::cout << "Elementos: " << n_elements_ << std::endl;
	std::cout << "Displacements: " << displacements_.size() << std::endl; 
	std::cout << "Strains: " << strains_.size() << std::endl;
	std::cout << "Stress: " << stress_.size() << std::endl;
	std::cout << "Vonmises: " << vonmises_.size() << std::endl;
}






























