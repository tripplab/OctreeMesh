#include "geometry.h"

double min_rand[ 26 ][ 3 ] = { { -1.0 ,  0.0 ,  0.0 } ,//0
															 {  1.0 ,  0.0 ,  0.0 } ,
															 {  0.0 , -1.0 ,  0.0 } ,
															 {  0.0 ,  1.0 ,  0.0 } ,
															 {  0.0 ,  0.0 , -1.0 } ,
															 {  0.0 ,  0.0 ,  1.0 } ,//5
															 { -1.0 , -1.0 ,  0.0 } ,
															 {  1.0 , -1.0 ,  0.0 } ,
															 { -1.0 ,  1.0 ,  0.0 } ,
															 {  1.0 ,  1.0 ,  0.0 } ,
															 { -1.0 ,  0.0 , -1.0 } ,//10
															 {  1.0 ,  0.0 , -1.0 } ,
															 { -1.0 ,  0.0 ,  1.0 } ,
															 {  1.0 ,  0.0 ,  1.0 } ,
															 {  0.0 , -1.0 , -1.0 } ,
															 {  0.0 ,  1.0 , -1.0 } ,//15
															 {  0.0 , -1.0 ,  1.0 } ,
															 {  0.0 ,  1.0 ,  1.0 } ,
															 { -1.0 , -1.0 , -1.0 } ,
															 {  1.0 , -1.0 , -1.0 } ,
															 {  1.0 ,  1.0 , -1.0 } ,//20
															 { -1.0 ,  1.0 , -1.0 } ,
															 { -1.0 , -1.0 ,  1.0 } ,
															 {  1.0 , -1.0 ,  1.0 } ,
															 {  1.0 ,  1.0 ,  1.0 } ,
															 { -1.0 ,  1.0 ,  1.0 } };//25

//AUXILIAR
/**
 *Calculating vector P0-P1
 *@param[out] dest Is the result vector from P0-P1
 *@param[in] P0 Is the begining point of vector
 *@param[in] P1 Is the end point of vector 
 */
void CalcVector( double* dest , double* P0 , double* P1 ){
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		dest[ i_dim ] = P1[ i_dim ] - P0[ i_dim ];
	}
}

/**
 *Substracting an amount from a vector
 *@param[out] dest Is the result of the operation
 *@param[in] source Is the initial vector value
 *@param[in] subtra Is the amopunt to substract from vector
 */
void SubstractFromVector( double* dest , double* source , double* subtra ){
	for(  int i_dim = 0  ; i_dim < 3  ;  i_dim++  ){
		dest[ i_dim ] = source[ i_dim ] - subtra[ i_dim ];
	}
}

/**
 *Dot product of two vectors
 *@param[in] vector1 Is the first vector
 *@param[in] vector2 Is the second vector
 *@return A double value with the dot product
 */
double Dot( double* vector1 , double* vector2 ){
	double sum = 0.0;
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		sum += ( vector1[ i_dim ] * vector2[ i_dim ] );
	}
	return sum;
}

/**
 *Calculates cross product between vectors
 *@param[out] dest Resultant vector from operation
 *@param[in] vector1 Is the first vector 
 *@param[in] vector2 Is the second vector
 */
void Cross( double* dest , double* vector1 , double* vector2 ){
	dest[ 0 ] = ( vector1[ 1 ] * vector2[ 2 ] ) - ( vector1[ 2 ] * vector2[ 1 ] );
	dest[ 1 ] = ( vector1[ 2 ] * vector2[ 0 ] ) - ( vector1[ 0 ] * vector2[ 2 ] );
	dest[ 2 ] = ( vector1[ 0 ] * vector2[ 1 ] ) - ( vector1[ 1 ] * vector2[ 0 ] ); 
}

/**
 *Finds min and max from 3 values
 *@param[in] x0 Is the first value
 *@param[in] x1 Is the second value
 *@param[in] x2 Is the third value
 *@param[out] min Is the minimum between x0, x1 and x2
 *@param[out] max Is the maximum between x0, x1 and x2
 */
void FindMinMax( double& x0 , double& x1 , double& x2 , double& min , double& max ){
	min = max = x0;
	if( x1 < min ){ 
		min = x1;
	}
	if( x1 > max ){
		max = x1;
	}
	if( x2 < min ){
		min = x2;
	}
	if( x2 > max ){
		max = x2;
	}
}

/**
 *Normalizing a vector
 *@param[out] vec Is the input vector and in it is stored the normalized vector
 */
void Normalize( double* vec ){
    double val;
    val = Magnitude( vec );
		assert( !( val == 0.0 ) );
    for(  int i_element = 0  ;  i_element < 3  ;  i_element++  ){
        vec[ i_element ] = vec[ i_element ] / val;
    }
}

/**
 *Calculating norm of a vector
 *@param[in] vec Is the vector pointer to calculate its magnitude
 *@return A double value with the vector magnitude 
 */
double Magnitude( double* vec ){
    //||Vec||=sqrt(sum(pow(vec[1],2)))
    double val = 0.0;
    for(  int i_element = 0  ;  i_element < 3  ;  i_element++  ){
        val += ( vec[ i_element ] * vec[ i_element ] );
    }
		assert( val >= 0.0 );
    val = pow( val , 0.5 );
    return val;
}

/**
 *Generating a random point outside the octree
 *@param[out] P1 is the point generated
 */
void GenerateOuterOctreePoint( double* P ){
	int rnd = rand() % 26;
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		P[ i_dim ] = ( (double)rand()/(double)RAND_MAX ) + min_rand[ rnd ][ i_dim ]; 
	}
}

/**
 *Calculating distance between two points
 *@param[in] P is the first point
 *@param[in] Q Is the second point
 *@return A double value indicating the distance between two points
 */
double DistancePointPoint( double *P , double *Q ){
    double dist = 0.0;
    for(  int i_element = 0  ;  i_element < 3  ;  i_element++  ){
        dist += pow( P[ i_element ] - Q[ i_element ] , 2 );
    }
    dist = pow( dist , 0.5 );
    return dist;
}

//TESTS
/**
 *Testing intersection triangle-AABB(Axis aligned bounding box)
 * AABB-triangle overlap test code                      
 * by Tomas Akenine-Möller                              
 * Function: int triBoxOverlap(float boxcenter[3],      
 *          float boxhalfsize[3],float triverts[3][3]); 
 * History:                                             
 *   2001-03-05: released the code in its first version 
 *   2001-06-18: changed the order of the tests, faster 
 *                                                      
 * Acknowledgement: Many thanks to Pierre Terdiman for  
 * suggestions and discussions on how to optimize code. 
 * Thanks to David Hunt for finding a ">="-bug!         
 *@param[in] box_center Is the center of the AABB
 *@param[in] box_half_size Is the length from the AABB center to the sides of AABB
 *@param[in] triverts Are the triangle vertex by rows
 */
bool TriangleIntersectsBox( double* box_center , double* box_half_size , double** triverts ){
	double v0[ 3 ],v1[ 3 ],v2[ 3 ];
	double min,max,p0,p1,p2,rad,fex,fey,fez;
	double normal[ 3 ],e0[ 3 ],e1[ 3 ],e2[ 3 ];

	

	//Moving everything so the box center is on {0,0,0}
	SubstractFromVector( v0 , triverts[ 0 ] , box_center );
	SubstractFromVector( v1 , triverts[ 1 ] , box_center );
	SubstractFromVector( v2 , triverts[ 2 ] , box_center );
	
	//Getting triangle edges
	CalcVector( e0 , v0 , v1 );
	CalcVector( e1 , v1 , v2 );
	CalcVector( e2 , v2 , v0 );

	//Testing 9 test first
  fex = fabs( e0[ 0 ] );
  fey = fabs( e0[ 1 ] );
  fez = fabs( e0[ 2 ] );
  if(  !AxisTest_X01( e0[ 2 ] , e0[ 1 ] , fez , fey , p0 , p2 , min , max , rad , v0 , v2 , box_half_size )  ){
		return false;
	}
  if(  !AxisTest_Y02( e0[ 2 ] , e0[ 0 ] , fez , fex , p0 , p2 , min , max , rad , v0 , v2 , box_half_size )  ){
		return false;
	}
  if(  !AxisTest_Z12( e0[ 1 ] , e0[ 0 ] , fey , fex , p1 , p2 , min , max , rad , v1 , v2 , box_half_size ) ){
		return false;
	}

	fex = fabs( e1[ 0 ] );
	fey = fabs( e1[ 1 ] );
	fez = fabs( e1[ 2 ] );
  if(  !AxisTest_X01( e1[ 2 ] , e1[ 1 ] , fez , fey , p0 , p2 , min , max , rad , v0 , v2 , box_half_size )  ){
		return false;
	}
  if(  !AxisTest_Y02( e1[ 2 ] , e1[ 0 ] , fez , fex , p0 , p2 , min , max , rad , v0 , v2 , box_half_size )  ){
		return false;
	}
  if(  !AxisTest_Z0( e1[ 1 ] , e1[ 0 ] , fey , fex , p0 , p1 , min , max , rad , v0 , v1 , box_half_size )  ){
		return false;
	}

	fex = fabs( e2[ 0 ] );
	fey = fabs( e2[ 1 ] );
	fez = fabs( e2[ 2 ] );
  if(  !AxisTest_X2( e2[ 2 ] , e2[ 1 ] , fez , fey , p0 , p1 , min , max , rad , v0 , v1 , box_half_size )  ){
		return false;
	}
  if(  !AxisTest_Y1( e2[ 2 ] , e2[ 0 ] , fez , fex , p0 , p1 , min , max , rad , v0 , v1 , box_half_size )  ){
		return false;
	}
  if(  !AxisTest_Z12( e2[ 1 ] , e2[ 0 ] , fey , fex , p1 , p2 , min , max , rad , v1 , v2 , box_half_size )  ){
		return false;
	}

	// test in X-direction
	FindMinMax( v0[ 0 ] , v1[ 0 ] , v2[ 0 ] , min , max );
	if(  ( min > box_half_size[ 0 ] ) || ( max < -box_half_size[ 0 ] )  ){ 
		return false;
	}

	//test in Y-direction
	FindMinMax( v0[ 1 ] , v1[ 1 ] , v2[ 1 ] , min , max );
	if( ( min > box_half_size[ 1 ] ) || ( max < -box_half_size[ 1 ] )  ){
		return false;
	}

	//test in Z-direction
	FindMinMax( v0[ 2 ] , v1[ 2 ] , v2[ 2 ] , min , max );
	if(  ( min > box_half_size[ 2 ] ) || ( max < -box_half_size[ 2 ] )  ){
		return false;
	}


	Cross( normal , e0 , e1 );
	if(  !PlaneIntersectsBox( normal , v0 , box_half_size )  ){
		return false;
	}
	return true;
}

/**
 *Testing if plane intersects AABB
 *@param[in] normal Is the normal of the plane
 *@param[in] vertex Is a point of the plane
 *@param[in] box_half_size Is the size og the AABB in each dimension
 *@return a bool value indicating if intersects or not
 */
bool PlaneIntersectsBox( double* normal , double* vertex , double* box_half_size ){
  double vmin[ 3 ],vmax[ 3 ],v;
  for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
    v = vertex[ i_dim ];
    if(  normal[ i_dim ] > 0.0  ){
      vmin[ i_dim ] = -box_half_size[ i_dim ] - v;
      vmax[ i_dim ] =  box_half_size[ i_dim ] - v;
    }else{
      vmin[ i_dim ] =  box_half_size[ i_dim ] - v;
      vmax[ i_dim ] = -box_half_size[ i_dim ] - v;
    }
  }
  if(  Dot( normal , vmin ) > 0.0  ){
 		return false;
	}
  if(  Dot( normal , vmax ) >= 0.0  ){
		return true;
	}
  return false;
}

/**
 *Testing for line-triangle intersection
 *This method test if a line intersects on a triangle
 *@param[in] P0 Is the beginning point of the line
 *@param[in] P1 Is the end point of the line
 *@param[in] T0 Is the first node of the triangle
 *@param[in] T1 Is the second node of the triangle
 *@param[in] T2 Is the third node of the triangle
 *@return A bool value indicating if the triangle intersects or not
 */
bool LineIntersectsTriangle( double* P0 , double* P1 , double* T0 , double* T1 , double* T2 ){

	double dist1 , dist2 , val , eps , Q[ 3 ] , BAR[ 3 ];
	//calculating signed distance from points Q1 Y Q2 to the plane defined by triangle P1-P2-P3
	dist1 = SignedDistanceFromFaceToPoint( P0 , T0 , T1 , T2 );
	dist2 = SignedDistanceFromFaceToPoint( P1 , T0 , T1 , T2 );
	val = dist1 * dist2;
	//if the signed distance product is negative then, line intersects the plane of triangular face
	if( val < 0.0 ){
		eps = ( dist1 ) / ( dist1 - dist2 );
	}else{
		return false;
	}
	//Calculating intersection point from line with plane
	for(  int i_element = 0  ;  i_element < 3  ;  i_element++  ){
		Q[ i_element ] = P0[ i_element ] + ( eps * ( P1[ i_element ] - P0[ i_element ] ) );
	}
	//calculating area coordinates from intersection point to verify if it is inside triangle P1-P2-P3
	CalculateBaricentricCoordinates( Q , T0 , T1 , T2 , BAR );
	if( ( BAR[ 0 ] >= -1e-10 ) && ( BAR[ 1 ] >= -1e-10 ) && ( BAR[ 2 ] >= -1e-10 ) ){
		return true;
	}
	return false;
}	

/**
 *Performing AABB vs AABB intersection test
 *This methos performs AABB vs AABB test by using the center and halfsize of each AABB
 *@param[in] A_center Is the center of the box A
 *@param[in] A_radius is the radius for cell A
 *@param[in] B_center Is the center of the box B
 *@param[in] B_radius is the radius for cell B
 *@return A bool value indicating if exist or not the collision
 */
bool AABBvsAABB( double* A_center , double A_radius , double* B_center , double B_radius ){
	double value = A_radius + B_radius;

	if(  fabs( A_center[ 0 ] - B_center[ 0 ] ) > value ) return false;
	if(  fabs( A_center[ 1 ] - B_center[ 1 ] ) > value ) return false;
	if(  fabs( A_center[ 2 ] - B_center[ 2 ] ) > value ) return false;

	return true;
}

/**
 *Sphere intersects cell
 *
 *This functions calculates if a sphere intersects a cell
 *
 *@param[in] min_coord Is the lower frontal inferior point of a AABB
 *@param[in] max_coord Is the upper back superior point of a AABB
 *@param[in] center Is the center coordinate from the sphere
 *@param[in] radius Is the radius of the sphere.
 *@return A bool value indicating if intersection exists or not
 */
bool SphereIntersectsCell(  double *min_coord , double *max_coord , double *center , double radius  ){
  double dist = 0.0;
  for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
    double edge = center[ i_dim ] - min_coord[ i_dim ];
    if(  edge < 0.0  ){
      if(  edge < -radius  ){
        return false;
      }else{
        dist += ( edge * edge );
      }
    }else{
      edge = center[ i_dim ] - max_coord[ i_dim ];;
      if(  edge > 0.0  ){
        if(  edge > radius  ){
          return false; 
        }else{
          dist += ( edge * edge );
        }
      }
    }
  }
  if(  dist <= ( radius * radius )  ){
    return true;
  }
  return false;
}

/**
 *Testing if a point is contained in an AABB
 *@param[in] P Is the point to be tested
 *@param[in] min_coord Is the AABB minimum coordinate (inferior,left,frontal)
 *@param[in] max_coord Is the AABB maximum coordinate (superior,right,back)
 *@return A bool value indicating if point is contained on cell or not
 */
bool PointIsContainedInAABB( double* P , double* min_coord , double* max_coord ){

	for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
		if(  ( P[ i_pos ] < min_coord[ i_pos ] ) || ( P[ i_pos ] > max_coord[ i_pos ] )  ) return false;
	}
	return true;
}

/**
 *Function to testt the intersection between 2 boxes with different length from each dimention
 *@param[in] A_center Is the center of box A.
 *@param[in] A_width Is half of the length of the box dimention from box A.
 *@param[in] B_center Is the center of box B.
 *@param[in] B_width Is half of the length of the box dimention from box B.
 *@return A bool value indicating if intersection exist or not
 */
bool AABBvsAABBDifferentSize( double* A_center , double* A_width , double* B_center , double* B_width ){
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		if( fabs( A_center[ i_dim ] - B_center[ i_dim ] ) > ( A_width[ i_dim ] + B_width[ i_dim ] ) ){
			return false;	
		}	
	}
	return true;
}

//X-TESTS
bool AxisTest_X01( double& a, double& b, double& fa, double& fb, double& p0, double& p2,
                   double& min, double& max, double& rad, double* v0, double* v2, 
									 double* box_half_size ){
	p0 = ( a * v0[ 1 ] ) - ( b * v0[ 2 ] );
	p2 = ( a * v2[ 1 ] ) - ( b * v2[ 2 ] );
	if(  p0 < p2  ){
		min = p0;
		max = p2;
	}else{
		min = p2;
		max = p0;
	}
  rad = ( fa * box_half_size[ 1 ] ) + ( fb * box_half_size[ 2 ] );
  if(  ( min > rad ) || ( max < -rad )  ){
		return false;
	}else{
		return true;
	}
}

bool AxisTest_X2( double& a, double& b, double& fa, double& fb, double& p0, double& p1,
									double& min, double& max, double& rad, double* v0, double* v1, 
									double* box_half_size ){
	p0 = ( a * v0[ 1 ] ) - ( b * v0[ 2 ] );
	p1 = ( a * v1[ 1 ] ) - ( b * v1[ 2 ] );
	if(  p0 < p1  ){
		min = p0;
		max = p1;
	}else{
		min = p1;
		max = p0;
	}
	rad = ( fa * box_half_size[ 1 ] ) + ( fb * box_half_size[ 2 ] );
	if(  ( min > rad ) || ( max < -rad )  ){ 
		return false;
	}else{
		return true;
	}
}

//Y-TESTS
bool AxisTest_Y02( double& a, double& b, double& fa,  double& fb, double& p0,  double& p2,
									 double& min, double& max, double& rad, double* v0, double* v2, 
									 double* box_half_size ){
	p0 = ( -a * v0[ 0 ] ) + ( b * v0[ 2 ] );
	p2 = ( -a * v2[ 0 ] ) + ( b * v2[ 2 ] );
	if(  p0 < p2  ){
		min = p0;
		max = p2;
	}else{
		min = p2;
		max = p0;
	}
	rad = ( fa * box_half_size[ 0 ] ) + ( fb * box_half_size[ 2 ] );
	if(  ( min > rad ) || ( max < -rad )  ){ 
		return false;
	}else{ 
		return true;
	}
}

bool AxisTest_Y1( double& a, double& b, double& fa,  double& fb, double& p0, double& p1,
									double& min, double& max, double& rad, double* v0, double* v1, 
									double* box_half_size ){
	p0 = ( -a * v0[ 0 ] ) + ( b * v0[ 2 ] );
	p1 = ( -a * v1[ 0 ] ) + ( b * v1[ 2 ] );
	if(  p0 < p1  ){
		min = p0;
		max = p1;
	}else{
		min = p1;
		max = p0;
	}
	rad = ( fa * box_half_size[ 0 ] ) + ( fb * box_half_size[ 2 ] );
	if(  ( min > rad ) || ( max < -rad )  ){ 
		return false;
	}else{ 
		return true;
	}
}

//Z-TESTS
bool AxisTest_Z12( double& a, double& b, double& fa, double& fb, double& p1, double& p2,
									 double& min, double& max, double& rad, double* v1, double* v2, 
									 double* box_half_size ){
	p1 = ( a * v1[ 0 ] ) - ( b * v1[ 1 ] );
	p2 = ( a * v2[ 0 ] ) - ( b * v2[ 1 ] );
	if(  p2 < p1  ){
		min = p2;
		max = p1;
	}else{
		min = p1;
		max = p2;
	}
	rad = ( fa * box_half_size[ 0 ] ) + ( fb * box_half_size[ 1 ] );
	if(  ( min > rad ) || ( max < -rad )  ){ 
		return false;
	}else{ 
		return true;
	}
}

bool AxisTest_Z0( double& a, double& b, double& fa, double& fb, double& p0,  double& p1,
									double& min, double& max, double& rad, double* v0, double* v1, 
									double* box_half_size ){
	p0 = ( a * v0[ 0 ] ) - ( b * v0[ 1 ] );
	p1 = ( a * v1[ 0 ] ) - ( b * v1[ 1 ] );
	if(  p0 < p1  ){
		min = p0;
		max = p1;
	}else{
		min = p1;
		max = p0;
	}
	rad = ( fa * box_half_size[ 0 ] ) + ( fb * box_half_size[ 1 ] );
	if(  ( min > rad ) || ( max < -rad )  ){
		return false;
	}else{ 
		return true;
	}
}


//UTILITIES
/**
 *This method calculates the signed distance from a point to a triangular face
 *The signed distance means that calculates the distance from the point to the plane and 
 *see the sign of the projection to the plane containing the triangular face
 *@param[in] P Is the point to be used in the distance calculation
 *@param[in] A Is the first coordinate of the triangle 
 *@param[in] B Is the second coordinate of the triangle
 *@param[in] C Is the third coordinate of the triangle
 */
double SignedDistanceFromFaceToPoint(double *P,double *A,double *B,double *C){
    double dist;
    //Ccalculating normal vector to traingle ABC
    double AB[ 3 ] , AC[ 3 ] , N[ 3 ];
    CalcVector( AB , A , B );
    CalcVector( AC , A , C );
    Cross( N , AB , AC );
    Normalize( N );
		//calculating alfa to calculate the projection from point P to triangle ABC
    double alfa,AP[ 3 ];
    CalcVector( AP , A , P );
    alfa = Dot( AP , N );
		//calculating point p projected to plane ABC, this new point is called Q
    double AQ[ 3 ] , Q[ 3 ], QP[ 3 ];
		for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
	    AQ[ i_dim ] = AP[ i_dim ] - ( alfa * N[ i_dim ] );
	    Q[ i_dim ] = AQ[ i_dim ] + A[ i_dim ];
			QP[ i_dim ] = P[ i_dim ] - Q[ i_dim ];
		}
    dist = Dot( N , QP );
    return dist;
}

/**
 *Calculating baricentric coordinates
 *This function calculates the coordinates of area and store it on a parameter received
 *@param[in] Q Is the point located on the plate that contains the triangle
 *@param[in] A Is the first coordinate of the triangle
 *@param[in] B Is the second coordinate of the triangle
 *@param[in] C Is the third coordinate of the triangle
 *@param[out] coordinates Is the vector where the area coordinates are stored 
 */
void CalculateBaricentricCoordinates(double* Q,double* A,double* B,double* C,double* coordinates){
    //calculating vector normal to triangle ABC
    double AB[ 3 ] , AC[ 3 ] , N[ 3 ];
    CalcVector( AB , A , B );
    CalcVector( AC , A , C );
    Cross( N , AB , AC );
    //normal of QAB area coordinate a
    double QA[ 3 ] , QB[ 3 ] , NA[ 3 ];
    CalcVector( QA , Q , A );
    CalcVector( QB , Q , B );
    Cross( NA , QA , QB );
    //normal of QBC area coordinate b
    double  QC[ 3 ] , NB[ 3 ];
    CalcVector( QC , Q , C );
    Cross( NB , QB , QC );
    //normal of QCA area coordinate c
    double NC[3];
    Cross( NC , QC , QA );
		//calculating each one of the subtriangle areas    
		double areaT;
		areaT = Dot( N , N );		
    coordinates[ 0 ] = Dot( N , NA )/areaT;
    coordinates[ 1 ] = Dot( N , NB )/areaT;
    coordinates[ 2 ] = Dot( N , NC )/areaT;
}

/**
 *Looking for the intersection between a triangle and a line that are numerically parallel
 *@param[in] P0 Is the initial node of the line 
 *@param[in] P1 Is the final node of the line
 *@param[in] T0 Is the first node of the triangle
 *@param[in] T1 Is the second node of the triangle
 *@param[in] T2 Is the third node of the triangle
 *@param[out] X Is the vector where the intersection is stored
 *@param[in] dim Is the dimention where the ray is traced
 */
int FindStartAndEndOfIntersection( double* P0 , double* P1 , double* T0 , double* T1 , 
																	 double* T2 , double* X , int dim ){
/*CASES
  CASE 0: Over        CASE 1: Over 2         CASE2: Over        CASE3: No 
     vertex              lines                 axis            intersection
  ______\/______                                                 /
        /\              \   /\                  /\              /   /\
       /  \              \ /  \                /  \            /   /  \
      /    \              X    \              /    \          /   /    \
     /      \            / \    \            /      \        /   /      \
    /________\          /___X____\      ____XXXXXXXXXX____  /   /________\
                             \
                              \                                               */
	double p_inter;
	X[ 0 ] = 1e20;
	X[ 1 ] = 1e20;
	int flag[ 2 ];
	flag[ 0 ] = 0;    flag[ 1 ] = 0;
	if(  LineIntersectLine( T0 , T1 , P0 , P1 , &p_inter , dim )  ){
		flag[ 0 ] = 1;
		X[ 0 ] = p_inter;
	}
	if(  LineIntersectLine( T0 , T2 , P0 , P1 , &p_inter , dim )  ){
		if(  flag[ 0 ] == 0  ){
			flag[ 0 ] = 1;
			X[ 0 ] = p_inter;
		}else{
			flag[ 1 ] = 1;
			X[ 1 ] = p_inter;
		}
	}
	if(  LineIntersectLine( T1 , T2 , P0 , P1 , &p_inter ,dim )  ){
		if(  flag[ 0 ] == 0  ){
			flag[ 0 ] = 1;
			X[ 0 ] = p_inter;
		}else{
			if(  flag[ 1 ] == 0   ){
				flag[ 1 ] = 1;
				X[ 1 ] = p_inter;
			}
		}
	}
	return ( flag[ 0 ] && flag[ 1 ] );
}

/**
 *Finding intersection between line and triangular face
 *@param[in] P0 Is the initial node of the line
 *@param[in] P1 Is the final node of the line
 *@param[in] T0 Is the first node of the triangle
 *@param[in] T1 Is the second node of the triangle
 *@param[in] T2 Is the third node of the triangle
 *@param[out] X Is where the intersection is stored
 *@param[in] dim Is the direction where the ray is traced
 */
bool FindIntersectionLineTriangle( double* P0 , double* P1 , double* T0 , double* T1 , 
																	double* T2 , double* X , double eps , int dim ){
/*
  CASE 0: over        CASE 1: Over              CASE 2:           case 3: No 
     vertex                axis                 Interior          intersection
  ______\/______                                                 ------+------
        /\                  /\                     /\                  /\
       /  \         FRONT  /  \  BACK      FRONT  /  \ BACK    FRONT  /  \  BACK
      /    \        ------+----\----       ------/--+-\----          /    \
     /      \            /      \               /      \            /      \
    /________\          /________\             /________\          /________\
                             
                                                                             */
	double Q[ 3 ] , BAR[ 3 ];
	X[ 0 ] = 1e20;
	X[ 1 ] = 1e20;
	Q[ 0 ] = P0[ 0 ] + ( eps * ( P1[ 0 ] - P0[ 0 ] ) );
	Q[ 1 ] = P0[ 1 ] + ( eps * ( P1[ 1 ] - P0[ 1 ] ) );
	Q[ 2 ] = P0[ 2 ] + ( eps * ( P1[ 2 ] - P0[ 2 ] ) );
	CalculateBaricentricCoordinates( Q , T0 , T1 , T2 , BAR );
	if( ( BAR[ 0 ] >= -1e-10 ) && ( BAR[ 1 ] >= -1e-10 ) && ( BAR[ 2 ] >= -1e-10 ) ){
		X[ 0 ] = Q[ dim ];
		return true;
	}
	return false;
}

/**
 *Searching for intersection of line vs line
 *Based on F.S. Hill Jr., Graphics Gems IV, PP.138-148.
 *NOTE: This function works only for cartesian rays traced on dimention dim
 *@param[in] P0 Is the initial point of line one
 *@param[in] P1 Is the final point of line one
 *@param[in] Q0 Is the initial point of line two
 *@param[in] Q1 Is the final point of line two
 *@param[out] p Is the pointer to the memory where the intersection point is stored
 *@param[in] dim Is the dimention in which the intersection point is searched
 *@return A bool value indicating if a intersection has been founded or not
 */
bool LineIntersectLine( double* P0 , double* P1 , double* Q0 , double* Q1 , 
                       double* p , int dim){
	bool flag = false;
	double tol = 1e-10;
	//If lines are parallel then does not have sense to check if the intersection exist
	int indexes[ 3 ][ 2 ] = { { 1 , 2 } , { 0 , 2 } , { 0 , 1 } };
	if(  ( fabs( P0[ indexes[ dim ][ 0 ] ] - P1[ indexes[ dim ][ 0 ] ] ) < tol ) && 
			 ( fabs( P0[ indexes[ dim ][ 1 ] ] - P1[ indexes[ dim ][ 1 ] ] ) < tol )  ){
		return flag;
	}
	double A[ 3 ] , B[ 3 ] , C[ 3 ] , S , cross_ab[ 3 ] , cross_cb[ 3 ] , X[ 3 ];
	CalcVector( A , P0 , P1 ); /*P0-P1*/
	CalcVector( B , Q0 , Q1 ); /*Q0-Q1*/
	CalcVector( C , P0 , Q0 ); /*P0-Q0*/
	Cross( cross_cb , C , B );
	Cross( cross_ab , A , B );
	S = ( Dot( cross_cb , cross_ab ) / pow( Magnitude( cross_ab ) , 2 ) ); 
	X[ 0 ] = P0[ 0 ] + ( A[ 0 ] * S );
	X[ 1 ] = P0[ 1 ] + ( A[ 1 ] * S );
	X[ 2 ] = P0[ 2 ] + ( A[ 2 ] * S );
	/*Intersection between two lines has beencalculated, now is checked if intersection  
	point is between the initial and final point of line Q0-Q1, so, ar calculated the 
	distances d(X-Q0), d(X-Q1), d(Q0-Q1) and if: 
	d(Q0-Q1) = d(X-Q0) + d(X-Q1), then the intersection is between Q0-Q1*/
	double dist_1 , dist_2 , dist_T;
	dist_1 = DistancePointPoint( X , P0 );
	dist_2 = DistancePointPoint( X , P1 );
	dist_T = DistancePointPoint( P0 , P1 );
	if(  ( fabs( dist_T - ( dist_1 + dist_2 ) ) ) < tol  ){
		flag = true;
		(*p) = X[ dim ];
	}
	return flag;
}





