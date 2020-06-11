#pragma once

//C system files
#include <math.h>
#include <assert.h>

//C++ system files
#include <vector>
#include <iostream>

//Other libraries
using namespace std;

//AUXILIAR
void CalcVector( double* dest , double* P0 , double* P1 );
void SubstractFromVector( double* dest , double* source , double* subtra );
double Dot( double* vector1 , double* vector2 );
void Cross( double* dest , double* vector1 , double* vector2 );
void FindMinMax( double& x0 , double& x1 , double& x2 , double& min , double& max );
void Normalize( double *vec );
double Magnitude( double *vec );
void GenerateOuterOctreePoint( double* P1 );
double DistancePointPoint( double *P , double *Q );

//TESTS
bool TriangleIntersectsBox( double* box_center , double* box_half_size , double** triverts );
bool PlaneIntersectsBox( double* normal , double* vertex , double* box_half_size );
bool LineIntersectsTriangle( double* P0 , double* P1 , double* T0 , double* T1 , double* T2 );
bool AABBvsAABB( double* A_center , double A_radius , double* B_center , double B_radius );
bool SphereIntersectsCell(  double *min_coord , double *max_coord , double *center , double radius  );
bool PointIsContainedInAABB( double* P , double* min_coord , double* max_coord );
bool AABBvsAABBDifferentSize( double* A_center , double* A_width , double* B_center , double* B_width );

//X-TESTS
bool AxisTest_X01( double& a, double& b, double& fa, double& fb, double& p0, double& p2,
                   double& min, double& max, double& rad, double* v0, double* v2, 
									 double* box_half_size );
bool AxisTest_X2( double& a, double& b, double& fa, double& fb, double& p0, double& p1,
									double& min, double& max, double& rad, double* v0, double* v1, 
									double* box_half_size );
//Y-TESTS
bool AxisTest_Y02( double& a, double& b, double& fa,  double& fb, double& p0,  double& p2,
									 double& min, double& max, double& rad, double* v0, double* v2, 
									 double* box_half_size );
bool AxisTest_Y1( double& a, double& b, double& fa,  double& fb, double& p0, double& p1,
									double& min, double& max, double& rad, double* v0, double* v1, 
									double* box_half_size );
//Z-TESTS
bool AxisTest_Z12( double& a, double& b, double& fa, double& fb, double& p1, double& p2,
									 double& min, double& max, double& rad, double* v1, double* v2, 
									 double* box_half_size );
bool AxisTest_Z0( double& a, double& b, double& fa, double& fb, double& p0,  double& p1,
									double& min, double& max, double& rad, double* v0, double* v1, 
									double* box_half_size );


//UTILITIES
double SignedDistanceFromFaceToPoint(double *P,double *A,double *B,double *C);
void CalculateBaricentricCoordinates(double *Q,double *A,double *B,double *C,double *coordinates);
int FindStartAndEndOfIntersection( double* P0 , double* P1 , double* T0 , double* T1 , 
																	 double* T2 , double* X , int dim );
bool FindIntersectionLineTriangle( double* P0 , double* P1 , double* T0 , double* T1 , 
																	double* T2 , double* X , double eps , int dim );
bool LineIntersectLine( double* P0 , double* P1 , double* Q0 , double* Q1 , 
                       double* p , int dim);

