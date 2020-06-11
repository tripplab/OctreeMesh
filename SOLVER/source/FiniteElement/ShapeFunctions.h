// ShapeFunctions.h
// Copyright (C) 2011 Miguel Vargas (miguel.vargas@gmail.com)
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef _ShapeFunctions_h_
#define _ShapeFunctions_h_

#include <Basic/Memory.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>
#include <FiniteElement/Mesh.h>
#include <FiniteElement/Nodes.h>
#include <math.h>

#ifndef sgn
#define sgn(x) (((x) < 0) ? -1 : ((x) > 0))
#endif


template <typename T>
class ShapeFunctions
{
	public:

		class Exception : public ::Exception {};

		class ExceptionInvalidShape : public Exception {};

		static Exception exception;

		static ExceptionInvalidShape exception_invalid_shape;

		const Mesh& mesh;
		const Nodes<T>& nodes;


	private:

		void (ShapeFunctions::* element_shape_functions)(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw();
		bool (ShapeFunctions::* element_shape_functions_global)(int element_id, const Vector<T>& coordinate, Vector<T>& N) const throw();
		void (ShapeFunctions::* facet_shape_functions)(const Vector<int>& facet_nodes, const Vector<T>& point, Vector<T>& N, T& det_J) const throw();
		T (ShapeFunctions::* element_det_j)(int element_id, const Vector<T>& point) const throw();


	public:

		ShapeFunctions(const Mesh& mesh, const Nodes<T>& nodes) throw(ExceptionInvalidShape)
		:	mesh(mesh),
			nodes(nodes)
		{
			// Set shape function
			switch (mesh.element_type)
			{
				case shape_triangle:
				{
					switch (mesh.nodes_per_element)
					{
						case 3:
						{
							element_shape_functions = &ShapeFunctions::Triangle3;
							element_shape_functions_global = &ShapeFunctions::Triangle3Global;
							facet_shape_functions = &ShapeFunctions::FacetLinear2;
							element_det_j = &ShapeFunctions::Triangle3DetJ;
							break;
						}
						case 6:
						{
							element_shape_functions = &ShapeFunctions::Triangle6;
							element_shape_functions_global = &ShapeFunctions::Triangle6Global;
							facet_shape_functions = &ShapeFunctions::FacetLinear3;
							element_det_j = &ShapeFunctions::Triangle6DetJ;
							break;
						}
						default:
						{
							Throw(ShapeFunctions::exception_invalid_shape);
						}
					}
					break;
				}
				case shape_quadrilateral:
				{
					switch (mesh.nodes_per_element)
					{
						case 4:
						{
							element_shape_functions = &ShapeFunctions::Quadrilateral4;
							element_shape_functions_global = &ShapeFunctions::Quadrilateral4Global;
							facet_shape_functions = &ShapeFunctions::FacetLinear2;
							element_det_j = &ShapeFunctions::Quadrilateral4DetJ;
							break;
						}
						case 8:
						{
							element_shape_functions = &ShapeFunctions::Quadrilateral8;
							element_shape_functions_global = &ShapeFunctions::Quadrilateral8Global;
							facet_shape_functions = &ShapeFunctions::FacetLinear3;
							element_det_j = &ShapeFunctions::Quadrilateral8DetJ;
							break;
						}
						case 9:
						{
							element_shape_functions = &ShapeFunctions::Quadrilateral9;
							element_shape_functions_global = &ShapeFunctions::Quadrilateral9Global;
							facet_shape_functions = &ShapeFunctions::FacetLinear3;
							element_det_j = &ShapeFunctions::Quadrilateral9DetJ;
							break;
						}
						default:
						{
							Throw(ShapeFunctions::exception_invalid_shape);
						}
					}
					break;
				}
				case shape_tetrahedron:
				{
					switch (mesh.nodes_per_element)
					{
						case 4:
						{
							element_shape_functions = &ShapeFunctions::Tetrahedron4;
							element_shape_functions_global = &ShapeFunctions::Tetrahedron4Global;
							facet_shape_functions = &ShapeFunctions::FacetTriangle3;
							element_det_j = &ShapeFunctions::Tetrahedron4DetJ;
							break;
						}
						case 10:
						{
							element_shape_functions = &ShapeFunctions::Tetrahedron10;
							element_shape_functions_global = &ShapeFunctions::Tetrahedron10Global;
							facet_shape_functions = &ShapeFunctions::FacetTriangle6;
							element_det_j = &ShapeFunctions::Tetrahedron10DetJ;
							break;
						}
						default:
						{
							Throw(ShapeFunctions::exception_invalid_shape);
						}
					}
					break;
				}
				case shape_hexahedron:
				{
					switch (mesh.nodes_per_element)
					{
						case 8:
						{
							element_shape_functions = &ShapeFunctions::Hexahedron8;
							element_shape_functions_global = &ShapeFunctions::Hexahedron8Global;
							facet_shape_functions = &ShapeFunctions::FacetQuadrilateral4;
							element_det_j = &ShapeFunctions::Hexahedron8DetJ;
							break;
						}
						case 20:
						{
							element_shape_functions = &ShapeFunctions::Hexahedron20;
							element_shape_functions_global = &ShapeFunctions::Hexahedron20Global;
							facet_shape_functions = &ShapeFunctions::FacetQuadrilateral8;
							element_det_j = &ShapeFunctions::Hexahedron20DetJ;
							break;
						}
						case 27:
						{
							element_shape_functions = &ShapeFunctions::Hexahedron27;
							element_shape_functions_global = &ShapeFunctions::Hexahedron27Global;
							facet_shape_functions = &ShapeFunctions::FacetQuadrilateral9;
							element_det_j = &ShapeFunctions::Hexahedron27DetJ;
							break;
						}
						default:
						{
							Throw(ShapeFunctions::exception_invalid_shape);
						}
					}
					break;
				}
				default:
				{
					Throw(ShapeFunctions::exception_invalid_shape);
				}
			}
		}


		void ElementShapeFunctions(const int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			(this->*element_shape_functions)(element_id, point, N, dN, det_J);
		}


		bool ElementShapeFunctionsGlobal(int element_id, const Vector<T>& coordinate, Vector<T>& N) const throw()
		{
			return (this->*element_shape_functions_global)(element_id, coordinate, N);
		}


		void FacetShapeFunctions(const Vector<int>& facet_nodes, const Vector<T>& point, Vector<T>& N, T& det_J) const throw()
		{
			(this->*facet_shape_functions)(facet_nodes, point, N, det_J);
		}


		T ElementDetJ(const int element_id, const Vector<T>& point) const throw()
		{
			return (this->*element_det_j)(element_id, point);
		}


	protected:


		T Triangle3DetJ(int element_id, const Vector<T>&) const throw()
		{
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];

			// Local coordinates derivative
			T dxdrho = x2 - x1;
			T dxdeta = x3 - x1;
			T dydrho = y2 - y1;
			T dydeta = y3 - y1;

			// Jacobian determinant
			return dxdrho*dydeta - dxdeta*dydrho;
		}


		T Triangle6DetJ(int element_id, const Vector<T>& point) const throw()
		{
			register T rho = point.entry[1];
			register T eta = point.entry[2];
			register T xi = 1 - rho - eta;

			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			int node5 = mesh.connectivity.entry[element_id][5];
			int node6 = mesh.connectivity.entry[element_id][6];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];
			T x5 = nodes.coordinate.entry[node5][1];
			T y5 = nodes.coordinate.entry[node5][2];
			T x6 = nodes.coordinate.entry[node6][1];
			T y6 = nodes.coordinate.entry[node6][2];

			// Local coordinates derivative
			T dxdrho = -4*eta*x6 + 4*eta*x5 + 4*(xi - rho)*x4 + (4*rho - 1)*x2 - (4*xi - 1)*x1;
			T dxdeta =  4*(xi - eta)*x6 + 4*rho*x5 - 4*rho*x4 + (4*eta - 1)*x3 - (4*xi - 1)*x1;
			T dydrho = -4*eta*y6 + 4*eta*y5 + 4*(xi - rho)*y4 + (4*rho - 1)*y2 - (4*xi - 1)*y1;
			T dydeta =  4*(xi - eta)*y6 + 4*rho*y5 - 4*rho*y4 + (4*eta - 1)*y3 - (4*xi - 1)*y1;

			// Jacobian determinant
			return dxdrho*dydeta - dxdeta*dydrho;
		}

		
		T Quadrilateral4DetJ(int element_id, const Vector<T>& point) const throw()
		{
			register T rho = point.entry[1];
			register T eta = point.entry[2];

			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];

			// Local coordinates derivative
			T dxdrho = T(-0.25)*((1 + eta)*x4 - (1 + eta)*x3 - (1 - eta)*x2 + (1 - eta)*x1);
			T dxdeta = T( 0.25)*((1 - rho)*x4 + (1 + rho)*x3 - (1 + rho)*x2 - (1 - rho)*x1);
			T dydrho = T(-0.25)*((1 + eta)*y4 - (1 + eta)*y3 - (1 - eta)*y2 + (1 - eta)*y1);
			T dydeta = T( 0.25)*((1 - rho)*y4 + (1 + rho)*y3 - (1 + rho)*y2 - (1 - rho)*y1);

			// Jacobian determinant
			return dxdrho*dydeta - dxdeta*dydrho;
		}


		T Quadrilateral8DetJ(int element_id, const Vector<T>& point) const throw()
		{
			register T rho = point.entry[1];
			register T eta = point.entry[2];

			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			int node5 = mesh.connectivity.entry[element_id][5];
			int node6 = mesh.connectivity.entry[element_id][6];
			int node7 = mesh.connectivity.entry[element_id][7];
			int node8 = mesh.connectivity.entry[element_id][8];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];
			T x5 = nodes.coordinate.entry[node5][1];
			T y5 = nodes.coordinate.entry[node5][2];
			T x6 = nodes.coordinate.entry[node6][1];
			T y6 = nodes.coordinate.entry[node6][2];
			T x7 = nodes.coordinate.entry[node7][1];
			T y7 = nodes.coordinate.entry[node7][2];
			T x8 = nodes.coordinate.entry[node8][1];
			T y8 = nodes.coordinate.entry[node8][2];

			// Local coordinates derivative
			T dxdrho = T(0.25)*((2*eta*eta - 2)*x8 + (-4*eta - 4)*rho*x7 + (2 - 2*eta*eta)*x6 + (4*eta - 4)*rho*x5 + ((2*eta + 2)*rho - eta*eta - eta)*x4 + ((2*eta + 2)*rho + eta*eta + eta)*x3 + ((2 - 2*eta)*rho + eta*eta - eta)*x2 + ((2 - 2*eta)*rho - eta*eta + eta)*x1);
			T dxdeta = T(0.25)*((4*eta*rho - 4*eta)*x8 + (2 - 2*rho*rho)*x7 + (-4*eta*rho - 4*eta)*x6 + (2*rho*rho - 2)*x5 + (rho*rho + (-2*eta - 1)*rho + 2*eta)*x4 + (rho*rho + (2*eta + 1)*rho + 2*eta)*x3 + (-rho*rho + (2*eta - 1)*rho + 2*eta)*x2 + (-rho*rho + (1 - 2*eta)*rho + 2*eta)*x1);
			T dydrho = T(0.25)*((2*eta*eta - 2)*y8 + (-4*eta - 4)*rho*y7 + (2 - 2*eta*eta)*y6 + (4*eta - 4)*rho*y5 + ((2*eta + 2)*rho - eta*eta - eta)*y4 + ((2*eta + 2)*rho + eta*eta + eta)*y3 + ((2 - 2*eta)*rho + eta*eta - eta)*y2 + ((2 - 2*eta)*rho - eta*eta + eta)*y1);
			T dydeta = T(0.25)*((4*eta*rho - 4*eta)*y8 + (2 - 2*rho*rho)*y7 + (-4*eta*rho - 4*eta)*y6 + (2*rho*rho - 2)*y5 + (rho*rho + (-2*eta - 1)*rho + 2*eta)*y4 + (rho*rho + (2*eta + 1)*rho + 2*eta)*y3 + (-rho*rho + (2*eta - 1)*rho + 2*eta)*y2 + (-rho*rho + (1 - 2*eta)*rho + 2*eta)*y1);

			// Jacobian determinant
			return dxdrho*dydeta - dxdeta*dydrho;
		}


		T Quadrilateral9DetJ(int element_id, const Vector<T>& point) const throw()
		{
			register T rho = point.entry[1];
			register T eta = point.entry[2];

			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			int node5 = mesh.connectivity.entry[element_id][5];
			int node6 = mesh.connectivity.entry[element_id][6];
			int node7 = mesh.connectivity.entry[element_id][7];
			int node8 = mesh.connectivity.entry[element_id][8];
			int node9 = mesh.connectivity.entry[element_id][9];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];
			T x5 = nodes.coordinate.entry[node5][1];
			T y5 = nodes.coordinate.entry[node5][2];
			T x6 = nodes.coordinate.entry[node6][1];
			T y6 = nodes.coordinate.entry[node6][2];
			T x7 = nodes.coordinate.entry[node7][1];
			T y7 = nodes.coordinate.entry[node7][2];
			T x8 = nodes.coordinate.entry[node8][1];
			T y8 = nodes.coordinate.entry[node8][2];
			T x9 = nodes.coordinate.entry[node9][1];
			T y9 = nodes.coordinate.entry[node9][2];

			// Local coordinates derivative
			T dxdrho = T(0.25)*(8*(eta*eta - 1)*rho*x9 + (4*(1 - eta*eta)*rho + 2*eta*eta - 2)*x8 - 4*( eta + eta*eta)*rho*x7 + (4*(1 - eta*eta)*rho - 2*eta*eta + 2)*x6 + 4*(eta - eta*eta)*rho*x5 + (2*(eta*eta + eta)*rho - eta*eta - eta)*x4 + (2*(eta*eta + eta)*rho + eta*eta + eta)*x3 + (2*(eta*eta - eta)*rho + eta*eta - eta)*x2 + (2*(eta*eta - eta)*rho - eta*eta + eta)*x1);
			T dxdeta = T(0.25)*(8*(eta*rho*rho - eta)*x9 + 4*(eta*rho - eta*rho*rho)*x8 + (( - 4*eta - 2)*rho*rho + 4*eta + 2)*x7 - 4*(eta*rho*rho + eta*rho)*x6 + ((2 - 4*eta)*rho*rho + 4*eta - 2)*x5 + ((2*eta + 1)*rho*rho + ( - 2*eta - 1)*rho)*x4 + ((2*eta + 1)*rho*rho + (2*eta + 1)*rho)*x3 + ((2*eta - 1)*rho*rho + (2*eta - 1)*rho)*x2 + ((2*eta - 1)*rho*rho + (1 - 2*eta)*rho)*x1);
			T dydrho = T(0.25)*(8*(eta*eta - 1)*rho*y9 + (4*(1 - eta*eta)*rho + 2*eta*eta - 2)*y8 - 4*( eta + eta*eta)*rho*y7 + (4*(1 - eta*eta)*rho - 2*eta*eta + 2)*y6 + 4*(eta - eta*eta)*rho*y5 + (2*(eta*eta + eta)*rho - eta*eta - eta)*y4 + (2*(eta*eta + eta)*rho + eta*eta + eta)*y3 + (2*(eta*eta - eta)*rho + eta*eta - eta)*y2 + (2*(eta*eta - eta)*rho - eta*eta + eta)*y1);
			T dydeta = T(0.25)*(8*(eta*rho*rho - eta)*y9 + 4*(eta*rho - eta*rho*rho)*y8 + (( - 4*eta - 2)*rho*rho + 4*eta + 2)*y7 - 4*(eta*rho*rho + eta*rho)*y6 + ((2 - 4*eta)*rho*rho + 4*eta - 2)*y5 + ((2*eta + 1)*rho*rho + ( - 2*eta - 1)*rho)*y4 + ((2*eta + 1)*rho*rho + (2*eta + 1)*rho)*y3 + ((2*eta - 1)*rho*rho + (2*eta - 1)*rho)*y2 + ((2*eta - 1)*rho*rho + (1 - 2*eta)*rho)*y1);

			// Jacobian determinant
			return dxdrho*dydeta - dxdeta*dydrho;
		}


		T Tetrahedron4DetJ(int element_id, const Vector<T>&) const throw()
		{
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T z1 = nodes.coordinate.entry[node1][3];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T z2 = nodes.coordinate.entry[node2][3];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T z3 = nodes.coordinate.entry[node3][3];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];
			T z4 = nodes.coordinate.entry[node4][3];

			// Local coordinates derivative
			T dxdrho  = x2 - x1;
			T dxdeta  = x3 - x1;
			T dxdzeta = x4 - x1;
			T dydrho  = y2 - y1;
			T dydeta  = y3 - y1;
			T dydzeta = y4 - y1;
			T dzdrho  = z2 - z1;
			T dzdeta  = z3 - z1;
			T dzdzeta = z4 - z1;

			// Jacobian determinant
			return dxdrho*(dydeta*dzdzeta - dydzeta*dzdeta) - dydrho*(dxdeta*dzdzeta - dxdzeta*dzdeta) + (dxdeta*dydzeta - dxdzeta*dydeta)*dzdrho;
		}


		T Tetrahedron10DetJ(int element_id, const Vector<T>& point) const throw()
		{
			register T rho  = point.entry[1];
			register T eta  = point.entry[2];
			register T zeta = point.entry[3];

			int node1  = mesh.connectivity.entry[element_id][1];
			int node2  = mesh.connectivity.entry[element_id][2];
			int node3  = mesh.connectivity.entry[element_id][3];
			int node4  = mesh.connectivity.entry[element_id][4];
			int node5  = mesh.connectivity.entry[element_id][5];
			int node6  = mesh.connectivity.entry[element_id][6];
			int node7  = mesh.connectivity.entry[element_id][7];
			int node8  = mesh.connectivity.entry[element_id][8];
			int node9  = mesh.connectivity.entry[element_id][9];
			int node10 = mesh.connectivity.entry[element_id][10];
			T x1  = nodes.coordinate.entry[node1][1];
			T y1  = nodes.coordinate.entry[node1][2];
			T z1  = nodes.coordinate.entry[node1][3];
			T x2  = nodes.coordinate.entry[node2][1];
			T y2  = nodes.coordinate.entry[node2][2];
			T z2  = nodes.coordinate.entry[node2][3];
			T x3  = nodes.coordinate.entry[node3][1];
			T y3  = nodes.coordinate.entry[node3][2];
			T z3  = nodes.coordinate.entry[node3][3];
			T x4  = nodes.coordinate.entry[node4][1];
			T y4  = nodes.coordinate.entry[node4][2];
			T z4  = nodes.coordinate.entry[node4][3];
			T x5  = nodes.coordinate.entry[node5][1];
			T y5  = nodes.coordinate.entry[node5][2];
			T z5  = nodes.coordinate.entry[node5][3];
			T x6  = nodes.coordinate.entry[node6][1];
			T y6  = nodes.coordinate.entry[node6][2];
			T z6  = nodes.coordinate.entry[node6][3];
			T x7  = nodes.coordinate.entry[node7][1];
			T y7  = nodes.coordinate.entry[node7][2];
			T z7  = nodes.coordinate.entry[node7][3];
			T x8  = nodes.coordinate.entry[node8][1];
			T y8  = nodes.coordinate.entry[node8][2];
			T z8  = nodes.coordinate.entry[node8][3];
			T x9  = nodes.coordinate.entry[node9][1];
			T y9  = nodes.coordinate.entry[node9][2];
			T z9  = nodes.coordinate.entry[node9][3];
			T x10 = nodes.coordinate.entry[node10][1];
			T y10 = nodes.coordinate.entry[node10][2];
			T z10 = nodes.coordinate.entry[node10][3];

			// Local coordinates derivative
			T dxdrho  = x1*(4*(zeta + rho + eta) - 3) - 4*x5*(zeta + 2*rho + eta - 1) + 4*x9*zeta - 4*x8*zeta - 4*eta*x7 + 4*eta*x6 + (4*rho - 1)*x2;
			T dxdeta  = x1*(4*(zeta + rho + eta) - 3) - 4*x7*(zeta + rho + 2*eta - 1) - 4*x8*zeta + 4*x10*zeta + 4*rho*x6 - 4*rho*x5 + (4*eta - 1)*x3;
			T dxdzeta = x1*(4*(zeta + rho + eta) - 3) + x4*(4*zeta - 1) - 4*x8*(2*zeta + rho + eta - 1) + 4*rho*x9 - 4*eta*x7 - 4*rho*x5 + 4*eta*x10;
			T dydrho  = y1*(4*(zeta + rho + eta) - 3) - 4*y5*(zeta + 2*rho + eta - 1) + 4*y9*zeta - 4*y8*zeta - 4*eta*y7 + 4*eta*y6 + (4*rho - 1)*y2;
			T dydeta  = y1*(4*(zeta + rho + eta) - 3) - 4*y7*(zeta + rho + 2*eta - 1) - 4*y8*zeta + 4*y10*zeta + 4*rho*y6 - 4*rho*y5 + (4*eta - 1)*y3;
			T dydzeta = y1*(4*(zeta + rho + eta) - 3) + y4*(4*zeta - 1) - 4*y8*(2*zeta + rho + eta - 1) + 4*rho*y9 - 4*eta*y7 - 4*rho*y5 + 4*eta*y10;
			T dzdrho  = z1*(4*(zeta + rho + eta) - 3) - 4*z5*(zeta + 2*rho + eta - 1) + 4*z9*zeta - 4*z8*zeta - 4*eta*z7 + 4*eta*z6 + (4*rho - 1)*z2;
			T dzdeta  = z1*(4*(zeta + rho + eta) - 3) - 4*z7*(zeta + rho + 2*eta - 1) - 4*z8*zeta + 4*z10*zeta + 4*rho*z6 - 4*rho*z5 + (4*eta - 1)*z3;
			T dzdzeta = z1*(4*(zeta + rho + eta) - 3) + z4*(4*zeta - 1) - 4*z8*(2*zeta + rho + eta - 1) + 4*rho*z9 - 4*eta*z7 - 4*rho*z5 + 4*eta*z10;

			// Jacobian determinant
			return dxdrho*(dydeta*dzdzeta - dydzeta*dzdeta) - dydrho*(dxdeta*dzdzeta - dxdzeta*dzdeta) + (dxdeta*dydzeta - dxdzeta*dydeta)*dzdrho;
		}


		T Hexahedron8DetJ(int element_id, const Vector<T>& point) const throw()
		{
			register T rho  = point.entry[1];
			register T eta  = point.entry[2];
			register T zeta = point.entry[3];

			int node1  = mesh.connectivity.entry[element_id][1];
			int node2  = mesh.connectivity.entry[element_id][2];
			int node3  = mesh.connectivity.entry[element_id][3];
			int node4  = mesh.connectivity.entry[element_id][4];
			int node5  = mesh.connectivity.entry[element_id][5];
			int node6  = mesh.connectivity.entry[element_id][6];
			int node7  = mesh.connectivity.entry[element_id][7];
			int node8  = mesh.connectivity.entry[element_id][8];
			T x1  = nodes.coordinate.entry[node1][1];
			T y1  = nodes.coordinate.entry[node1][2];
			T z1  = nodes.coordinate.entry[node1][3];
			T x2  = nodes.coordinate.entry[node2][1];
			T y2  = nodes.coordinate.entry[node2][2];
			T z2  = nodes.coordinate.entry[node2][3];
			T x3  = nodes.coordinate.entry[node3][1];
			T y3  = nodes.coordinate.entry[node3][2];
			T z3  = nodes.coordinate.entry[node3][3];
			T x4  = nodes.coordinate.entry[node4][1];
			T y4  = nodes.coordinate.entry[node4][2];
			T z4  = nodes.coordinate.entry[node4][3];
			T x5  = nodes.coordinate.entry[node5][1];
			T y5  = nodes.coordinate.entry[node5][2];
			T z5  = nodes.coordinate.entry[node5][3];
			T x6  = nodes.coordinate.entry[node6][1];
			T y6  = nodes.coordinate.entry[node6][2];
			T z6  = nodes.coordinate.entry[node6][3];
			T x7  = nodes.coordinate.entry[node7][1];
			T y7  = nodes.coordinate.entry[node7][2];
			T z7  = nodes.coordinate.entry[node7][3];
			T x8  = nodes.coordinate.entry[node8][1];
			T y8  = nodes.coordinate.entry[node8][2];
			T z8  = nodes.coordinate.entry[node8][3];

			// Local coordinates derivative
			T dxdrho  = (-(eta + 1)*(zeta + 1)*x8 + (eta + 1)*(zeta + 1)*x7 - (eta - 1)*(zeta + 1)*x6 + (eta - 1)*(zeta + 1)*x5 + (eta + 1)*(zeta - 1)*x4 - (eta + 1)*(zeta - 1)*x3 + (eta - 1)*(zeta - 1)*x2 - (eta - 1)*(zeta - 1)*x1)*T(0.125);
			T dxdeta  = (-(rho - 1)*(zeta + 1)*x8 + (rho + 1)*(zeta + 1)*x7 - (rho + 1)*(zeta + 1)*x6 + (rho - 1)*(zeta + 1)*x5 + (rho - 1)*(zeta - 1)*x4 - (rho + 1)*(zeta - 1)*x3 + (rho + 1)*(zeta - 1)*x2 - (rho - 1)*(zeta - 1)*x1)*T(0.125);
			T dxdzeta = (-(eta + 1)*(rho - 1)*x8 + (eta + 1)*(rho + 1)*x7 - (eta - 1)*(rho + 1)*x6 + (eta - 1)*(rho - 1)*x5 + (eta + 1)*(rho - 1)*x4 - (eta + 1)*(rho + 1)*x3 + (eta - 1)*(rho + 1)*x2 - (eta - 1)*(rho - 1)*x1)*T(0.125);
			T dydrho  = (-(eta + 1)*(zeta + 1)*y8 + (eta + 1)*(zeta + 1)*y7 - (eta - 1)*(zeta + 1)*y6 + (eta - 1)*(zeta + 1)*y5 + (eta + 1)*(zeta - 1)*y4 - (eta + 1)*(zeta - 1)*y3 + (eta - 1)*(zeta - 1)*y2 - (eta - 1)*(zeta - 1)*y1)*T(0.125);
			T dydeta  = (-(rho - 1)*(zeta + 1)*y8 + (rho + 1)*(zeta + 1)*y7 - (rho + 1)*(zeta + 1)*y6 + (rho - 1)*(zeta + 1)*y5 + (rho - 1)*(zeta - 1)*y4 - (rho + 1)*(zeta - 1)*y3 + (rho + 1)*(zeta - 1)*y2 - (rho - 1)*(zeta - 1)*y1)*T(0.125);
			T dydzeta = (-(eta + 1)*(rho - 1)*y8 + (eta + 1)*(rho + 1)*y7 - (eta - 1)*(rho + 1)*y6 + (eta - 1)*(rho - 1)*y5 + (eta + 1)*(rho - 1)*y4 - (eta + 1)*(rho + 1)*y3 + (eta - 1)*(rho + 1)*y2 - (eta - 1)*(rho - 1)*y1)*T(0.125);
			T dzdrho  = (-(eta + 1)*(zeta + 1)*z8 + (eta + 1)*(zeta + 1)*z7 - (eta - 1)*(zeta + 1)*z6 + (eta - 1)*(zeta + 1)*z5 + (eta + 1)*(zeta - 1)*z4 - (eta + 1)*(zeta - 1)*z3 + (eta - 1)*(zeta - 1)*z2 - (eta - 1)*(zeta - 1)*z1)*T(0.125);
			T dzdeta  = (-(rho - 1)*(zeta + 1)*z8 + (rho + 1)*(zeta + 1)*z7 - (rho + 1)*(zeta + 1)*z6 + (rho - 1)*(zeta + 1)*z5 + (rho - 1)*(zeta - 1)*z4 - (rho + 1)*(zeta - 1)*z3 + (rho + 1)*(zeta - 1)*z2 - (rho - 1)*(zeta - 1)*z1)*T(0.125);
			T dzdzeta = (-(eta + 1)*(rho - 1)*z8 + (eta + 1)*(rho + 1)*z7 - (eta - 1)*(rho + 1)*z6 + (eta - 1)*(rho - 1)*z5 + (eta + 1)*(rho - 1)*z4 - (eta + 1)*(rho + 1)*z3 + (eta - 1)*(rho + 1)*z2 - (eta - 1)*(rho - 1)*z1)*T(0.125);

			// Jacobian determinant
			return dxdrho*(dydeta*dzdzeta - dydzeta*dzdeta) - dydrho*(dxdeta*dzdzeta - dxdzeta*dzdeta) + (dxdeta*dydzeta - dxdzeta*dydeta)*dzdrho;
		}


		T Hexahedron20DetJ(int element_id, const Vector<T>& point) const throw()
		{
			register T rho  = point.entry[1];
			register T eta  = point.entry[2];
			register T zeta = point.entry[3];

			int node1  = mesh.connectivity.entry[element_id][1];
			int node2  = mesh.connectivity.entry[element_id][2];
			int node3  = mesh.connectivity.entry[element_id][3];
			int node4  = mesh.connectivity.entry[element_id][4];
			int node5  = mesh.connectivity.entry[element_id][5];
			int node6  = mesh.connectivity.entry[element_id][6];
			int node7  = mesh.connectivity.entry[element_id][7];
			int node8  = mesh.connectivity.entry[element_id][8];
			int node9  = mesh.connectivity.entry[element_id][9];
			int node10 = mesh.connectivity.entry[element_id][10];
			int node11 = mesh.connectivity.entry[element_id][11];
			int node12 = mesh.connectivity.entry[element_id][12];
			int node13 = mesh.connectivity.entry[element_id][13];
			int node14 = mesh.connectivity.entry[element_id][14];
			int node15 = mesh.connectivity.entry[element_id][15];
			int node16 = mesh.connectivity.entry[element_id][16];
			int node17 = mesh.connectivity.entry[element_id][17];
			int node18 = mesh.connectivity.entry[element_id][18];
			int node19 = mesh.connectivity.entry[element_id][19];
			int node20 = mesh.connectivity.entry[element_id][20];
			T x1  = nodes.coordinate.entry[node1][1];
			T y1  = nodes.coordinate.entry[node1][2];
			T z1  = nodes.coordinate.entry[node1][3];
			T x2  = nodes.coordinate.entry[node2][1];
			T y2  = nodes.coordinate.entry[node2][2];
			T z2  = nodes.coordinate.entry[node2][3];
			T x3  = nodes.coordinate.entry[node3][1];
			T y3  = nodes.coordinate.entry[node3][2];
			T z3  = nodes.coordinate.entry[node3][3];
			T x4  = nodes.coordinate.entry[node4][1];
			T y4  = nodes.coordinate.entry[node4][2];
			T z4  = nodes.coordinate.entry[node4][3];
			T x5  = nodes.coordinate.entry[node5][1];
			T y5  = nodes.coordinate.entry[node5][2];
			T z5  = nodes.coordinate.entry[node5][3];
			T x6  = nodes.coordinate.entry[node6][1];
			T y6  = nodes.coordinate.entry[node6][2];
			T z6  = nodes.coordinate.entry[node6][3];
			T x7  = nodes.coordinate.entry[node7][1];
			T y7  = nodes.coordinate.entry[node7][2];
			T z7  = nodes.coordinate.entry[node7][3];
			T x8  = nodes.coordinate.entry[node8][1];
			T y8  = nodes.coordinate.entry[node8][2];
			T z8  = nodes.coordinate.entry[node8][3];
			T x9  = nodes.coordinate.entry[node9][1];
			T y9  = nodes.coordinate.entry[node9][2];
			T z9  = nodes.coordinate.entry[node9][3];
			T x10 = nodes.coordinate.entry[node10][1];
			T y10 = nodes.coordinate.entry[node10][2];
			T z10 = nodes.coordinate.entry[node10][3];
			T x11 = nodes.coordinate.entry[node11][1];
			T y11 = nodes.coordinate.entry[node11][2];
			T z11 = nodes.coordinate.entry[node11][3];
			T x12 = nodes.coordinate.entry[node12][1];
			T y12 = nodes.coordinate.entry[node12][2];
			T z12 = nodes.coordinate.entry[node12][3];
			T x13 = nodes.coordinate.entry[node13][1];
			T y13 = nodes.coordinate.entry[node13][2];
			T z13 = nodes.coordinate.entry[node13][3];
			T x14 = nodes.coordinate.entry[node14][1];
			T y14 = nodes.coordinate.entry[node14][2];
			T z14 = nodes.coordinate.entry[node14][3];
			T x15 = nodes.coordinate.entry[node15][1];
			T y15 = nodes.coordinate.entry[node15][2];
			T z15 = nodes.coordinate.entry[node15][3];
			T x16 = nodes.coordinate.entry[node16][1];
			T y16 = nodes.coordinate.entry[node16][2];
			T z16 = nodes.coordinate.entry[node16][3];
			T x17 = nodes.coordinate.entry[node17][1];
			T y17 = nodes.coordinate.entry[node17][2];
			T z17 = nodes.coordinate.entry[node17][3];
			T x18 = nodes.coordinate.entry[node18][1];
			T y18 = nodes.coordinate.entry[node18][2];
			T z18 = nodes.coordinate.entry[node18][3];
			T x19 = nodes.coordinate.entry[node19][1];
			T y19 = nodes.coordinate.entry[node19][2];
			T z19 = nodes.coordinate.entry[node19][3];
			T x20 = nodes.coordinate.entry[node20][1];
			T y20 = nodes.coordinate.entry[node20][2];
			T z20 = nodes.coordinate.entry[node20][3];

			// Local coordinates derivative
			T dxdrho  = ((eta - 1)*x1*(zeta - 1)*(zeta + 2*rho + eta + 1) + (eta + 1)*x7*(zeta + 1)*(zeta + 2*rho + eta - 1) - (eta + 1)*x4*(zeta - 1)*(zeta + 2*rho - eta + 1) - (eta - 1)*x6*(zeta + 1)*(zeta + 2*rho - eta - 1) - (eta - 1)*x2*(zeta - 1)*(zeta - 2*rho + eta + 1) - (eta + 1)*x8*(zeta + 1)*(zeta - 2*rho + eta - 1) + (eta + 1)*x3*(zeta - 1)*(zeta - 2*rho - eta + 1) + (eta - 1)*x5*(zeta + 1)*(zeta - 2*rho - eta - 1) + 2*(eta + 1)*x16*(zeta - 1)*(zeta + 1) - 2*(eta + 1)*x15*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*x14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*x13*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*(eta + 1)*x20*(zeta + 1) - 4*(eta + 1)*rho*x19*(zeta + 1) - 2*(eta - 1)*(eta + 1)*x18*(zeta + 1) + 4*(eta - 1)*rho*x17*(zeta + 1) - 4*(eta - 1)*rho*x9*(zeta - 1) - 2*(eta - 1)*(eta + 1)*x12*(zeta - 1) + 4*(eta + 1)*rho*x11*(zeta - 1) + 2*(eta - 1)*(eta + 1)*x10*(zeta - 1))*T(0.125);
			T dxdeta  = ((rho - 1)*x1*(zeta - 1)*(zeta + rho + 2*eta + 1) + (rho + 1)*x7*(zeta + 1)*(zeta + rho + 2*eta - 1) - (rho - 1)*x4*(zeta - 1)*(zeta + rho - 2*eta + 1) - (rho + 1)*x6*(zeta + 1)*(zeta + rho - 2*eta - 1) - (rho + 1)*x2*(zeta - 1)*(zeta - rho + 2*eta + 1) - (rho - 1)*x8*(zeta + 1)*(zeta - rho + 2*eta - 1) + (rho + 1)*x3*(zeta - 1)*(zeta - rho - 2*eta + 1) + (rho - 1)*x5*(zeta + 1)*(zeta - rho - 2*eta - 1) + 2*(rho - 1)*x16*(zeta - 1)*(zeta + 1) - 2*(rho + 1)*x15*(zeta - 1)*(zeta + 1) + 2*(rho + 1)*x14*(zeta - 1)*(zeta + 1) - 2*(rho - 1)*x13*(zeta - 1)*(zeta + 1) + 4*eta*(rho - 1)*x20*(zeta + 1) - 2*(rho - 1)*(rho + 1)*x19*(zeta + 1) - 4*eta*(rho + 1)*x18*(zeta + 1) + 2*(rho - 1)*(rho + 1)*x17*(zeta + 1) - 2*(rho - 1)*(rho + 1)*x9*(zeta - 1) - 4*eta*(rho - 1)*x12*(zeta - 1) + 2*(rho - 1)*(rho + 1)*x11*(zeta - 1) + 4*eta*(rho + 1)*x10*(zeta - 1))*T(0.125);
			T dxdzeta = ((eta - 1)*(rho - 1)*x1*(2*zeta + rho + eta + 1) + (eta + 1)*(rho + 1)*x7*(2*zeta + rho + eta - 1) - (eta + 1)*(rho - 1)*x4*(2*zeta + rho - eta + 1) - (eta - 1)*(rho + 1)*x6*(2*zeta + rho - eta - 1) - (eta - 1)*(rho + 1)*x2*(2*zeta - rho + eta + 1) - (eta + 1)*(rho - 1)*x8*(2*zeta - rho + eta - 1) + (eta + 1)*(rho + 1)*x3*(2*zeta - rho - eta + 1) + (eta - 1)*(rho - 1)*x5*(2*zeta - rho - eta - 1) + 4*(eta + 1)*(rho - 1)*x16*zeta - 4*(eta + 1)*(rho + 1)*x15*zeta + 4*(eta - 1)*(rho + 1)*x14*zeta - 4*(eta - 1)*(rho - 1)*x13*zeta - 2*(eta - 1)*(rho - 1)*(rho + 1)*x9 + 2*(eta - 1)*(eta + 1)*(rho - 1)*x20 - 2*(eta + 1)*(rho - 1)*(rho + 1)*x19 - 2*(eta - 1)*(eta + 1)*(rho + 1)*x18 + 2*(eta - 1)*(rho - 1)*(rho + 1)*x17 - 2*(eta - 1)*(eta + 1)*(rho - 1)*x12 + 2*(eta + 1)*(rho - 1)*(rho + 1)*x11 + 2*(eta - 1)*(eta + 1)*(rho + 1)*x10)*T(0.125);
			T dydrho  = ((eta - 1)*y1*(zeta - 1)*(zeta + 2*rho + eta + 1) + (eta + 1)*y7*(zeta + 1)*(zeta + 2*rho + eta - 1) - (eta + 1)*y4*(zeta - 1)*(zeta + 2*rho - eta + 1) - (eta - 1)*y6*(zeta + 1)*(zeta + 2*rho - eta - 1) - (eta - 1)*y2*(zeta - 1)*(zeta - 2*rho + eta + 1) - (eta + 1)*y8*(zeta + 1)*(zeta - 2*rho + eta - 1) + (eta + 1)*y3*(zeta - 1)*(zeta - 2*rho - eta + 1) + (eta - 1)*y5*(zeta + 1)*(zeta - 2*rho - eta - 1) + 2*(eta + 1)*y16*(zeta - 1)*(zeta + 1) - 2*(eta + 1)*y15*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*y14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*y13*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*(eta + 1)*y20*(zeta + 1) - 4*(eta + 1)*rho*y19*(zeta + 1) - 2*(eta - 1)*(eta + 1)*y18*(zeta + 1) + 4*(eta - 1)*rho*y17*(zeta + 1) - 4*(eta - 1)*rho*y9*(zeta - 1) - 2*(eta - 1)*(eta + 1)*y12*(zeta - 1) + 4*(eta + 1)*rho*y11*(zeta - 1) + 2*(eta - 1)*(eta + 1)*y10*(zeta - 1))*T(0.125);
			T dydeta  = ((rho - 1)*y1*(zeta - 1)*(zeta + rho + 2*eta + 1) + (rho + 1)*y7*(zeta + 1)*(zeta + rho + 2*eta - 1) - (rho - 1)*y4*(zeta - 1)*(zeta + rho - 2*eta + 1) - (rho + 1)*y6*(zeta + 1)*(zeta + rho - 2*eta - 1) - (rho + 1)*y2*(zeta - 1)*(zeta - rho + 2*eta + 1) - (rho - 1)*y8*(zeta + 1)*(zeta - rho + 2*eta - 1) + (rho + 1)*y3*(zeta - 1)*(zeta - rho - 2*eta + 1) + (rho - 1)*y5*(zeta + 1)*(zeta - rho - 2*eta - 1) + 2*(rho - 1)*y16*(zeta - 1)*(zeta + 1) - 2*(rho + 1)*y15*(zeta - 1)*(zeta + 1) + 2*(rho + 1)*y14*(zeta - 1)*(zeta + 1) - 2*(rho - 1)*y13*(zeta - 1)*(zeta + 1) + 4*eta*(rho - 1)*y20*(zeta + 1) - 2*(rho - 1)*(rho + 1)*y19*(zeta + 1) - 4*eta*(rho + 1)*y18*(zeta + 1) + 2*(rho - 1)*(rho + 1)*y17*(zeta + 1) - 2*(rho - 1)*(rho + 1)*y9*(zeta - 1) - 4*eta*(rho - 1)*y12*(zeta - 1) + 2*(rho - 1)*(rho + 1)*y11*(zeta - 1) + 4*eta*(rho + 1)*y10*(zeta - 1))*T(0.125);
			T dydzeta = ((eta - 1)*(rho - 1)*y1*(2*zeta + rho + eta + 1) + (eta + 1)*(rho + 1)*y7*(2*zeta + rho + eta - 1) - (eta + 1)*(rho - 1)*y4*(2*zeta + rho - eta + 1) - (eta - 1)*(rho + 1)*y6*(2*zeta + rho - eta - 1) - (eta - 1)*(rho + 1)*y2*(2*zeta - rho + eta + 1) - (eta + 1)*(rho - 1)*y8*(2*zeta - rho + eta - 1) + (eta + 1)*(rho + 1)*y3*(2*zeta - rho - eta + 1) + (eta - 1)*(rho - 1)*y5*(2*zeta - rho - eta - 1) + 4*(eta + 1)*(rho - 1)*y16*zeta - 4*(eta + 1)*(rho + 1)*y15*zeta + 4*(eta - 1)*(rho + 1)*y14*zeta - 4*(eta - 1)*(rho - 1)*y13*zeta - 2*(eta - 1)*(rho - 1)*(rho + 1)*y9 + 2*(eta - 1)*(eta + 1)*(rho - 1)*y20 - 2*(eta + 1)*(rho - 1)*(rho + 1)*y19 - 2*(eta - 1)*(eta + 1)*(rho + 1)*y18 + 2*(eta - 1)*(rho - 1)*(rho + 1)*y17 - 2*(eta - 1)*(eta + 1)*(rho - 1)*y12 + 2*(eta + 1)*(rho - 1)*(rho + 1)*y11 + 2*(eta - 1)*(eta + 1)*(rho + 1)*y10)*T(0.125);
			T dzdrho  = ((eta - 1)*z1*(zeta - 1)*(zeta + 2*rho + eta + 1) + (eta + 1)*z7*(zeta + 1)*(zeta + 2*rho + eta - 1) - (eta + 1)*z4*(zeta - 1)*(zeta + 2*rho - eta + 1) - (eta - 1)*z6*(zeta + 1)*(zeta + 2*rho - eta - 1) - (eta - 1)*z2*(zeta - 1)*(zeta - 2*rho + eta + 1) - (eta + 1)*z8*(zeta + 1)*(zeta - 2*rho + eta - 1) + (eta + 1)*z3*(zeta - 1)*(zeta - 2*rho - eta + 1) + (eta - 1)*z5*(zeta + 1)*(zeta - 2*rho - eta - 1) + 2*(eta + 1)*z16*(zeta - 1)*(zeta + 1) - 2*(eta + 1)*z15*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*z14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*z13*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*(eta + 1)*z20*(zeta + 1) - 4*(eta + 1)*rho*z19*(zeta + 1) - 2*(eta - 1)*(eta + 1)*z18*(zeta + 1) + 4*(eta - 1)*rho*z17*(zeta + 1) - 4*(eta - 1)*rho*z9*(zeta - 1) - 2*(eta - 1)*(eta + 1)*z12*(zeta - 1) + 4*(eta + 1)*rho*z11*(zeta - 1) + 2*(eta - 1)*(eta + 1)*z10*(zeta - 1))*T(0.125);
			T dzdeta  = ((rho - 1)*z1*(zeta - 1)*(zeta + rho + 2*eta + 1) + (rho + 1)*z7*(zeta + 1)*(zeta + rho + 2*eta - 1) - (rho - 1)*z4*(zeta - 1)*(zeta + rho - 2*eta + 1) - (rho + 1)*z6*(zeta + 1)*(zeta + rho - 2*eta - 1) - (rho + 1)*z2*(zeta - 1)*(zeta - rho + 2*eta + 1) - (rho - 1)*z8*(zeta + 1)*(zeta - rho + 2*eta - 1) + (rho + 1)*z3*(zeta - 1)*(zeta - rho - 2*eta + 1) + (rho - 1)*z5*(zeta + 1)*(zeta - rho - 2*eta - 1) + 2*(rho - 1)*z16*(zeta - 1)*(zeta + 1) - 2*(rho + 1)*z15*(zeta - 1)*(zeta + 1) + 2*(rho + 1)*z14*(zeta - 1)*(zeta + 1) - 2*(rho - 1)*z13*(zeta - 1)*(zeta + 1) + 4*eta*(rho - 1)*z20*(zeta + 1) - 2*(rho - 1)*(rho + 1)*z19*(zeta + 1) - 4*eta*(rho + 1)*z18*(zeta + 1) + 2*(rho - 1)*(rho + 1)*z17*(zeta + 1) - 2*(rho - 1)*(rho + 1)*z9*(zeta - 1) - 4*eta*(rho - 1)*z12*(zeta - 1) + 2*(rho - 1)*(rho + 1)*z11*(zeta - 1) + 4*eta*(rho + 1)*z10*(zeta - 1))*T(0.125);
			T dzdzeta = ((eta - 1)*(rho - 1)*z1*(2*zeta + rho + eta + 1) + (eta + 1)*(rho + 1)*z7*(2*zeta + rho + eta - 1) - (eta + 1)*(rho - 1)*z4*(2*zeta + rho - eta + 1) - (eta - 1)*(rho + 1)*z6*(2*zeta + rho - eta - 1) - (eta - 1)*(rho + 1)*z2*(2*zeta - rho + eta + 1) - (eta + 1)*(rho - 1)*z8*(2*zeta - rho + eta - 1) + (eta + 1)*(rho + 1)*z3*(2*zeta - rho - eta + 1) + (eta - 1)*(rho - 1)*z5*(2*zeta - rho - eta - 1) + 4*(eta + 1)*(rho - 1)*z16*zeta - 4*(eta + 1)*(rho + 1)*z15*zeta + 4*(eta - 1)*(rho + 1)*z14*zeta - 4*(eta - 1)*(rho - 1)*z13*zeta - 2*(eta - 1)*(rho - 1)*(rho + 1)*z9 + 2*(eta - 1)*(eta + 1)*(rho - 1)*z20 - 2*(eta + 1)*(rho - 1)*(rho + 1)*z19 - 2*(eta - 1)*(eta + 1)*(rho + 1)*z18 + 2*(eta - 1)*(rho - 1)*(rho + 1)*z17 - 2*(eta - 1)*(eta + 1)*(rho - 1)*z12 + 2*(eta + 1)*(rho - 1)*(rho + 1)*z11 + 2*(eta - 1)*(eta + 1)*(rho + 1)*z10)*T(0.125);

			// Jacobian determinant
			return dxdrho*(dydeta*dzdzeta - dydzeta*dzdeta) - dydrho*(dxdeta*dzdzeta - dxdzeta*dzdeta) + (dxdeta*dydzeta - dxdzeta*dydeta)*dzdrho;
		}


		T Hexahedron27DetJ(int element_id, const Vector<T>& point) const throw()
		{
			register T rho  = point.entry[1];
			register T eta  = point.entry[2];
			register T zeta = point.entry[3];

			int node1  = mesh.connectivity.entry[element_id][1];
			int node2  = mesh.connectivity.entry[element_id][2];
			int node3  = mesh.connectivity.entry[element_id][3];
			int node4  = mesh.connectivity.entry[element_id][4];
			int node5  = mesh.connectivity.entry[element_id][5];
			int node6  = mesh.connectivity.entry[element_id][6];
			int node7  = mesh.connectivity.entry[element_id][7];
			int node8  = mesh.connectivity.entry[element_id][8];
			int node9  = mesh.connectivity.entry[element_id][9];
			int node10 = mesh.connectivity.entry[element_id][10];
			int node11 = mesh.connectivity.entry[element_id][11];
			int node12 = mesh.connectivity.entry[element_id][12];
			int node13 = mesh.connectivity.entry[element_id][13];
			int node14 = mesh.connectivity.entry[element_id][14];
			int node15 = mesh.connectivity.entry[element_id][15];
			int node16 = mesh.connectivity.entry[element_id][16];
			int node17 = mesh.connectivity.entry[element_id][17];
			int node18 = mesh.connectivity.entry[element_id][18];
			int node19 = mesh.connectivity.entry[element_id][19];
			int node20 = mesh.connectivity.entry[element_id][20];
			int node21 = mesh.connectivity.entry[element_id][21];
			int node22 = mesh.connectivity.entry[element_id][22];
			int node23 = mesh.connectivity.entry[element_id][23];
			int node24 = mesh.connectivity.entry[element_id][24];
			int node25 = mesh.connectivity.entry[element_id][25];
			int node26 = mesh.connectivity.entry[element_id][26];
			int node27 = mesh.connectivity.entry[element_id][27];
			T x1  = nodes.coordinate.entry[node1][1];
			T y1  = nodes.coordinate.entry[node1][2];
			T z1  = nodes.coordinate.entry[node1][3];
			T x2  = nodes.coordinate.entry[node2][1];
			T y2  = nodes.coordinate.entry[node2][2];
			T z2  = nodes.coordinate.entry[node2][3];
			T x3  = nodes.coordinate.entry[node3][1];
			T y3  = nodes.coordinate.entry[node3][2];
			T z3  = nodes.coordinate.entry[node3][3];
			T x4  = nodes.coordinate.entry[node4][1];
			T y4  = nodes.coordinate.entry[node4][2];
			T z4  = nodes.coordinate.entry[node4][3];
			T x5  = nodes.coordinate.entry[node5][1];
			T y5  = nodes.coordinate.entry[node5][2];
			T z5  = nodes.coordinate.entry[node5][3];
			T x6  = nodes.coordinate.entry[node6][1];
			T y6  = nodes.coordinate.entry[node6][2];
			T z6  = nodes.coordinate.entry[node6][3];
			T x7  = nodes.coordinate.entry[node7][1];
			T y7  = nodes.coordinate.entry[node7][2];
			T z7  = nodes.coordinate.entry[node7][3];
			T x8  = nodes.coordinate.entry[node8][1];
			T y8  = nodes.coordinate.entry[node8][2];
			T z8  = nodes.coordinate.entry[node8][3];
			T x9  = nodes.coordinate.entry[node9][1];
			T y9  = nodes.coordinate.entry[node9][2];
			T z9  = nodes.coordinate.entry[node9][3];
			T x10 = nodes.coordinate.entry[node10][1];
			T y10 = nodes.coordinate.entry[node10][2];
			T z10 = nodes.coordinate.entry[node10][3];
			T x11 = nodes.coordinate.entry[node11][1];
			T y11 = nodes.coordinate.entry[node11][2];
			T z11 = nodes.coordinate.entry[node11][3];
			T x12 = nodes.coordinate.entry[node12][1];
			T y12 = nodes.coordinate.entry[node12][2];
			T z12 = nodes.coordinate.entry[node12][3];
			T x13 = nodes.coordinate.entry[node13][1];
			T y13 = nodes.coordinate.entry[node13][2];
			T z13 = nodes.coordinate.entry[node13][3];
			T x14 = nodes.coordinate.entry[node14][1];
			T y14 = nodes.coordinate.entry[node14][2];
			T z14 = nodes.coordinate.entry[node14][3];
			T x15 = nodes.coordinate.entry[node15][1];
			T y15 = nodes.coordinate.entry[node15][2];
			T z15 = nodes.coordinate.entry[node15][3];
			T x16 = nodes.coordinate.entry[node16][1];
			T y16 = nodes.coordinate.entry[node16][2];
			T z16 = nodes.coordinate.entry[node16][3];
			T x17 = nodes.coordinate.entry[node17][1];
			T y17 = nodes.coordinate.entry[node17][2];
			T z17 = nodes.coordinate.entry[node17][3];
			T x18 = nodes.coordinate.entry[node18][1];
			T y18 = nodes.coordinate.entry[node18][2];
			T z18 = nodes.coordinate.entry[node18][3];
			T x19 = nodes.coordinate.entry[node19][1];
			T y19 = nodes.coordinate.entry[node19][2];
			T z19 = nodes.coordinate.entry[node19][3];
			T x20 = nodes.coordinate.entry[node20][1];
			T y20 = nodes.coordinate.entry[node20][2];
			T z20 = nodes.coordinate.entry[node20][3];
			T x21 = nodes.coordinate.entry[node21][1];
			T y21 = nodes.coordinate.entry[node21][2];
			T z21 = nodes.coordinate.entry[node21][3];
			T x22 = nodes.coordinate.entry[node22][1];
			T y22 = nodes.coordinate.entry[node22][2];
			T z22 = nodes.coordinate.entry[node22][3];
			T x23 = nodes.coordinate.entry[node23][1];
			T y23 = nodes.coordinate.entry[node23][2];
			T z23 = nodes.coordinate.entry[node23][3];
			T x24 = nodes.coordinate.entry[node24][1];
			T y24 = nodes.coordinate.entry[node24][2];
			T z24 = nodes.coordinate.entry[node24][3];
			T x25 = nodes.coordinate.entry[node25][1];
			T y25 = nodes.coordinate.entry[node25][2];
			T z25 = nodes.coordinate.entry[node25][3];
			T x26 = nodes.coordinate.entry[node26][1];
			T y26 = nodes.coordinate.entry[node26][2];
			T z26 = nodes.coordinate.entry[node26][3];
			T x27 = nodes.coordinate.entry[node27][1];
			T y27 = nodes.coordinate.entry[node27][2];
			T z27 = nodes.coordinate.entry[node27][3];

			// Local coordinates derivative
			T dxdrho  = (eta*(eta + 1)*(2*rho - 1)*x8*zeta*(zeta + 1) + eta*(eta + 1)*(2*rho + 1)*x7*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho + 1)*x6*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho - 1)*x5*zeta*(zeta + 1) + 8*(eta - 1)*(eta + 1)*rho*x26*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*x20*zeta*(zeta + 1) - 4*eta*(eta + 1)*rho*x19*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*x18*zeta*(zeta + 1) - 4*(eta - 1)*eta*rho*x17*zeta*(zeta + 1) - 16*(eta - 1)*(eta + 1)*rho*x27*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho - 1)*x25*(zeta - 1)*(zeta + 1) + 8*eta*(eta + 1)*rho*x24*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho + 1)*x23*(zeta - 1)*(zeta + 1) + 8*(eta - 1)*eta*rho*x22*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho - 1)*x16*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho + 1)*x15*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho + 1)*x14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho - 1)*x13*(zeta - 1)*(zeta + 1) - 4*(eta - 1)*eta*rho*x9*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho - 1)*x4*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho + 1)*x3*(zeta - 1)*zeta + 8*(eta - 1)*(eta + 1)*rho*x21*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho + 1)*x2*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*x12*(zeta - 1)*zeta - 4*eta*(eta + 1)*rho*x11*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*x10*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho - 1)*x1*(zeta - 1)*zeta)*T(0.125);
			T dxdeta  = ((2*eta + 1)*(rho - 1)*rho*x8*zeta*(zeta + 1) + (2*eta + 1)*rho*(rho + 1)*x7*zeta*(zeta + 1) + (2*eta - 1)*rho*(rho + 1)*x6*zeta*(zeta + 1) + (2*eta - 1)*(rho - 1)*rho*x5*zeta*(zeta + 1) + 8*eta*(rho - 1)*(rho + 1)*x26*zeta*(zeta + 1) - 4*eta*(rho - 1)*rho*x20*zeta*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*x19*zeta*(zeta + 1) - 4*eta*rho*(rho + 1)*x18*zeta*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*x17*zeta*(zeta + 1) - 16*eta*(rho - 1)*(rho + 1)*x27*(zeta - 1)*(zeta + 1) + 8*eta*(rho - 1)*rho*x25*(zeta - 1)*(zeta + 1) + 4*(2*eta + 1)*(rho - 1)*(rho + 1)*x24*(zeta - 1)*(zeta + 1) + 8*eta*rho*(rho + 1)*x23*(zeta - 1)*(zeta + 1) + 4*(2*eta - 1)*(rho - 1)*(rho + 1)*x22*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*rho*x16*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*rho*(rho + 1)*x15*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*rho*(rho + 1)*x14*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*rho*x13*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*x9*(zeta - 1)*zeta + (2*eta + 1)*(rho - 1)*rho*x4*(zeta - 1)*zeta + (2*eta + 1)*rho*(rho + 1)*x3*(zeta - 1)*zeta + 8*eta*(rho - 1)*(rho + 1)*x21*(zeta - 1)*zeta + (2*eta - 1)*rho*(rho + 1)*x2*(zeta - 1)*zeta - 4*eta*(rho - 1)*rho*x12*(zeta - 1)*zeta - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*x11*(zeta - 1)*zeta - 4*eta*rho*(rho + 1)*x10*(zeta - 1)*zeta + (2*eta - 1)*(rho - 1)*rho*x1*(zeta - 1)*zeta)*T(0.125);
			T dxdzeta = (eta*(eta + 1)*(rho - 1)*rho*x8*(2*zeta + 1) + eta*(eta + 1)*rho*(rho + 1)*x7*(2*zeta + 1) + (eta - 1)*eta*rho*(rho + 1)*x6*(2*zeta + 1) + (eta - 1)*eta*(rho - 1)*rho*x5*(2*zeta + 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*x26*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*x20*(2*zeta + 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*x19*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*x18*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*x17*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*x9*(2*zeta - 1) + eta*(eta + 1)*(rho - 1)*rho*x4*(2*zeta - 1) + eta*(eta + 1)*rho*(rho + 1)*x3*(2*zeta - 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*x21*(2*zeta - 1) + (eta - 1)*eta*rho*(rho + 1)*x2*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*x12*(2*zeta - 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*x11*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*x10*(2*zeta - 1) + (eta - 1)*eta*(rho - 1)*rho*x1*(2*zeta - 1) - 16*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*x27*zeta + 8*(eta - 1)*(eta + 1)*(rho - 1)*rho*x25*zeta + 8*eta*(eta + 1)*(rho - 1)*(rho + 1)*x24*zeta + 8*(eta - 1)*(eta + 1)*rho*(rho + 1)*x23*zeta + 8*(eta - 1)*eta*(rho - 1)*(rho + 1)*x22*zeta - 4*eta*(eta + 1)*(rho - 1)*rho*x16*zeta - 4*eta*(eta + 1)*rho*(rho + 1)*x15*zeta - 4*(eta - 1)*eta*rho*(rho + 1)*x14*zeta - 4*(eta - 1)*eta*(rho - 1)*rho*x13*zeta)*T(0.125);
			T dydrho  = (eta*(eta + 1)*(2*rho - 1)*y8*zeta*(zeta + 1) + eta*(eta + 1)*(2*rho + 1)*y7*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho + 1)*y6*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho - 1)*y5*zeta*(zeta + 1) + 8*(eta - 1)*(eta + 1)*rho*y26*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*y20*zeta*(zeta + 1) - 4*eta*(eta + 1)*rho*y19*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*y18*zeta*(zeta + 1) - 4*(eta - 1)*eta*rho*y17*zeta*(zeta + 1) - 16*(eta - 1)*(eta + 1)*rho*y27*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho - 1)*y25*(zeta - 1)*(zeta + 1) + 8*eta*(eta + 1)*rho*y24*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho + 1)*y23*(zeta - 1)*(zeta + 1) + 8*(eta - 1)*eta*rho*y22*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho - 1)*y16*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho + 1)*y15*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho + 1)*y14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho - 1)*y13*(zeta - 1)*(zeta + 1) - 4*(eta - 1)*eta*rho*y9*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho - 1)*y4*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho + 1)*y3*(zeta - 1)*zeta + 8*(eta - 1)*(eta + 1)*rho*y21*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho + 1)*y2*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*y12*(zeta - 1)*zeta - 4*eta*(eta + 1)*rho*y11*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*y10*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho - 1)*y1*(zeta - 1)*zeta)*T(0.125);
			T dydeta  = ((2*eta + 1)*(rho - 1)*rho*y8*zeta*(zeta + 1) + (2*eta + 1)*rho*(rho + 1)*y7*zeta*(zeta + 1) + (2*eta - 1)*rho*(rho + 1)*y6*zeta*(zeta + 1) + (2*eta - 1)*(rho - 1)*rho*y5*zeta*(zeta + 1) + 8*eta*(rho - 1)*(rho + 1)*y26*zeta*(zeta + 1) - 4*eta*(rho - 1)*rho*y20*zeta*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*y19*zeta*(zeta + 1) - 4*eta*rho*(rho + 1)*y18*zeta*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*y17*zeta*(zeta + 1) - 16*eta*(rho - 1)*(rho + 1)*y27*(zeta - 1)*(zeta + 1) + 8*eta*(rho - 1)*rho*y25*(zeta - 1)*(zeta + 1) + 4*(2*eta + 1)*(rho - 1)*(rho + 1)*y24*(zeta - 1)*(zeta + 1) + 8*eta*rho*(rho + 1)*y23*(zeta - 1)*(zeta + 1) + 4*(2*eta - 1)*(rho - 1)*(rho + 1)*y22*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*rho*y16*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*rho*(rho + 1)*y15*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*rho*(rho + 1)*y14*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*rho*y13*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*y9*(zeta - 1)*zeta + (2*eta + 1)*(rho - 1)*rho*y4*(zeta - 1)*zeta + (2*eta + 1)*rho*(rho + 1)*y3*(zeta - 1)*zeta + 8*eta*(rho - 1)*(rho + 1)*y21*(zeta - 1)*zeta + (2*eta - 1)*rho*(rho + 1)*y2*(zeta - 1)*zeta - 4*eta*(rho - 1)*rho*y12*(zeta - 1)*zeta - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*y11*(zeta - 1)*zeta - 4*eta*rho*(rho + 1)*y10*(zeta - 1)*zeta + (2*eta - 1)*(rho - 1)*rho*y1*(zeta - 1)*zeta)*T(0.125);
			T dydzeta = (eta*(eta + 1)*(rho - 1)*rho*y8*(2*zeta + 1) + eta*(eta + 1)*rho*(rho + 1)*y7*(2*zeta + 1) + (eta - 1)*eta*rho*(rho + 1)*y6*(2*zeta + 1) + (eta - 1)*eta*(rho - 1)*rho*y5*(2*zeta + 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*y26*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*y20*(2*zeta + 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*y19*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*y18*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*y17*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*y9*(2*zeta - 1) + eta*(eta + 1)*(rho - 1)*rho*y4*(2*zeta - 1) + eta*(eta + 1)*rho*(rho + 1)*y3*(2*zeta - 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*y21*(2*zeta - 1) + (eta - 1)*eta*rho*(rho + 1)*y2*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*y12*(2*zeta - 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*y11*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*y10*(2*zeta - 1) + (eta - 1)*eta*(rho - 1)*rho*y1*(2*zeta - 1) - 16*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*y27*zeta + 8*(eta - 1)*(eta + 1)*(rho - 1)*rho*y25*zeta + 8*eta*(eta + 1)*(rho - 1)*(rho + 1)*y24*zeta + 8*(eta - 1)*(eta + 1)*rho*(rho + 1)*y23*zeta + 8*(eta - 1)*eta*(rho - 1)*(rho + 1)*y22*zeta - 4*eta*(eta + 1)*(rho - 1)*rho*y16*zeta - 4*eta*(eta + 1)*rho*(rho + 1)*y15*zeta - 4*(eta - 1)*eta*rho*(rho + 1)*y14*zeta - 4*(eta - 1)*eta*(rho - 1)*rho*y13*zeta)*T(0.125);
			T dzdrho  = (eta*(eta + 1)*(2*rho - 1)*z8*zeta*(zeta + 1) + eta*(eta + 1)*(2*rho + 1)*z7*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho + 1)*z6*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho - 1)*z5*zeta*(zeta + 1) + 8*(eta - 1)*(eta + 1)*rho*z26*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*z20*zeta*(zeta + 1) - 4*eta*(eta + 1)*rho*z19*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*z18*zeta*(zeta + 1) - 4*(eta - 1)*eta*rho*z17*zeta*(zeta + 1) - 16*(eta - 1)*(eta + 1)*rho*z27*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho - 1)*z25*(zeta - 1)*(zeta + 1) + 8*eta*(eta + 1)*rho*z24*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho + 1)*z23*(zeta - 1)*(zeta + 1) + 8*(eta - 1)*eta*rho*z22*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho - 1)*z16*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho + 1)*z15*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho + 1)*z14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho - 1)*z13*(zeta - 1)*(zeta + 1) - 4*(eta - 1)*eta*rho*z9*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho - 1)*z4*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho + 1)*z3*(zeta - 1)*zeta + 8*(eta - 1)*(eta + 1)*rho*z21*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho + 1)*z2*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*z12*(zeta - 1)*zeta - 4*eta*(eta + 1)*rho*z11*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*z10*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho - 1)*z1*(zeta - 1)*zeta)*T(0.125);
			T dzdeta  = ((2*eta + 1)*(rho - 1)*rho*z8*zeta*(zeta + 1) + (2*eta + 1)*rho*(rho + 1)*z7*zeta*(zeta + 1) + (2*eta - 1)*rho*(rho + 1)*z6*zeta*(zeta + 1) + (2*eta - 1)*(rho - 1)*rho*z5*zeta*(zeta + 1) + 8*eta*(rho - 1)*(rho + 1)*z26*zeta*(zeta + 1) - 4*eta*(rho - 1)*rho*z20*zeta*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*z19*zeta*(zeta + 1) - 4*eta*rho*(rho + 1)*z18*zeta*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*z17*zeta*(zeta + 1) - 16*eta*(rho - 1)*(rho + 1)*z27*(zeta - 1)*(zeta + 1) + 8*eta*(rho - 1)*rho*z25*(zeta - 1)*(zeta + 1) + 4*(2*eta + 1)*(rho - 1)*(rho + 1)*z24*(zeta - 1)*(zeta + 1) + 8*eta*rho*(rho + 1)*z23*(zeta - 1)*(zeta + 1) + 4*(2*eta - 1)*(rho - 1)*(rho + 1)*z22*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*rho*z16*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*rho*(rho + 1)*z15*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*rho*(rho + 1)*z14*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*rho*z13*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*z9*(zeta - 1)*zeta + (2*eta + 1)*(rho - 1)*rho*z4*(zeta - 1)*zeta + (2*eta + 1)*rho*(rho + 1)*z3*(zeta - 1)*zeta + 8*eta*(rho - 1)*(rho + 1)*z21*(zeta - 1)*zeta + (2*eta - 1)*rho*(rho + 1)*z2*(zeta - 1)*zeta - 4*eta*(rho - 1)*rho*z12*(zeta - 1)*zeta - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*z11*(zeta - 1)*zeta - 4*eta*rho*(rho + 1)*z10*(zeta - 1)*zeta + (2*eta - 1)*(rho - 1)*rho*z1*(zeta - 1)*zeta)*T(0.125);
			T dzdzeta = (eta*(eta + 1)*(rho - 1)*rho*z8*(2*zeta + 1) + eta*(eta + 1)*rho*(rho + 1)*z7*(2*zeta + 1) + (eta - 1)*eta*rho*(rho + 1)*z6*(2*zeta + 1) + (eta - 1)*eta*(rho - 1)*rho*z5*(2*zeta + 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*z26*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*z20*(2*zeta + 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*z19*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*z18*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*z17*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*z9*(2*zeta - 1) + eta*(eta + 1)*(rho - 1)*rho*z4*(2*zeta - 1) + eta*(eta + 1)*rho*(rho + 1)*z3*(2*zeta - 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*z21*(2*zeta - 1) + (eta - 1)*eta*rho*(rho + 1)*z2*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*z12*(2*zeta - 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*z11*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*z10*(2*zeta - 1) + (eta - 1)*eta*(rho - 1)*rho*z1*(2*zeta - 1) - 16*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*z27*zeta + 8*(eta - 1)*(eta + 1)*(rho - 1)*rho*z25*zeta + 8*eta*(eta + 1)*(rho - 1)*(rho + 1)*z24*zeta + 8*(eta - 1)*(eta + 1)*rho*(rho + 1)*z23*zeta + 8*(eta - 1)*eta*(rho - 1)*(rho + 1)*z22*zeta - 4*eta*(eta + 1)*(rho - 1)*rho*z16*zeta - 4*eta*(eta + 1)*rho*(rho + 1)*z15*zeta - 4*(eta - 1)*eta*rho*(rho + 1)*z14*zeta - 4*(eta - 1)*eta*(rho - 1)*rho*z13*zeta)*T(0.125);

			// Jacobian determinant
			return dxdrho*(dydeta*dzdzeta - dydzeta*dzdeta) - dydrho*(dxdeta*dzdzeta - dxdzeta*dzdeta) + (dxdeta*dydzeta - dxdzeta*dydeta)*dzdrho;
		}


		void Triangle3(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			register T rho = point.entry[1];
			register T eta = point.entry[2];

			// Shape functions
			N.entry[1] = 1 - rho - eta;
			N.entry[2] = rho;
			N.entry[3] = eta;

			// Local coordinates derivative
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];

			T dxdrho = x2 - x1;
			T dxdeta = x3 - x1;
			T dydrho = y2 - y1;
			T dydeta = y3 - y1;

			// Jacobian determinant
			det_J = dxdrho*dydeta - dxdeta*dydrho;

			// Global coordinates derivative
			T dN1drho = -1;
			T dN1deta = -1;
			T dN2drho =  1;
			T dN2deta =  0;
			T dN3drho =  0;
			T dN3deta =  1;

			dN.entry[1][1] = (dydeta*dN1drho - dydrho*dN1deta)/det_J;
			dN.entry[1][2] = (dxdrho*dN1deta - dxdeta*dN1drho)/det_J;
			dN.entry[2][1] = (dydeta*dN2drho - dydrho*dN2deta)/det_J;
			dN.entry[2][2] = (dxdrho*dN2deta - dxdeta*dN2drho)/det_J;
			dN.entry[3][1] = (dydeta*dN3drho - dydrho*dN3deta)/det_J;
			dN.entry[3][2] = (dxdrho*dN3deta - dxdeta*dN3drho)/det_J;
		}


		void Triangle6(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			register T rho = point.entry[1];
			register T eta = point.entry[2];
			register T xi = 1 - rho - eta;

			// Shape functions
			N.entry[1] = (2*xi - 1)*xi;;
			N.entry[2] = rho*(2*rho - 1);
			N.entry[3] = eta*(2*eta - 1);
			N.entry[4] = 4*rho*xi;
			N.entry[5] = 4*rho*eta;
			N.entry[6] = 4*eta*xi;

			// Local coordinates derivative
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			int node5 = mesh.connectivity.entry[element_id][5];
			int node6 = mesh.connectivity.entry[element_id][6];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];
			T x5 = nodes.coordinate.entry[node5][1];
			T y5 = nodes.coordinate.entry[node5][2];
			T x6 = nodes.coordinate.entry[node6][1];
			T y6 = nodes.coordinate.entry[node6][2];

			T dxdrho = -4*eta*x6 + 4*eta*x5 + 4*(xi - rho)*x4 + (4*rho - 1)*x2 - (4*xi - 1)*x1;
			T dxdeta =  4*(xi - eta)*x6 + 4*rho*x5 - 4*rho*x4 + (4*eta - 1)*x3 - (4*xi - 1)*x1;
			T dydrho = -4*eta*y6 + 4*eta*y5 + 4*(xi - rho)*y4 + (4*rho - 1)*y2 - (4*xi - 1)*y1;
			T dydeta =  4*(xi - eta)*y6 + 4*rho*y5 - 4*rho*y4 + (4*eta - 1)*y3 - (4*xi - 1)*y1;

			// Jacobian determinant
			det_J = dxdrho*dydeta - dxdeta*dydrho;

			// Global coordinates derivative
			T dN1drho =  1 - 4*xi;
			T dN1deta =  1 - 4*xi;
			T dN2drho =  4*rho - 1;
			T dN2deta =  0;
			T dN3drho =  0;
			T dN3deta =  4*eta - 1;
			T dN4drho =  4*xi - 4*rho;
			T dN4deta = -4*rho;
			T dN5drho =  4*eta;
			T dN5deta =  4*rho;
			T dN6drho = -4*eta;
			T dN6deta =  4*xi - 4*eta;

			dN.entry[1][1] = (dydeta*dN1drho - dydrho*dN1deta)/det_J;
			dN.entry[1][2] = (dxdrho*dN1deta - dxdeta*dN1drho)/det_J;
			dN.entry[2][1] = (dydeta*dN2drho - dydrho*dN2deta)/det_J;
			dN.entry[2][2] = (dxdrho*dN2deta - dxdeta*dN2drho)/det_J;
			dN.entry[3][1] = (dydeta*dN3drho - dydrho*dN3deta)/det_J;
			dN.entry[3][2] = (dxdrho*dN3deta - dxdeta*dN3drho)/det_J;
			dN.entry[4][1] = (dydeta*dN4drho - dydrho*dN4deta)/det_J;
			dN.entry[4][2] = (dxdrho*dN4deta - dxdeta*dN4drho)/det_J;
			dN.entry[5][1] = (dydeta*dN5drho - dydrho*dN5deta)/det_J;
			dN.entry[5][2] = (dxdrho*dN5deta - dxdeta*dN5drho)/det_J;
			dN.entry[6][1] = (dydeta*dN6drho - dydrho*dN6deta)/det_J;
			dN.entry[6][2] = (dxdrho*dN6deta - dxdeta*dN6drho)/det_J;
		}


		void Quadrilateral4(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			register T rho = point.entry[1];
			register T eta = point.entry[2];

			// Shape functions
			N.entry[1] = T(0.25)*(1 - rho)*(1 - eta);
			N.entry[2] = T(0.25)*(1 + rho)*(1 - eta);
			N.entry[3] = T(0.25)*(1 + rho)*(1 + eta);
			N.entry[4] = T(0.25)*(1 - rho)*(1 + eta);

			// Local coordinates derivative
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];

			T dxdrho = T(-0.25)*((1 + eta)*x4 - (1 + eta)*x3 - (1 - eta)*x2 + (1 - eta)*x1);
			T dxdeta = T( 0.25)*((1 - rho)*x4 + (1 + rho)*x3 - (1 + rho)*x2 - (1 - rho)*x1);
			T dydrho = T(-0.25)*((1 + eta)*y4 - (1 + eta)*y3 - (1 - eta)*y2 + (1 - eta)*y1);
			T dydeta = T( 0.25)*((1 - rho)*y4 + (1 + rho)*y3 - (1 + rho)*y2 - (1 - rho)*y1);

			// Jacobian determinant
			det_J = dxdrho*dydeta - dxdeta*dydrho;

			// Global coordinates derivative
			T dN1drho = T(-0.25)*(1 - eta);
			T dN1deta = T(-0.25)*(1 - rho);
			T dN2drho = T( 0.25)*(1 - eta);
			T dN2deta = T(-0.25)*(1 + rho);
			T dN3drho = T( 0.25)*(1 + eta);
			T dN3deta = T( 0.25)*(1 + rho);
			T dN4drho = T(-0.25)*(1 + eta);
			T dN4deta = T( 0.25)*(1 - rho);

			dN.entry[1][1] = (dydeta*dN1drho - dydrho*dN1deta)/det_J;
			dN.entry[1][2] = (dxdrho*dN1deta - dxdeta*dN1drho)/det_J;
			dN.entry[2][1] = (dydeta*dN2drho - dydrho*dN2deta)/det_J;
			dN.entry[2][2] = (dxdrho*dN2deta - dxdeta*dN2drho)/det_J;
			dN.entry[3][1] = (dydeta*dN3drho - dydrho*dN3deta)/det_J;
			dN.entry[3][2] = (dxdrho*dN3deta - dxdeta*dN3drho)/det_J;
			dN.entry[4][1] = (dydeta*dN4drho - dydrho*dN4deta)/det_J;
			dN.entry[4][2] = (dxdrho*dN4deta - dxdeta*dN4drho)/det_J;
		}


		void Quadrilateral8(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			register T rho = point.entry[1];
			register T eta = point.entry[2];

			// Shape functions
			N.entry[1] = T(0.25)*(1 - rho)*(1 - eta)*(-1 - rho - eta);
			N.entry[2] = T(0.25)*(1 + rho)*(1 - eta)*(-1 + rho - eta);
			N.entry[3] = T(0.25)*(1 + rho)*(1 + eta)*(-1 + rho + eta);
			N.entry[4] = T(0.25)*(1 - rho)*(1 + eta)*(-1 - rho + eta);
			N.entry[5] = T(0.5)*(1 - rho*rho)*(1 - eta);
			N.entry[6] = T(0.5)*(1 - eta*eta)*(1 + rho);
			N.entry[7] = T(0.5)*(1 - rho*rho)*(1 + eta);
			N.entry[8] = T(0.5)*(1 - eta*eta)*(1 - rho);

			// Local coordinates derivative
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			int node5 = mesh.connectivity.entry[element_id][5];
			int node6 = mesh.connectivity.entry[element_id][6];
			int node7 = mesh.connectivity.entry[element_id][7];
			int node8 = mesh.connectivity.entry[element_id][8];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];
			T x5 = nodes.coordinate.entry[node5][1];
			T y5 = nodes.coordinate.entry[node5][2];
			T x6 = nodes.coordinate.entry[node6][1];
			T y6 = nodes.coordinate.entry[node6][2];
			T x7 = nodes.coordinate.entry[node7][1];
			T y7 = nodes.coordinate.entry[node7][2];
			T x8 = nodes.coordinate.entry[node8][1];
			T y8 = nodes.coordinate.entry[node8][2];

			T dxdrho = T(0.25)*((2*eta*eta - 2)*x8 + (-4*eta - 4)*rho*x7 + (2 - 2*eta*eta)*x6 + (4*eta - 4)*rho*x5 + ((2*eta + 2)*rho - eta*eta - eta)*x4 + ((2*eta + 2)*rho + eta*eta + eta)*x3 + ((2 - 2*eta)*rho + eta*eta - eta)*x2 + ((2 - 2*eta)*rho - eta*eta + eta)*x1);
			T dxdeta = T(0.25)*((4*eta*rho - 4*eta)*x8 + (2 - 2*rho*rho)*x7 + (-4*eta*rho - 4*eta)*x6 + (2*rho*rho - 2)*x5 + (rho*rho + (-2*eta - 1)*rho + 2*eta)*x4 + (rho*rho + (2*eta + 1)*rho + 2*eta)*x3 + (-rho*rho + (2*eta - 1)*rho + 2*eta)*x2 + (-rho*rho + (1 - 2*eta)*rho + 2*eta)*x1);
			T dydrho = T(0.25)*((2*eta*eta - 2)*y8 + (-4*eta - 4)*rho*y7 + (2 - 2*eta*eta)*y6 + (4*eta - 4)*rho*y5 + ((2*eta + 2)*rho - eta*eta - eta)*y4 + ((2*eta + 2)*rho + eta*eta + eta)*y3 + ((2 - 2*eta)*rho + eta*eta - eta)*y2 + ((2 - 2*eta)*rho - eta*eta + eta)*y1);
			T dydeta = T(0.25)*((4*eta*rho - 4*eta)*y8 + (2 - 2*rho*rho)*y7 + (-4*eta*rho - 4*eta)*y6 + (2*rho*rho - 2)*y5 + (rho*rho + (-2*eta - 1)*rho + 2*eta)*y4 + (rho*rho + (2*eta + 1)*rho + 2*eta)*y3 + (-rho*rho + (2*eta - 1)*rho + 2*eta)*y2 + (-rho*rho + (1 - 2*eta)*rho + 2*eta)*y1);

			// Jacobian determinant
			det_J = dxdrho*dydeta - dxdeta*dydrho;

			// Global coordinates derivative
			T dN1drho = T(-0.25)*((eta - 1)*(2*rho + eta));
			T dN1deta = T(-0.25)*((rho - 1)*(rho + 2*eta));
			T dN2drho = T(-0.25)*((eta - 1)*(2*rho - eta));
			T dN2deta = T(-0.25)*((rho + 1)*(rho - 2*eta));
			T dN3drho = T( 0.25)*((eta + 1)*(2*rho + eta));
			T dN3deta = T( 0.25)*((rho + 1)*(rho + 2*eta));
			T dN4drho = T( 0.25)*((eta + 1)*(2*rho - eta));
			T dN4deta = T( 0.25)*((rho - 1)*(rho - 2*eta));
			T dN5drho =  (eta - 1)*rho;
			T dN5deta = T( 0.5)*((rho - 1)*(rho + 1));
			T dN6drho = T(-0.5)*((eta - 1)*(eta + 1));
			T dN6deta = -eta*(rho + 1);
			T dN7drho = -(eta + 1)*rho;
			T dN7deta = T(-0.5)*((rho - 1)*(rho + 1));
			T dN8drho = T( 0.5)*((eta - 1)*(eta + 1));
			T dN8deta =  eta*(rho - 1);

			dN.entry[1][1] = (dydeta*dN1drho - dydrho*dN1deta)/det_J;
			dN.entry[1][2] = (dxdrho*dN1deta - dxdeta*dN1drho)/det_J;
			dN.entry[2][1] = (dydeta*dN2drho - dydrho*dN2deta)/det_J;
			dN.entry[2][2] = (dxdrho*dN2deta - dxdeta*dN2drho)/det_J;
			dN.entry[3][1] = (dydeta*dN3drho - dydrho*dN3deta)/det_J;
			dN.entry[3][2] = (dxdrho*dN3deta - dxdeta*dN3drho)/det_J;
			dN.entry[4][1] = (dydeta*dN4drho - dydrho*dN4deta)/det_J;
			dN.entry[4][2] = (dxdrho*dN4deta - dxdeta*dN4drho)/det_J;
			dN.entry[5][1] = (dydeta*dN5drho - dydrho*dN5deta)/det_J;
			dN.entry[5][2] = (dxdrho*dN5deta - dxdeta*dN5drho)/det_J;
			dN.entry[6][1] = (dydeta*dN6drho - dydrho*dN6deta)/det_J;
			dN.entry[6][2] = (dxdrho*dN6deta - dxdeta*dN6drho)/det_J;
			dN.entry[7][1] = (dydeta*dN7drho - dydrho*dN7deta)/det_J;
			dN.entry[7][2] = (dxdrho*dN7deta - dxdeta*dN7drho)/det_J;
			dN.entry[8][1] = (dydeta*dN8drho - dydrho*dN8deta)/det_J;
			dN.entry[8][2] = (dxdrho*dN8deta - dxdeta*dN8drho)/det_J;
		}


		void Quadrilateral9(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			register T rho = point.entry[1];
			register T eta = point.entry[2];

			// Shape functions
			N.entry[1] = T( 0.25)*(1 - rho)*rho*(1 - eta)*eta;
			N.entry[2] = T(-0.25)*(1 + rho)*rho*(1 - eta)*eta;
			N.entry[3] = T( 0.25)*(1 + rho)*rho*(1 + eta)*eta;
			N.entry[4] = T(-0.25)*(1 - rho)*rho*(1 + eta)*eta;
			N.entry[5] = T(-0.5)*(1 - rho*rho)*(1 - eta)*eta;
			N.entry[6] = T( 0.5)*(1 - eta*eta)*(1 + rho)*rho;
			N.entry[7] = T( 0.5)*(1 - rho*rho)*(1 + eta)*eta;
			N.entry[8] = T(-0.5)*(1 - eta*eta)*(1 - rho)*rho;
			N.entry[9] = (1 - rho*rho)*(1 - eta*eta);

			// Local coordinates derivative
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			int node5 = mesh.connectivity.entry[element_id][5];
			int node6 = mesh.connectivity.entry[element_id][6];
			int node7 = mesh.connectivity.entry[element_id][7];
			int node8 = mesh.connectivity.entry[element_id][8];
			int node9 = mesh.connectivity.entry[element_id][9];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];
			T x5 = nodes.coordinate.entry[node5][1];
			T y5 = nodes.coordinate.entry[node5][2];
			T x6 = nodes.coordinate.entry[node6][1];
			T y6 = nodes.coordinate.entry[node6][2];
			T x7 = nodes.coordinate.entry[node7][1];
			T y7 = nodes.coordinate.entry[node7][2];
			T x8 = nodes.coordinate.entry[node8][1];
			T y8 = nodes.coordinate.entry[node8][2];
			T x9 = nodes.coordinate.entry[node9][1];
			T y9 = nodes.coordinate.entry[node9][2];

			T dxdrho = T(0.25)*(8*(eta*eta - 1)*rho*x9 + (4*(1 - eta*eta)*rho + 2*eta*eta - 2)*x8 - 4*( eta + eta*eta)*rho*x7 + (4*(1 - eta*eta)*rho - 2*eta*eta + 2)*x6 + 4*(eta - eta*eta)*rho*x5 + (2*(eta*eta + eta)*rho - eta*eta - eta)*x4 + (2*(eta*eta + eta)*rho + eta*eta + eta)*x3 + (2*(eta*eta - eta)*rho + eta*eta - eta)*x2 + (2*(eta*eta - eta)*rho - eta*eta + eta)*x1);
			T dxdeta = T(0.25)*(8*(eta*rho*rho - eta)*x9 + 4*(eta*rho - eta*rho*rho)*x8 + (( - 4*eta - 2)*rho*rho + 4*eta + 2)*x7 - 4*(eta*rho*rho + eta*rho)*x6 + ((2 - 4*eta)*rho*rho + 4*eta - 2)*x5 + ((2*eta + 1)*rho*rho + ( - 2*eta - 1)*rho)*x4 + ((2*eta + 1)*rho*rho + (2*eta + 1)*rho)*x3 + ((2*eta - 1)*rho*rho + (2*eta - 1)*rho)*x2 + ((2*eta - 1)*rho*rho + (1 - 2*eta)*rho)*x1);
			T dydrho = T(0.25)*(8*(eta*eta - 1)*rho*y9 + (4*(1 - eta*eta)*rho + 2*eta*eta - 2)*y8 - 4*( eta + eta*eta)*rho*y7 + (4*(1 - eta*eta)*rho - 2*eta*eta + 2)*y6 + 4*(eta - eta*eta)*rho*y5 + (2*(eta*eta + eta)*rho - eta*eta - eta)*y4 + (2*(eta*eta + eta)*rho + eta*eta + eta)*y3 + (2*(eta*eta - eta)*rho + eta*eta - eta)*y2 + (2*(eta*eta - eta)*rho - eta*eta + eta)*y1);
			T dydeta = T(0.25)*(8*(eta*rho*rho - eta)*y9 + 4*(eta*rho - eta*rho*rho)*y8 + (( - 4*eta - 2)*rho*rho + 4*eta + 2)*y7 - 4*(eta*rho*rho + eta*rho)*y6 + ((2 - 4*eta)*rho*rho + 4*eta - 2)*y5 + ((2*eta + 1)*rho*rho + ( - 2*eta - 1)*rho)*y4 + ((2*eta + 1)*rho*rho + (2*eta + 1)*rho)*y3 + ((2*eta - 1)*rho*rho + (2*eta - 1)*rho)*y2 + ((2*eta - 1)*rho*rho + (1 - 2*eta)*rho)*y1);

			// Jacobian determinant
			det_J = dxdrho*dydeta - dxdeta*dydrho;

			// Global coordinates derivative
			T dN1drho = T( 0.25)*((eta - 1)*eta*(2*rho - 1));
			T dN1deta = T( 0.25)*((2*eta - 1)*(rho - 1)*rho);
			T dN2drho = T( 0.25)*((eta - 1)*eta*(2*rho + 1));
			T dN2deta = T( 0.25)*((2*eta - 1)*rho*(rho + 1));
			T dN3drho = T( 0.25)*(eta*(eta + 1)*(2*rho + 1));
			T dN3deta = T( 0.25)*((2*eta + 1)*rho*(rho + 1));
			T dN4drho = T( 0.25)*(eta*(eta + 1)*(2*rho - 1));
			T dN4deta = T( 0.25)*((2*eta + 1)*(rho - 1)*rho);
			T dN5drho = -(eta - 1)*eta*rho;
			T dN5deta = T(-0.5)*((2*eta - 1)*(rho - 1)*(rho + 1));
			T dN6drho = T(-0.5)*((eta - 1)*(eta + 1)*(2*rho + 1));
			T dN6deta = -eta*rho*(rho + 1);
			T dN7drho = -eta*(eta + 1)*rho;
			T dN7deta = T(-0.5)*((2*eta + 1)*(rho - 1)*(rho + 1));
			T dN8drho = T(-0.5)*((eta - 1)*(eta + 1)*(2*rho - 1));
			T dN8deta = -eta*(rho - 1)*rho;
			T dN9drho =  2*(eta - 1)*(eta + 1)*rho;
			T dN9deta =  2*eta*(rho - 1)*(rho + 1);

			dN.entry[1][1] = (dydeta*dN1drho - dydrho*dN1deta)/det_J;
			dN.entry[1][2] = (dxdrho*dN1deta - dxdeta*dN1drho)/det_J;
			dN.entry[2][1] = (dydeta*dN2drho - dydrho*dN2deta)/det_J;
			dN.entry[2][2] = (dxdrho*dN2deta - dxdeta*dN2drho)/det_J;
			dN.entry[3][1] = (dydeta*dN3drho - dydrho*dN3deta)/det_J;
			dN.entry[3][2] = (dxdrho*dN3deta - dxdeta*dN3drho)/det_J;
			dN.entry[4][1] = (dydeta*dN4drho - dydrho*dN4deta)/det_J;
			dN.entry[4][2] = (dxdrho*dN4deta - dxdeta*dN4drho)/det_J;
			dN.entry[5][1] = (dydeta*dN5drho - dydrho*dN5deta)/det_J;
			dN.entry[5][2] = (dxdrho*dN5deta - dxdeta*dN5drho)/det_J;
			dN.entry[6][1] = (dydeta*dN6drho - dydrho*dN6deta)/det_J;
			dN.entry[6][2] = (dxdrho*dN6deta - dxdeta*dN6drho)/det_J;
			dN.entry[7][1] = (dydeta*dN7drho - dydrho*dN7deta)/det_J;
			dN.entry[7][2] = (dxdrho*dN7deta - dxdeta*dN7drho)/det_J;
			dN.entry[8][1] = (dydeta*dN8drho - dydrho*dN8deta)/det_J;
			dN.entry[8][2] = (dxdrho*dN8deta - dxdeta*dN8drho)/det_J;
			dN.entry[8][1] = (dydeta*dN8drho - dydrho*dN8deta)/det_J;
			dN.entry[8][2] = (dxdrho*dN8deta - dxdeta*dN8drho)/det_J;
			dN.entry[9][1] = (dydeta*dN9drho - dydrho*dN9deta)/det_J;
			dN.entry[9][2] = (dxdrho*dN9deta - dxdeta*dN9drho)/det_J;
		}


		void Tetrahedron4(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			register T rho  = point.entry[1];
			register T eta  = point.entry[2];
			register T zeta = point.entry[3];

			// Shape functions
			N.entry[1] = 1 - rho - eta - zeta;
			N.entry[2] = rho;
			N.entry[3] = eta;
			N.entry[4] = zeta;

			// Local coordinates derivative
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T z1 = nodes.coordinate.entry[node1][3];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T z2 = nodes.coordinate.entry[node2][3];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T z3 = nodes.coordinate.entry[node3][3];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];
			T z4 = nodes.coordinate.entry[node4][3];

			T dxdrho  = x2 - x1;
			T dxdeta  = x3 - x1;
			T dxdzeta = x4 - x1;
			T dydrho  = y2 - y1;
			T dydeta  = y3 - y1;
			T dydzeta = y4 - y1;
			T dzdrho  = z2 - z1;
			T dzdeta  = z3 - z1;
			T dzdzeta = z4 - z1;

			// Jacobian determinant
			det_J = dxdrho*(dydeta*dzdzeta - dydzeta*dzdeta) - dydrho*(dxdeta*dzdzeta - dxdzeta*dzdeta) + (dxdeta*dydzeta - dxdzeta*dydeta)*dzdrho;

			// Inverse of the Jacobian
			T invJ11 = (dydeta*dzdzeta - dydzeta*dzdeta)/det_J;
			T invJ12 = (dydzeta*dzdrho - dydrho*dzdzeta)/det_J;
			T invJ13 = (dydrho*dzdeta - dydeta*dzdrho)/det_J;
			T invJ21 = (dxdzeta*dzdeta - dxdeta*dzdzeta)/det_J;
			T invJ22 = (dxdrho*dzdzeta - dxdzeta*dzdrho)/det_J;
			T invJ23 = (dxdeta*dzdrho - dxdrho*dzdeta)/det_J;
			T invJ31 = (dxdeta*dydzeta - dxdzeta*dydeta)/det_J;
			T invJ32 = (dxdzeta*dydrho - dxdrho*dydzeta)/det_J;
			T invJ33 = (dxdrho*dydeta - dxdeta*dydrho)/det_J;

			// Global coordinates derivative
			T dN1drho  = -1;
			T dN1deta  = -1;
			T dN1dzeta = -1;
			T dN2drho  =  1;
			T dN2deta  =  0;
			T dN2dzeta =  0;
			T dN3drho  =  0;
			T dN3deta  =  1;
			T dN3dzeta =  0;
			T dN4drho  =  0;
			T dN4deta  =  0;
			T dN4dzeta =  1;

			dN.entry[1][1]  = invJ11*dN1drho + invJ12*dN1deta + invJ13*dN1dzeta;
			dN.entry[1][2]  = invJ21*dN1drho + invJ22*dN1deta + invJ23*dN1dzeta;
			dN.entry[1][3]  = invJ31*dN1drho + invJ32*dN1deta + invJ33*dN1dzeta;
			dN.entry[2][1]  = invJ11*dN2drho + invJ12*dN2deta + invJ13*dN2dzeta;
			dN.entry[2][2]  = invJ21*dN2drho + invJ22*dN2deta + invJ23*dN2dzeta;
			dN.entry[2][3]  = invJ31*dN2drho + invJ32*dN2deta + invJ33*dN2dzeta;
			dN.entry[3][1]  = invJ11*dN3drho + invJ12*dN3deta + invJ13*dN3dzeta;
			dN.entry[3][2]  = invJ21*dN3drho + invJ22*dN3deta + invJ23*dN3dzeta;
			dN.entry[3][3]  = invJ31*dN3drho + invJ32*dN3deta + invJ33*dN3dzeta;
			dN.entry[4][1]  = invJ11*dN4drho + invJ12*dN4deta + invJ13*dN4dzeta;
			dN.entry[4][2]  = invJ21*dN4drho + invJ22*dN4deta + invJ23*dN4dzeta;
			dN.entry[4][3]  = invJ31*dN4drho + invJ32*dN4deta + invJ33*dN4dzeta;
		}


		void Tetrahedron10(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			register T rho  = point.entry[1];
			register T eta  = point.entry[2];
			register T zeta = point.entry[3];

			// Shape functions
			N.entry[1]  = (1 - rho - eta - zeta)*(1 - 2*rho - 2*eta - 2*zeta);
			N.entry[2]  = rho*(2*rho - 1);
			N.entry[3]  = eta*(2*eta - 1);
			N.entry[4]  = zeta*(2*zeta - 1);
			N.entry[5]  = 4*rho*(1 - rho - eta - zeta);
			N.entry[6]  = 4*rho*eta;
			N.entry[7]  = 4*eta*(1 - rho - eta - zeta);
			N.entry[8]  = 4*zeta*(1 - rho - eta - zeta);
			N.entry[9]  = 4*rho*zeta;
			N.entry[10] = 4*eta*zeta;

			// Local coordinates derivative
			int node1  = mesh.connectivity.entry[element_id][1];
			int node2  = mesh.connectivity.entry[element_id][2];
			int node3  = mesh.connectivity.entry[element_id][3];
			int node4  = mesh.connectivity.entry[element_id][4];
			int node5  = mesh.connectivity.entry[element_id][5];
			int node6  = mesh.connectivity.entry[element_id][6];
			int node7  = mesh.connectivity.entry[element_id][7];
			int node8  = mesh.connectivity.entry[element_id][8];
			int node9  = mesh.connectivity.entry[element_id][9];
			int node10 = mesh.connectivity.entry[element_id][10];
			T x1  = nodes.coordinate.entry[node1][1];
			T y1  = nodes.coordinate.entry[node1][2];
			T z1  = nodes.coordinate.entry[node1][3];
			T x2  = nodes.coordinate.entry[node2][1];
			T y2  = nodes.coordinate.entry[node2][2];
			T z2  = nodes.coordinate.entry[node2][3];
			T x3  = nodes.coordinate.entry[node3][1];
			T y3  = nodes.coordinate.entry[node3][2];
			T z3  = nodes.coordinate.entry[node3][3];
			T x4  = nodes.coordinate.entry[node4][1];
			T y4  = nodes.coordinate.entry[node4][2];
			T z4  = nodes.coordinate.entry[node4][3];
			T x5  = nodes.coordinate.entry[node5][1];
			T y5  = nodes.coordinate.entry[node5][2];
			T z5  = nodes.coordinate.entry[node5][3];
			T x6  = nodes.coordinate.entry[node6][1];
			T y6  = nodes.coordinate.entry[node6][2];
			T z6  = nodes.coordinate.entry[node6][3];
			T x7  = nodes.coordinate.entry[node7][1];
			T y7  = nodes.coordinate.entry[node7][2];
			T z7  = nodes.coordinate.entry[node7][3];
			T x8  = nodes.coordinate.entry[node8][1];
			T y8  = nodes.coordinate.entry[node8][2];
			T z8  = nodes.coordinate.entry[node8][3];
			T x9  = nodes.coordinate.entry[node9][1];
			T y9  = nodes.coordinate.entry[node9][2];
			T z9  = nodes.coordinate.entry[node9][3];
			T x10 = nodes.coordinate.entry[node10][1];
			T y10 = nodes.coordinate.entry[node10][2];
			T z10 = nodes.coordinate.entry[node10][3];

			T dxdrho  = x1*(4*(zeta + rho + eta) - 3) - 4*x5*(zeta + 2*rho + eta - 1) + 4*x9*zeta - 4*x8*zeta - 4*eta*x7 + 4*eta*x6 + (4*rho - 1)*x2;
			T dxdeta  = x1*(4*(zeta + rho + eta) - 3) - 4*x7*(zeta + rho + 2*eta - 1) - 4*x8*zeta + 4*x10*zeta + 4*rho*x6 - 4*rho*x5 + (4*eta - 1)*x3;
			T dxdzeta = x1*(4*(zeta + rho + eta) - 3) + x4*(4*zeta - 1) - 4*x8*(2*zeta + rho + eta - 1) + 4*rho*x9 - 4*eta*x7 - 4*rho*x5 + 4*eta*x10;
			T dydrho  = y1*(4*(zeta + rho + eta) - 3) - 4*y5*(zeta + 2*rho + eta - 1) + 4*y9*zeta - 4*y8*zeta - 4*eta*y7 + 4*eta*y6 + (4*rho - 1)*y2;
			T dydeta  = y1*(4*(zeta + rho + eta) - 3) - 4*y7*(zeta + rho + 2*eta - 1) - 4*y8*zeta + 4*y10*zeta + 4*rho*y6 - 4*rho*y5 + (4*eta - 1)*y3;
			T dydzeta = y1*(4*(zeta + rho + eta) - 3) + y4*(4*zeta - 1) - 4*y8*(2*zeta + rho + eta - 1) + 4*rho*y9 - 4*eta*y7 - 4*rho*y5 + 4*eta*y10;
			T dzdrho  = z1*(4*(zeta + rho + eta) - 3) - 4*z5*(zeta + 2*rho + eta - 1) + 4*z9*zeta - 4*z8*zeta - 4*eta*z7 + 4*eta*z6 + (4*rho - 1)*z2;
			T dzdeta  = z1*(4*(zeta + rho + eta) - 3) - 4*z7*(zeta + rho + 2*eta - 1) - 4*z8*zeta + 4*z10*zeta + 4*rho*z6 - 4*rho*z5 + (4*eta - 1)*z3;
			T dzdzeta = z1*(4*(zeta + rho + eta) - 3) + z4*(4*zeta - 1) - 4*z8*(2*zeta + rho + eta - 1) + 4*rho*z9 - 4*eta*z7 - 4*rho*z5 + 4*eta*z10;

			// Jacobian determinant
			det_J = dxdrho*(dydeta*dzdzeta - dydzeta*dzdeta) - dydrho*(dxdeta*dzdzeta - dxdzeta*dzdeta) + (dxdeta*dydzeta - dxdzeta*dydeta)*dzdrho;

			// Inverse of the Jacobian
			T invJ11 = (dydeta*dzdzeta - dydzeta*dzdeta)/det_J;
			T invJ12 = (dydzeta*dzdrho - dydrho*dzdzeta)/det_J;
			T invJ13 = (dydrho*dzdeta - dydeta*dzdrho)/det_J;
			T invJ21 = (dxdzeta*dzdeta - dxdeta*dzdzeta)/det_J;
			T invJ22 = (dxdrho*dzdzeta - dxdzeta*dzdrho)/det_J;
			T invJ23 = (dxdeta*dzdrho - dxdrho*dzdeta)/det_J;
			T invJ31 = (dxdeta*dydzeta - dxdzeta*dydeta)/det_J;
			T invJ32 = (dxdzeta*dydrho - dxdrho*dydzeta)/det_J;
			T invJ33 = (dxdrho*dydeta - dxdeta*dydrho)/det_J;

			// Global coordinates derivative
			T dN1drho   =  4*(zeta + rho + eta) - 3;
			T dN1deta   =  4*(zeta + rho + eta) - 3;
			T dN1dzeta  =  4*(zeta + rho + eta) - 3;
			T dN2drho   =  4*rho - 1;
			T dN2deta   =  0;
			T dN2dzeta  =  0;
			T dN3drho   =  0;
			T dN3deta   =  4*eta - 1;
			T dN3dzeta  =  0;
			T dN4drho   =  0;
			T dN4deta   =  0;
			T dN4dzeta  =  4*zeta - 1;
			T dN5drho   = -4*(zeta + 2*rho + eta - 1);
			T dN5deta   = -4*rho;
			T dN5dzeta  = -4*rho;
			T dN6drho   =  4*eta;
			T dN6deta   =  4*rho;
			T dN6dzeta  =  0;
			T dN7drho   = -4*eta;
			T dN7deta   = -4*(zeta + rho + 2*eta - 1);
			T dN7dzeta  = -4*eta;
			T dN8drho   = -4*zeta;
			T dN8deta   = -4*zeta;
			T dN8dzeta  = -4*(2*zeta + rho + eta - 1);
			T dN9drho   =  4*zeta;
			T dN9deta   =  0;
			T dN9dzeta  =  4*rho;
			T dN10drho  =  0;
			T dN10deta  =  4*zeta;
			T dN10dzeta =  4*eta;

			dN.entry[1][1]  = invJ11*dN1drho + invJ12*dN1deta + invJ13*dN1dzeta;
			dN.entry[1][2]  = invJ21*dN1drho + invJ22*dN1deta + invJ23*dN1dzeta;
			dN.entry[1][3]  = invJ31*dN1drho + invJ32*dN1deta + invJ33*dN1dzeta;
			dN.entry[2][1]  = invJ11*dN2drho + invJ12*dN2deta + invJ13*dN2dzeta;
			dN.entry[2][2]  = invJ21*dN2drho + invJ22*dN2deta + invJ23*dN2dzeta;
			dN.entry[2][3]  = invJ31*dN2drho + invJ32*dN2deta + invJ33*dN2dzeta;
			dN.entry[3][1]  = invJ11*dN3drho + invJ12*dN3deta + invJ13*dN3dzeta;
			dN.entry[3][2]  = invJ21*dN3drho + invJ22*dN3deta + invJ23*dN3dzeta;
			dN.entry[3][3]  = invJ31*dN3drho + invJ32*dN3deta + invJ33*dN3dzeta;
			dN.entry[4][1]  = invJ11*dN4drho + invJ12*dN4deta + invJ13*dN4dzeta;
			dN.entry[4][2]  = invJ21*dN4drho + invJ22*dN4deta + invJ23*dN4dzeta;
			dN.entry[4][3]  = invJ31*dN4drho + invJ32*dN4deta + invJ33*dN4dzeta;
			dN.entry[5][1]  = invJ11*dN5drho + invJ12*dN5deta + invJ13*dN5dzeta;
			dN.entry[5][2]  = invJ21*dN5drho + invJ22*dN5deta + invJ23*dN5dzeta;
			dN.entry[5][3]  = invJ31*dN5drho + invJ32*dN5deta + invJ33*dN5dzeta;
			dN.entry[6][1]  = invJ11*dN6drho + invJ12*dN6deta + invJ13*dN6dzeta;
			dN.entry[6][2]  = invJ21*dN6drho + invJ22*dN6deta + invJ23*dN6dzeta;
			dN.entry[6][3]  = invJ31*dN6drho + invJ32*dN6deta + invJ33*dN6dzeta;
			dN.entry[7][1]  = invJ11*dN7drho + invJ12*dN7deta + invJ13*dN7dzeta;
			dN.entry[7][2]  = invJ21*dN7drho + invJ22*dN7deta + invJ23*dN7dzeta;
			dN.entry[7][3]  = invJ31*dN7drho + invJ32*dN7deta + invJ33*dN7dzeta;
			dN.entry[8][1]  = invJ11*dN8drho + invJ12*dN8deta + invJ13*dN8dzeta;
			dN.entry[8][2]  = invJ21*dN8drho + invJ22*dN8deta + invJ23*dN8dzeta;
			dN.entry[8][3]  = invJ31*dN8drho + invJ32*dN8deta + invJ33*dN8dzeta;
			dN.entry[9][1]  = invJ11*dN9drho + invJ12*dN9deta + invJ13*dN9dzeta;
			dN.entry[9][2]  = invJ21*dN9drho + invJ22*dN9deta + invJ23*dN9dzeta;
			dN.entry[9][3]  = invJ31*dN9drho + invJ32*dN9deta + invJ33*dN9dzeta;
			dN.entry[10][1] = invJ11*dN10drho + invJ12*dN10deta + invJ13*dN10dzeta;
			dN.entry[10][2] = invJ21*dN10drho + invJ22*dN10deta + invJ23*dN10dzeta;
			dN.entry[10][3] = invJ31*dN10drho + invJ32*dN10deta + invJ33*dN10dzeta;
		}


		void Hexahedron8(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			register T rho  = point.entry[1];
			register T eta  = point.entry[2];
			register T zeta = point.entry[3];

			// Shape functions
			N.entry[1] = T(0.125)*(1 - rho)*(1 - eta)*(1 - zeta);
			N.entry[2] = T(0.125)*(1 + rho)*(1 - eta)*(1 - zeta);
			N.entry[3] = T(0.125)*(1 + rho)*(1 + eta)*(1 - zeta);
			N.entry[4] = T(0.125)*(1 - rho)*(1 + eta)*(1 - zeta);
			N.entry[5] = T(0.125)*(1 - rho)*(1 - eta)*(1 + zeta);
			N.entry[6] = T(0.125)*(1 + rho)*(1 - eta)*(1 + zeta);
			N.entry[7] = T(0.125)*(1 + rho)*(1 + eta)*(1 + zeta);
			N.entry[8] = T(0.125)*(1 - rho)*(1 + eta)*(1 + zeta);

			// Local coordinates derivative
			int node1  = mesh.connectivity.entry[element_id][1];
			int node2  = mesh.connectivity.entry[element_id][2];
			int node3  = mesh.connectivity.entry[element_id][3];
			int node4  = mesh.connectivity.entry[element_id][4];
			int node5  = mesh.connectivity.entry[element_id][5];
			int node6  = mesh.connectivity.entry[element_id][6];
			int node7  = mesh.connectivity.entry[element_id][7];
			int node8  = mesh.connectivity.entry[element_id][8];
			T x1  = nodes.coordinate.entry[node1][1];
			T y1  = nodes.coordinate.entry[node1][2];
			T z1  = nodes.coordinate.entry[node1][3];
			T x2  = nodes.coordinate.entry[node2][1];
			T y2  = nodes.coordinate.entry[node2][2];
			T z2  = nodes.coordinate.entry[node2][3];
			T x3  = nodes.coordinate.entry[node3][1];
			T y3  = nodes.coordinate.entry[node3][2];
			T z3  = nodes.coordinate.entry[node3][3];
			T x4  = nodes.coordinate.entry[node4][1];
			T y4  = nodes.coordinate.entry[node4][2];
			T z4  = nodes.coordinate.entry[node4][3];
			T x5  = nodes.coordinate.entry[node5][1];
			T y5  = nodes.coordinate.entry[node5][2];
			T z5  = nodes.coordinate.entry[node5][3];
			T x6  = nodes.coordinate.entry[node6][1];
			T y6  = nodes.coordinate.entry[node6][2];
			T z6  = nodes.coordinate.entry[node6][3];
			T x7  = nodes.coordinate.entry[node7][1];
			T y7  = nodes.coordinate.entry[node7][2];
			T z7  = nodes.coordinate.entry[node7][3];
			T x8  = nodes.coordinate.entry[node8][1];
			T y8  = nodes.coordinate.entry[node8][2];
			T z8  = nodes.coordinate.entry[node8][3];

			T dxdrho  = (-(eta + 1)*(zeta + 1)*x8 + (eta + 1)*(zeta + 1)*x7 - (eta - 1)*(zeta + 1)*x6 + (eta - 1)*(zeta + 1)*x5 + (eta + 1)*(zeta - 1)*x4 - (eta + 1)*(zeta - 1)*x3 + (eta - 1)*(zeta - 1)*x2 - (eta - 1)*(zeta - 1)*x1)*T(0.125);
			T dxdeta  = (-(rho - 1)*(zeta + 1)*x8 + (rho + 1)*(zeta + 1)*x7 - (rho + 1)*(zeta + 1)*x6 + (rho - 1)*(zeta + 1)*x5 + (rho - 1)*(zeta - 1)*x4 - (rho + 1)*(zeta - 1)*x3 + (rho + 1)*(zeta - 1)*x2 - (rho - 1)*(zeta - 1)*x1)*T(0.125);
			T dxdzeta = (-(eta + 1)*(rho - 1)*x8 + (eta + 1)*(rho + 1)*x7 - (eta - 1)*(rho + 1)*x6 + (eta - 1)*(rho - 1)*x5 + (eta + 1)*(rho - 1)*x4 - (eta + 1)*(rho + 1)*x3 + (eta - 1)*(rho + 1)*x2 - (eta - 1)*(rho - 1)*x1)*T(0.125);
			T dydrho  = (-(eta + 1)*(zeta + 1)*y8 + (eta + 1)*(zeta + 1)*y7 - (eta - 1)*(zeta + 1)*y6 + (eta - 1)*(zeta + 1)*y5 + (eta + 1)*(zeta - 1)*y4 - (eta + 1)*(zeta - 1)*y3 + (eta - 1)*(zeta - 1)*y2 - (eta - 1)*(zeta - 1)*y1)*T(0.125);
			T dydeta  = (-(rho - 1)*(zeta + 1)*y8 + (rho + 1)*(zeta + 1)*y7 - (rho + 1)*(zeta + 1)*y6 + (rho - 1)*(zeta + 1)*y5 + (rho - 1)*(zeta - 1)*y4 - (rho + 1)*(zeta - 1)*y3 + (rho + 1)*(zeta - 1)*y2 - (rho - 1)*(zeta - 1)*y1)*T(0.125);
			T dydzeta = (-(eta + 1)*(rho - 1)*y8 + (eta + 1)*(rho + 1)*y7 - (eta - 1)*(rho + 1)*y6 + (eta - 1)*(rho - 1)*y5 + (eta + 1)*(rho - 1)*y4 - (eta + 1)*(rho + 1)*y3 + (eta - 1)*(rho + 1)*y2 - (eta - 1)*(rho - 1)*y1)*T(0.125);
			T dzdrho  = (-(eta + 1)*(zeta + 1)*z8 + (eta + 1)*(zeta + 1)*z7 - (eta - 1)*(zeta + 1)*z6 + (eta - 1)*(zeta + 1)*z5 + (eta + 1)*(zeta - 1)*z4 - (eta + 1)*(zeta - 1)*z3 + (eta - 1)*(zeta - 1)*z2 - (eta - 1)*(zeta - 1)*z1)*T(0.125);
			T dzdeta  = (-(rho - 1)*(zeta + 1)*z8 + (rho + 1)*(zeta + 1)*z7 - (rho + 1)*(zeta + 1)*z6 + (rho - 1)*(zeta + 1)*z5 + (rho - 1)*(zeta - 1)*z4 - (rho + 1)*(zeta - 1)*z3 + (rho + 1)*(zeta - 1)*z2 - (rho - 1)*(zeta - 1)*z1)*T(0.125);
			T dzdzeta = (-(eta + 1)*(rho - 1)*z8 + (eta + 1)*(rho + 1)*z7 - (eta - 1)*(rho + 1)*z6 + (eta - 1)*(rho - 1)*z5 + (eta + 1)*(rho - 1)*z4 - (eta + 1)*(rho + 1)*z3 + (eta - 1)*(rho + 1)*z2 - (eta - 1)*(rho - 1)*z1)*T(0.125);

			// Jacobian determinant
			det_J = dxdrho*(dydeta*dzdzeta - dydzeta*dzdeta) - dydrho*(dxdeta*dzdzeta - dxdzeta*dzdeta) + (dxdeta*dydzeta - dxdzeta*dydeta)*dzdrho;

			// Inverse of the Jacobian
			T invJ11 = (dydeta*dzdzeta - dydzeta*dzdeta)/det_J;
			T invJ12 = (dydzeta*dzdrho - dydrho*dzdzeta)/det_J;
			T invJ13 = (dydrho*dzdeta - dydeta*dzdrho)/det_J;
			T invJ21 = (dxdzeta*dzdeta - dxdeta*dzdzeta)/det_J;
			T invJ22 = (dxdrho*dzdzeta - dxdzeta*dzdrho)/det_J;
			T invJ23 = (dxdeta*dzdrho - dxdrho*dzdeta)/det_J;
			T invJ31 = (dxdeta*dydzeta - dxdzeta*dydeta)/det_J;
			T invJ32 = (dxdzeta*dydrho - dxdrho*dydzeta)/det_J;
			T invJ33 = (dxdrho*dydeta - dxdeta*dydrho)/det_J;

			// Global coordinates derivative
			T dN1drho  = -((eta - 1)*(zeta - 1))*T(0.125);
			T dN1deta  = -((rho - 1)*(zeta - 1))*T(0.125);
			T dN1dzeta = -((eta - 1)*(rho - 1))*T(0.125);
			T dN2drho  =  ((eta - 1)*(zeta - 1))*T(0.125);
			T dN2deta  =  ((rho + 1)*(zeta - 1))*T(0.125);
			T dN2dzeta =  ((eta - 1)*(rho + 1))*T(0.125);
			T dN3drho  = -((eta + 1)*(zeta - 1))*T(0.125);
			T dN3deta  = -((rho + 1)*(zeta - 1))*T(0.125);
			T dN3dzeta = -((eta + 1)*(rho + 1))*T(0.125);
			T dN4drho  =  ((eta + 1)*(zeta - 1))*T(0.125);
			T dN4deta  =  ((rho - 1)*(zeta - 1))*T(0.125);
			T dN4dzeta =  ((eta + 1)*(rho - 1))*T(0.125);
			T dN5drho  =  ((eta - 1)*(zeta + 1))*T(0.125);
			T dN5deta  =  ((rho - 1)*(zeta + 1))*T(0.125);
			T dN5dzeta =  ((eta - 1)*(rho - 1))*T(0.125);
			T dN6drho  = -((eta - 1)*(zeta + 1))*T(0.125);
			T dN6deta  = -((rho + 1)*(zeta + 1))*T(0.125);
			T dN6dzeta = -((eta - 1)*(rho + 1))*T(0.125);
			T dN7drho  =  ((eta + 1)*(zeta + 1))*T(0.125);
			T dN7deta  =  ((rho + 1)*(zeta + 1))*T(0.125);
			T dN7dzeta =  ((eta + 1)*(rho + 1))*T(0.125);
			T dN8drho  = -((eta + 1)*(zeta + 1))*T(0.125);
			T dN8deta  = -((rho - 1)*(zeta + 1))*T(0.125);
			T dN8dzeta = -((eta + 1)*(rho - 1))*T(0.125);

			dN.entry[1][1]  = invJ11*dN1drho + invJ12*dN1deta + invJ13*dN1dzeta;
			dN.entry[1][2]  = invJ21*dN1drho + invJ22*dN1deta + invJ23*dN1dzeta;
			dN.entry[1][3]  = invJ31*dN1drho + invJ32*dN1deta + invJ33*dN1dzeta;
			dN.entry[2][1]  = invJ11*dN2drho + invJ12*dN2deta + invJ13*dN2dzeta;
			dN.entry[2][2]  = invJ21*dN2drho + invJ22*dN2deta + invJ23*dN2dzeta;
			dN.entry[2][3]  = invJ31*dN2drho + invJ32*dN2deta + invJ33*dN2dzeta;
			dN.entry[3][1]  = invJ11*dN3drho + invJ12*dN3deta + invJ13*dN3dzeta;
			dN.entry[3][2]  = invJ21*dN3drho + invJ22*dN3deta + invJ23*dN3dzeta;
			dN.entry[3][3]  = invJ31*dN3drho + invJ32*dN3deta + invJ33*dN3dzeta;
			dN.entry[4][1]  = invJ11*dN4drho + invJ12*dN4deta + invJ13*dN4dzeta;
			dN.entry[4][2]  = invJ21*dN4drho + invJ22*dN4deta + invJ23*dN4dzeta;
			dN.entry[4][3]  = invJ31*dN4drho + invJ32*dN4deta + invJ33*dN4dzeta;
			dN.entry[5][1]  = invJ11*dN5drho + invJ12*dN5deta + invJ13*dN5dzeta;
			dN.entry[5][2]  = invJ21*dN5drho + invJ22*dN5deta + invJ23*dN5dzeta;
			dN.entry[5][3]  = invJ31*dN5drho + invJ32*dN5deta + invJ33*dN5dzeta;
			dN.entry[6][1]  = invJ11*dN6drho + invJ12*dN6deta + invJ13*dN6dzeta;
			dN.entry[6][2]  = invJ21*dN6drho + invJ22*dN6deta + invJ23*dN6dzeta;
			dN.entry[6][3]  = invJ31*dN6drho + invJ32*dN6deta + invJ33*dN6dzeta;
			dN.entry[7][1]  = invJ11*dN7drho + invJ12*dN7deta + invJ13*dN7dzeta;
			dN.entry[7][2]  = invJ21*dN7drho + invJ22*dN7deta + invJ23*dN7dzeta;
			dN.entry[7][3]  = invJ31*dN7drho + invJ32*dN7deta + invJ33*dN7dzeta;
			dN.entry[8][1]  = invJ11*dN8drho + invJ12*dN8deta + invJ13*dN8dzeta;
			dN.entry[8][2]  = invJ21*dN8drho + invJ22*dN8deta + invJ23*dN8dzeta;
			dN.entry[8][3]  = invJ31*dN8drho + invJ32*dN8deta + invJ33*dN8dzeta;
		}


		void Hexahedron20(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			register T rho  = point.entry[1];
			register T eta  = point.entry[2];
			register T zeta = point.entry[3];

			// Shape functions
			N.entry[1] = T(0.125)*(1 - rho)*(1 - eta)*(1 - zeta)*(-rho - eta - zeta - 2);
			N.entry[2] = T(0.125)*(1 + rho)*(1 - eta)*(1 - zeta)*(rho - eta - zeta - 2);
			N.entry[3] = T(0.125)*(1 + rho)*(1 + eta)*(1 - zeta)*(rho + eta - zeta - 2);
			N.entry[4] = T(0.125)*(1 - rho)*(1 + eta)*(1 - zeta)*(-rho + eta - zeta - 2);
			N.entry[5] = T(0.125)*(1 - rho)*(1 - eta)*(1 + zeta)*(-rho - eta + zeta - 2);
			N.entry[6] = T(0.125)*(1 + rho)*(1 - eta)*(1 + zeta)*(rho - eta + zeta - 2);
			N.entry[7] = T(0.125)*(1 + rho)*(1 + eta)*(1 + zeta)*(rho + eta + zeta - 2);
			N.entry[8] = T(0.125)*(1 - rho)*(1 + eta)*(1 + zeta)*(-rho + eta + zeta - 2);
			N.entry[9] = T(0.25)*(1 - rho*rho)*(1 - eta)*(1 - zeta);
			N.entry[10] = T(0.25)*(1 + rho)*(1 - eta*eta)*(1 - zeta);
			N.entry[11] = T(0.25)*(1 - rho*rho)*(1 + eta)*(1 - zeta);
			N.entry[12] = T(0.25)*(1 - rho)*(1 - eta*eta)*(1 - zeta);
			N.entry[13] = T(0.25)*(1 - rho)*(1 - eta)*(1 - zeta*zeta);
			N.entry[14] = T(0.25)*(1 + rho)*(1 - eta)*(1 - zeta*zeta);
			N.entry[15] = T(0.25)*(1 + rho)*(1 + eta)*(1 - zeta*zeta);
			N.entry[16] = T(0.25)*(1 - rho)*(1 + eta)*(1 - zeta*zeta);
			N.entry[17] = T(0.25)*(1 - rho*rho)*(1 - eta)*(1 + zeta);
			N.entry[18] = T(0.25)*(1 + rho)*(1 - eta*eta)*(1 + zeta);
			N.entry[19] = T(0.25)*(1 - rho*rho)*(1 + eta)*(1 + zeta);
			N.entry[20] = T(0.25)*(1 - rho)*(1 - eta*eta)*(1 + zeta);

			// Local coordinates derivative
			int node1  = mesh.connectivity.entry[element_id][1];
			int node2  = mesh.connectivity.entry[element_id][2];
			int node3  = mesh.connectivity.entry[element_id][3];
			int node4  = mesh.connectivity.entry[element_id][4];
			int node5  = mesh.connectivity.entry[element_id][5];
			int node6  = mesh.connectivity.entry[element_id][6];
			int node7  = mesh.connectivity.entry[element_id][7];
			int node8  = mesh.connectivity.entry[element_id][8];
			int node9  = mesh.connectivity.entry[element_id][9];
			int node10 = mesh.connectivity.entry[element_id][10];
			int node11 = mesh.connectivity.entry[element_id][11];
			int node12 = mesh.connectivity.entry[element_id][12];
			int node13 = mesh.connectivity.entry[element_id][13];
			int node14 = mesh.connectivity.entry[element_id][14];
			int node15 = mesh.connectivity.entry[element_id][15];
			int node16 = mesh.connectivity.entry[element_id][16];
			int node17 = mesh.connectivity.entry[element_id][17];
			int node18 = mesh.connectivity.entry[element_id][18];
			int node19 = mesh.connectivity.entry[element_id][19];
			int node20 = mesh.connectivity.entry[element_id][20];
			T x1  = nodes.coordinate.entry[node1][1];
			T y1  = nodes.coordinate.entry[node1][2];
			T z1  = nodes.coordinate.entry[node1][3];
			T x2  = nodes.coordinate.entry[node2][1];
			T y2  = nodes.coordinate.entry[node2][2];
			T z2  = nodes.coordinate.entry[node2][3];
			T x3  = nodes.coordinate.entry[node3][1];
			T y3  = nodes.coordinate.entry[node3][2];
			T z3  = nodes.coordinate.entry[node3][3];
			T x4  = nodes.coordinate.entry[node4][1];
			T y4  = nodes.coordinate.entry[node4][2];
			T z4  = nodes.coordinate.entry[node4][3];
			T x5  = nodes.coordinate.entry[node5][1];
			T y5  = nodes.coordinate.entry[node5][2];
			T z5  = nodes.coordinate.entry[node5][3];
			T x6  = nodes.coordinate.entry[node6][1];
			T y6  = nodes.coordinate.entry[node6][2];
			T z6  = nodes.coordinate.entry[node6][3];
			T x7  = nodes.coordinate.entry[node7][1];
			T y7  = nodes.coordinate.entry[node7][2];
			T z7  = nodes.coordinate.entry[node7][3];
			T x8  = nodes.coordinate.entry[node8][1];
			T y8  = nodes.coordinate.entry[node8][2];
			T z8  = nodes.coordinate.entry[node8][3];
			T x9  = nodes.coordinate.entry[node9][1];
			T y9  = nodes.coordinate.entry[node9][2];
			T z9  = nodes.coordinate.entry[node9][3];
			T x10 = nodes.coordinate.entry[node10][1];
			T y10 = nodes.coordinate.entry[node10][2];
			T z10 = nodes.coordinate.entry[node10][3];
			T x11 = nodes.coordinate.entry[node11][1];
			T y11 = nodes.coordinate.entry[node11][2];
			T z11 = nodes.coordinate.entry[node11][3];
			T x12 = nodes.coordinate.entry[node12][1];
			T y12 = nodes.coordinate.entry[node12][2];
			T z12 = nodes.coordinate.entry[node12][3];
			T x13 = nodes.coordinate.entry[node13][1];
			T y13 = nodes.coordinate.entry[node13][2];
			T z13 = nodes.coordinate.entry[node13][3];
			T x14 = nodes.coordinate.entry[node14][1];
			T y14 = nodes.coordinate.entry[node14][2];
			T z14 = nodes.coordinate.entry[node14][3];
			T x15 = nodes.coordinate.entry[node15][1];
			T y15 = nodes.coordinate.entry[node15][2];
			T z15 = nodes.coordinate.entry[node15][3];
			T x16 = nodes.coordinate.entry[node16][1];
			T y16 = nodes.coordinate.entry[node16][2];
			T z16 = nodes.coordinate.entry[node16][3];
			T x17 = nodes.coordinate.entry[node17][1];
			T y17 = nodes.coordinate.entry[node17][2];
			T z17 = nodes.coordinate.entry[node17][3];
			T x18 = nodes.coordinate.entry[node18][1];
			T y18 = nodes.coordinate.entry[node18][2];
			T z18 = nodes.coordinate.entry[node18][3];
			T x19 = nodes.coordinate.entry[node19][1];
			T y19 = nodes.coordinate.entry[node19][2];
			T z19 = nodes.coordinate.entry[node19][3];
			T x20 = nodes.coordinate.entry[node20][1];
			T y20 = nodes.coordinate.entry[node20][2];
			T z20 = nodes.coordinate.entry[node20][3];

			T dxdrho  = ((eta - 1)*x1*(zeta - 1)*(zeta + 2*rho + eta + 1) + (eta + 1)*x7*(zeta + 1)*(zeta + 2*rho + eta - 1) - (eta + 1)*x4*(zeta - 1)*(zeta + 2*rho - eta + 1) - (eta - 1)*x6*(zeta + 1)*(zeta + 2*rho - eta - 1) - (eta - 1)*x2*(zeta - 1)*(zeta - 2*rho + eta + 1) - (eta + 1)*x8*(zeta + 1)*(zeta - 2*rho + eta - 1) + (eta + 1)*x3*(zeta - 1)*(zeta - 2*rho - eta + 1) + (eta - 1)*x5*(zeta + 1)*(zeta - 2*rho - eta - 1) + 2*(eta + 1)*x16*(zeta - 1)*(zeta + 1) - 2*(eta + 1)*x15*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*x14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*x13*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*(eta + 1)*x20*(zeta + 1) - 4*(eta + 1)*rho*x19*(zeta + 1) - 2*(eta - 1)*(eta + 1)*x18*(zeta + 1) + 4*(eta - 1)*rho*x17*(zeta + 1) - 4*(eta - 1)*rho*x9*(zeta - 1) - 2*(eta - 1)*(eta + 1)*x12*(zeta - 1) + 4*(eta + 1)*rho*x11*(zeta - 1) + 2*(eta - 1)*(eta + 1)*x10*(zeta - 1))*T(0.125);
			T dxdeta  = ((rho - 1)*x1*(zeta - 1)*(zeta + rho + 2*eta + 1) + (rho + 1)*x7*(zeta + 1)*(zeta + rho + 2*eta - 1) - (rho - 1)*x4*(zeta - 1)*(zeta + rho - 2*eta + 1) - (rho + 1)*x6*(zeta + 1)*(zeta + rho - 2*eta - 1) - (rho + 1)*x2*(zeta - 1)*(zeta - rho + 2*eta + 1) - (rho - 1)*x8*(zeta + 1)*(zeta - rho + 2*eta - 1) + (rho + 1)*x3*(zeta - 1)*(zeta - rho - 2*eta + 1) + (rho - 1)*x5*(zeta + 1)*(zeta - rho - 2*eta - 1) + 2*(rho - 1)*x16*(zeta - 1)*(zeta + 1) - 2*(rho + 1)*x15*(zeta - 1)*(zeta + 1) + 2*(rho + 1)*x14*(zeta - 1)*(zeta + 1) - 2*(rho - 1)*x13*(zeta - 1)*(zeta + 1) + 4*eta*(rho - 1)*x20*(zeta + 1) - 2*(rho - 1)*(rho + 1)*x19*(zeta + 1) - 4*eta*(rho + 1)*x18*(zeta + 1) + 2*(rho - 1)*(rho + 1)*x17*(zeta + 1) - 2*(rho - 1)*(rho + 1)*x9*(zeta - 1) - 4*eta*(rho - 1)*x12*(zeta - 1) + 2*(rho - 1)*(rho + 1)*x11*(zeta - 1) + 4*eta*(rho + 1)*x10*(zeta - 1))*T(0.125);
			T dxdzeta = ((eta - 1)*(rho - 1)*x1*(2*zeta + rho + eta + 1) + (eta + 1)*(rho + 1)*x7*(2*zeta + rho + eta - 1) - (eta + 1)*(rho - 1)*x4*(2*zeta + rho - eta + 1) - (eta - 1)*(rho + 1)*x6*(2*zeta + rho - eta - 1) - (eta - 1)*(rho + 1)*x2*(2*zeta - rho + eta + 1) - (eta + 1)*(rho - 1)*x8*(2*zeta - rho + eta - 1) + (eta + 1)*(rho + 1)*x3*(2*zeta - rho - eta + 1) + (eta - 1)*(rho - 1)*x5*(2*zeta - rho - eta - 1) + 4*(eta + 1)*(rho - 1)*x16*zeta - 4*(eta + 1)*(rho + 1)*x15*zeta + 4*(eta - 1)*(rho + 1)*x14*zeta - 4*(eta - 1)*(rho - 1)*x13*zeta - 2*(eta - 1)*(rho - 1)*(rho + 1)*x9 + 2*(eta - 1)*(eta + 1)*(rho - 1)*x20 - 2*(eta + 1)*(rho - 1)*(rho + 1)*x19 - 2*(eta - 1)*(eta + 1)*(rho + 1)*x18 + 2*(eta - 1)*(rho - 1)*(rho + 1)*x17 - 2*(eta - 1)*(eta + 1)*(rho - 1)*x12 + 2*(eta + 1)*(rho - 1)*(rho + 1)*x11 + 2*(eta - 1)*(eta + 1)*(rho + 1)*x10)*T(0.125);
			T dydrho  = ((eta - 1)*y1*(zeta - 1)*(zeta + 2*rho + eta + 1) + (eta + 1)*y7*(zeta + 1)*(zeta + 2*rho + eta - 1) - (eta + 1)*y4*(zeta - 1)*(zeta + 2*rho - eta + 1) - (eta - 1)*y6*(zeta + 1)*(zeta + 2*rho - eta - 1) - (eta - 1)*y2*(zeta - 1)*(zeta - 2*rho + eta + 1) - (eta + 1)*y8*(zeta + 1)*(zeta - 2*rho + eta - 1) + (eta + 1)*y3*(zeta - 1)*(zeta - 2*rho - eta + 1) + (eta - 1)*y5*(zeta + 1)*(zeta - 2*rho - eta - 1) + 2*(eta + 1)*y16*(zeta - 1)*(zeta + 1) - 2*(eta + 1)*y15*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*y14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*y13*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*(eta + 1)*y20*(zeta + 1) - 4*(eta + 1)*rho*y19*(zeta + 1) - 2*(eta - 1)*(eta + 1)*y18*(zeta + 1) + 4*(eta - 1)*rho*y17*(zeta + 1) - 4*(eta - 1)*rho*y9*(zeta - 1) - 2*(eta - 1)*(eta + 1)*y12*(zeta - 1) + 4*(eta + 1)*rho*y11*(zeta - 1) + 2*(eta - 1)*(eta + 1)*y10*(zeta - 1))*T(0.125);
			T dydeta  = ((rho - 1)*y1*(zeta - 1)*(zeta + rho + 2*eta + 1) + (rho + 1)*y7*(zeta + 1)*(zeta + rho + 2*eta - 1) - (rho - 1)*y4*(zeta - 1)*(zeta + rho - 2*eta + 1) - (rho + 1)*y6*(zeta + 1)*(zeta + rho - 2*eta - 1) - (rho + 1)*y2*(zeta - 1)*(zeta - rho + 2*eta + 1) - (rho - 1)*y8*(zeta + 1)*(zeta - rho + 2*eta - 1) + (rho + 1)*y3*(zeta - 1)*(zeta - rho - 2*eta + 1) + (rho - 1)*y5*(zeta + 1)*(zeta - rho - 2*eta - 1) + 2*(rho - 1)*y16*(zeta - 1)*(zeta + 1) - 2*(rho + 1)*y15*(zeta - 1)*(zeta + 1) + 2*(rho + 1)*y14*(zeta - 1)*(zeta + 1) - 2*(rho - 1)*y13*(zeta - 1)*(zeta + 1) + 4*eta*(rho - 1)*y20*(zeta + 1) - 2*(rho - 1)*(rho + 1)*y19*(zeta + 1) - 4*eta*(rho + 1)*y18*(zeta + 1) + 2*(rho - 1)*(rho + 1)*y17*(zeta + 1) - 2*(rho - 1)*(rho + 1)*y9*(zeta - 1) - 4*eta*(rho - 1)*y12*(zeta - 1) + 2*(rho - 1)*(rho + 1)*y11*(zeta - 1) + 4*eta*(rho + 1)*y10*(zeta - 1))*T(0.125);
			T dydzeta = ((eta - 1)*(rho - 1)*y1*(2*zeta + rho + eta + 1) + (eta + 1)*(rho + 1)*y7*(2*zeta + rho + eta - 1) - (eta + 1)*(rho - 1)*y4*(2*zeta + rho - eta + 1) - (eta - 1)*(rho + 1)*y6*(2*zeta + rho - eta - 1) - (eta - 1)*(rho + 1)*y2*(2*zeta - rho + eta + 1) - (eta + 1)*(rho - 1)*y8*(2*zeta - rho + eta - 1) + (eta + 1)*(rho + 1)*y3*(2*zeta - rho - eta + 1) + (eta - 1)*(rho - 1)*y5*(2*zeta - rho - eta - 1) + 4*(eta + 1)*(rho - 1)*y16*zeta - 4*(eta + 1)*(rho + 1)*y15*zeta + 4*(eta - 1)*(rho + 1)*y14*zeta - 4*(eta - 1)*(rho - 1)*y13*zeta - 2*(eta - 1)*(rho - 1)*(rho + 1)*y9 + 2*(eta - 1)*(eta + 1)*(rho - 1)*y20 - 2*(eta + 1)*(rho - 1)*(rho + 1)*y19 - 2*(eta - 1)*(eta + 1)*(rho + 1)*y18 + 2*(eta - 1)*(rho - 1)*(rho + 1)*y17 - 2*(eta - 1)*(eta + 1)*(rho - 1)*y12 + 2*(eta + 1)*(rho - 1)*(rho + 1)*y11 + 2*(eta - 1)*(eta + 1)*(rho + 1)*y10)*T(0.125);
			T dzdrho  = ((eta - 1)*z1*(zeta - 1)*(zeta + 2*rho + eta + 1) + (eta + 1)*z7*(zeta + 1)*(zeta + 2*rho + eta - 1) - (eta + 1)*z4*(zeta - 1)*(zeta + 2*rho - eta + 1) - (eta - 1)*z6*(zeta + 1)*(zeta + 2*rho - eta - 1) - (eta - 1)*z2*(zeta - 1)*(zeta - 2*rho + eta + 1) - (eta + 1)*z8*(zeta + 1)*(zeta - 2*rho + eta - 1) + (eta + 1)*z3*(zeta - 1)*(zeta - 2*rho - eta + 1) + (eta - 1)*z5*(zeta + 1)*(zeta - 2*rho - eta - 1) + 2*(eta + 1)*z16*(zeta - 1)*(zeta + 1) - 2*(eta + 1)*z15*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*z14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*z13*(zeta - 1)*(zeta + 1) + 2*(eta - 1)*(eta + 1)*z20*(zeta + 1) - 4*(eta + 1)*rho*z19*(zeta + 1) - 2*(eta - 1)*(eta + 1)*z18*(zeta + 1) + 4*(eta - 1)*rho*z17*(zeta + 1) - 4*(eta - 1)*rho*z9*(zeta - 1) - 2*(eta - 1)*(eta + 1)*z12*(zeta - 1) + 4*(eta + 1)*rho*z11*(zeta - 1) + 2*(eta - 1)*(eta + 1)*z10*(zeta - 1))*T(0.125);
			T dzdeta  = ((rho - 1)*z1*(zeta - 1)*(zeta + rho + 2*eta + 1) + (rho + 1)*z7*(zeta + 1)*(zeta + rho + 2*eta - 1) - (rho - 1)*z4*(zeta - 1)*(zeta + rho - 2*eta + 1) - (rho + 1)*z6*(zeta + 1)*(zeta + rho - 2*eta - 1) - (rho + 1)*z2*(zeta - 1)*(zeta - rho + 2*eta + 1) - (rho - 1)*z8*(zeta + 1)*(zeta - rho + 2*eta - 1) + (rho + 1)*z3*(zeta - 1)*(zeta - rho - 2*eta + 1) + (rho - 1)*z5*(zeta + 1)*(zeta - rho - 2*eta - 1) + 2*(rho - 1)*z16*(zeta - 1)*(zeta + 1) - 2*(rho + 1)*z15*(zeta - 1)*(zeta + 1) + 2*(rho + 1)*z14*(zeta - 1)*(zeta + 1) - 2*(rho - 1)*z13*(zeta - 1)*(zeta + 1) + 4*eta*(rho - 1)*z20*(zeta + 1) - 2*(rho - 1)*(rho + 1)*z19*(zeta + 1) - 4*eta*(rho + 1)*z18*(zeta + 1) + 2*(rho - 1)*(rho + 1)*z17*(zeta + 1) - 2*(rho - 1)*(rho + 1)*z9*(zeta - 1) - 4*eta*(rho - 1)*z12*(zeta - 1) + 2*(rho - 1)*(rho + 1)*z11*(zeta - 1) + 4*eta*(rho + 1)*z10*(zeta - 1))*T(0.125);
			T dzdzeta = ((eta - 1)*(rho - 1)*z1*(2*zeta + rho + eta + 1) + (eta + 1)*(rho + 1)*z7*(2*zeta + rho + eta - 1) - (eta + 1)*(rho - 1)*z4*(2*zeta + rho - eta + 1) - (eta - 1)*(rho + 1)*z6*(2*zeta + rho - eta - 1) - (eta - 1)*(rho + 1)*z2*(2*zeta - rho + eta + 1) - (eta + 1)*(rho - 1)*z8*(2*zeta - rho + eta - 1) + (eta + 1)*(rho + 1)*z3*(2*zeta - rho - eta + 1) + (eta - 1)*(rho - 1)*z5*(2*zeta - rho - eta - 1) + 4*(eta + 1)*(rho - 1)*z16*zeta - 4*(eta + 1)*(rho + 1)*z15*zeta + 4*(eta - 1)*(rho + 1)*z14*zeta - 4*(eta - 1)*(rho - 1)*z13*zeta - 2*(eta - 1)*(rho - 1)*(rho + 1)*z9 + 2*(eta - 1)*(eta + 1)*(rho - 1)*z20 - 2*(eta + 1)*(rho - 1)*(rho + 1)*z19 - 2*(eta - 1)*(eta + 1)*(rho + 1)*z18 + 2*(eta - 1)*(rho - 1)*(rho + 1)*z17 - 2*(eta - 1)*(eta + 1)*(rho - 1)*z12 + 2*(eta + 1)*(rho - 1)*(rho + 1)*z11 + 2*(eta - 1)*(eta + 1)*(rho + 1)*z10)*T(0.125);

			// Jacobian determinant
			det_J = dxdrho*(dydeta*dzdzeta - dydzeta*dzdeta) - dydrho*(dxdeta*dzdzeta - dxdzeta*dzdeta) + (dxdeta*dydzeta - dxdzeta*dydeta)*dzdrho;

			// Inverse of the Jacobian
			T invJ11 = (dydeta*dzdzeta - dydzeta*dzdeta)/det_J;
			T invJ12 = (dydzeta*dzdrho - dydrho*dzdzeta)/det_J;
			T invJ13 = (dydrho*dzdeta - dydeta*dzdrho)/det_J;
			T invJ21 = (dxdzeta*dzdeta - dxdeta*dzdzeta)/det_J;
			T invJ22 = (dxdrho*dzdzeta - dxdzeta*dzdrho)/det_J;
			T invJ23 = (dxdeta*dzdrho - dxdrho*dzdeta)/det_J;
			T invJ31 = (dxdeta*dydzeta - dxdzeta*dydeta)/det_J;
			T invJ32 = (dxdzeta*dydrho - dxdrho*dydzeta)/det_J;
			T invJ33 = (dxdrho*dydeta - dxdeta*dydrho)/det_J;

			// Global coordinates derivative
			T dN1drho   =  ((eta - 1)*(zeta - 1)*(zeta + 2*rho + eta + 1))*T(0.125);
			T dN1deta   =  ((rho - 1)*(zeta - 1)*(zeta + rho + 2*eta + 1))*T(0.125);
			T dN1dzeta  =  ((eta - 1)*(rho - 1)*(2*zeta + rho + eta + 1))*T(0.125);
			T dN2drho   = -((eta - 1)*(zeta - 1)*(zeta - 2*rho + eta + 1))*T(0.125);
			T dN2deta   = -((rho + 1)*(zeta - 1)*(zeta - rho + 2*eta + 1))*T(0.125);
			T dN2dzeta  = -((eta - 1)*(rho + 1)*(2*zeta - rho + eta + 1))*T(0.125);
			T dN3drho   =  ((eta + 1)*(zeta - 1)*(zeta - 2*rho - eta + 1))*T(0.125);
			T dN3deta   =  ((rho + 1)*(zeta - 1)*(zeta - rho - 2*eta + 1))*T(0.125);
			T dN3dzeta  =  ((eta + 1)*(rho + 1)*(2*zeta - rho - eta + 1))*T(0.125);
			T dN4drho   = -((eta + 1)*(zeta - 1)*(zeta + 2*rho - eta + 1))*T(0.125);
			T dN4deta   = -((rho - 1)*(zeta - 1)*(zeta + rho - 2*eta + 1))*T(0.125);
			T dN4dzeta  = -((eta + 1)*(rho - 1)*(2*zeta + rho - eta + 1))*T(0.125);
			T dN5drho   =  ((eta - 1)*(zeta + 1)*(zeta - 2*rho - eta - 1))*T(0.125);
			T dN5deta   =  ((rho - 1)*(zeta + 1)*(zeta - rho - 2*eta - 1))*T(0.125);
			T dN5dzeta  =  ((eta - 1)*(rho - 1)*(2*zeta - rho - eta - 1))*T(0.125);
			T dN6drho   = -((eta - 1)*(zeta + 1)*(zeta + 2*rho - eta - 1))*T(0.125);
			T dN6deta   = -((rho + 1)*(zeta + 1)*(zeta + rho - 2*eta - 1))*T(0.125);
			T dN6dzeta  = -((eta - 1)*(rho + 1)*(2*zeta + rho - eta - 1))*T(0.125);
			T dN7drho   =  ((eta + 1)*(zeta + 1)*(zeta + 2*rho + eta - 1))*T(0.125);
			T dN7deta   =  ((rho + 1)*(zeta + 1)*(zeta + rho + 2*eta - 1))*T(0.125);
			T dN7dzeta  =  ((eta + 1)*(rho + 1)*(2*zeta + rho + eta - 1))*T(0.125);
			T dN8drho   = -((eta + 1)*(zeta + 1)*(zeta - 2*rho + eta - 1))*T(0.125);
			T dN8deta   = -((rho - 1)*(zeta + 1)*(zeta - rho + 2*eta - 1))*T(0.125);
			T dN8dzeta  = -((eta + 1)*(rho - 1)*(2*zeta - rho + eta - 1))*T(0.125);
			T dN9drho   = -((eta - 1)*rho*(zeta - 1))*T(0.5);
			T dN9deta   = -((rho - 1)*(rho + 1)*(zeta - 1))*T(0.25);
			T dN9dzeta  = -((eta - 1)*(rho - 1)*(rho + 1))*T(0.25);
			T dN10drho  =  ((eta - 1)*(eta + 1)*(zeta - 1))*T(0.25);
			T dN10deta  =  (eta*(rho + 1)*(zeta - 1))*T(0.5);
			T dN10dzeta =  ((eta - 1)*(eta + 1)*(rho + 1))*T(0.25);
			T dN11drho  =  ((eta + 1)*rho*(zeta - 1))*T(0.5);
			T dN11deta  =  ((rho - 1)*(rho + 1)*(zeta - 1))*T(0.25);
			T dN11dzeta =  ((eta + 1)*(rho - 1)*(rho + 1))*T(0.25);
			T dN12drho  = -((eta - 1)*(eta + 1)*(zeta - 1))*T(0.25);
			T dN12deta  = -(eta*(rho - 1)*(zeta - 1))*T(0.5);
			T dN12dzeta = -((eta - 1)*(eta + 1)*(rho - 1))*T(0.25);
			T dN13drho  = -((eta - 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN13deta  = -((rho - 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN13dzeta = -((eta - 1)*(rho - 1)*zeta)*T(0.5);
			T dN14drho  =  ((eta - 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN14deta  =  ((rho + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN14dzeta =  ((eta - 1)*(rho + 1)*zeta)*T(0.5);
			T dN15drho  = -((eta + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN15deta  = -((rho + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN15dzeta = -((eta + 1)*(rho + 1)*zeta)*T(0.5);
			T dN16drho  =  ((eta + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN16deta  =  ((rho - 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN16dzeta =  ((eta + 1)*(rho - 1)*zeta)*T(0.5);
			T dN17drho  =  ((eta - 1)*rho*(zeta + 1))*T(0.5);
			T dN17deta  =  ((rho - 1)*(rho + 1)*(zeta + 1))*T(0.25);
			T dN17dzeta =  ((eta - 1)*(rho - 1)*(rho + 1))*T(0.25);
			T dN18drho  = -((eta - 1)*(eta + 1)*(zeta + 1))*T(0.25);
			T dN18deta  = -(eta*(rho + 1)*(zeta + 1))*T(0.5);
			T dN18dzeta = -((eta - 1)*(eta + 1)*(rho + 1))*T(0.25);
			T dN19drho  = -((eta + 1)*rho*(zeta + 1))*T(0.5);
			T dN19deta  = -((rho - 1)*(rho + 1)*(zeta + 1))*T(0.25);
			T dN19dzeta = -((eta + 1)*(rho - 1)*(rho + 1))*T(0.25);
			T dN20drho  =  ((eta - 1)*(eta + 1)*(zeta + 1))*T(0.25);
			T dN20deta  =  (eta*(rho - 1)*(zeta + 1))*T(0.5);
			T dN20dzeta =  ((eta - 1)*(eta + 1)*(rho - 1))*T(0.25);

			dN.entry[1][1]  = invJ11*dN1drho + invJ12*dN1deta + invJ13*dN1dzeta;
			dN.entry[1][2]  = invJ21*dN1drho + invJ22*dN1deta + invJ23*dN1dzeta;
			dN.entry[1][3]  = invJ31*dN1drho + invJ32*dN1deta + invJ33*dN1dzeta;
			dN.entry[2][1]  = invJ11*dN2drho + invJ12*dN2deta + invJ13*dN2dzeta;
			dN.entry[2][2]  = invJ21*dN2drho + invJ22*dN2deta + invJ23*dN2dzeta;
			dN.entry[2][3]  = invJ31*dN2drho + invJ32*dN2deta + invJ33*dN2dzeta;
			dN.entry[3][1]  = invJ11*dN3drho + invJ12*dN3deta + invJ13*dN3dzeta;
			dN.entry[3][2]  = invJ21*dN3drho + invJ22*dN3deta + invJ23*dN3dzeta;
			dN.entry[3][3]  = invJ31*dN3drho + invJ32*dN3deta + invJ33*dN3dzeta;
			dN.entry[4][1]  = invJ11*dN4drho + invJ12*dN4deta + invJ13*dN4dzeta;
			dN.entry[4][2]  = invJ21*dN4drho + invJ22*dN4deta + invJ23*dN4dzeta;
			dN.entry[4][3]  = invJ31*dN4drho + invJ32*dN4deta + invJ33*dN4dzeta;
			dN.entry[5][1]  = invJ11*dN5drho + invJ12*dN5deta + invJ13*dN5dzeta;
			dN.entry[5][2]  = invJ21*dN5drho + invJ22*dN5deta + invJ23*dN5dzeta;
			dN.entry[5][3]  = invJ31*dN5drho + invJ32*dN5deta + invJ33*dN5dzeta;
			dN.entry[6][1]  = invJ11*dN6drho + invJ12*dN6deta + invJ13*dN6dzeta;
			dN.entry[6][2]  = invJ21*dN6drho + invJ22*dN6deta + invJ23*dN6dzeta;
			dN.entry[6][3]  = invJ31*dN6drho + invJ32*dN6deta + invJ33*dN6dzeta;
			dN.entry[7][1]  = invJ11*dN7drho + invJ12*dN7deta + invJ13*dN7dzeta;
			dN.entry[7][2]  = invJ21*dN7drho + invJ22*dN7deta + invJ23*dN7dzeta;
			dN.entry[7][3]  = invJ31*dN7drho + invJ32*dN7deta + invJ33*dN7dzeta;
			dN.entry[8][1]  = invJ11*dN8drho + invJ12*dN8deta + invJ13*dN8dzeta;
			dN.entry[8][2]  = invJ21*dN8drho + invJ22*dN8deta + invJ23*dN8dzeta;
			dN.entry[8][3]  = invJ31*dN8drho + invJ32*dN8deta + invJ33*dN8dzeta;
			dN.entry[9][1]  = invJ11*dN9drho + invJ12*dN9deta + invJ13*dN9dzeta;
			dN.entry[9][2]  = invJ21*dN9drho + invJ22*dN9deta + invJ23*dN9dzeta;
			dN.entry[9][3]  = invJ31*dN9drho + invJ32*dN9deta + invJ33*dN9dzeta;
			dN.entry[10][1] = invJ11*dN10drho + invJ12*dN10deta + invJ13*dN10dzeta;
			dN.entry[10][2] = invJ21*dN10drho + invJ22*dN10deta + invJ23*dN10dzeta;
			dN.entry[10][3] = invJ31*dN10drho + invJ32*dN10deta + invJ33*dN10dzeta;
			dN.entry[11][1] = invJ11*dN11drho + invJ12*dN11deta + invJ13*dN11dzeta;
			dN.entry[11][2] = invJ21*dN11drho + invJ22*dN11deta + invJ23*dN11dzeta;
			dN.entry[11][3] = invJ31*dN11drho + invJ32*dN11deta + invJ33*dN11dzeta;
			dN.entry[12][1] = invJ11*dN12drho + invJ12*dN12deta + invJ13*dN12dzeta;
			dN.entry[12][2] = invJ21*dN12drho + invJ22*dN12deta + invJ23*dN12dzeta;
			dN.entry[12][3] = invJ31*dN12drho + invJ32*dN12deta + invJ33*dN12dzeta;
			dN.entry[13][1] = invJ11*dN13drho + invJ12*dN13deta + invJ13*dN13dzeta;
			dN.entry[13][2] = invJ21*dN13drho + invJ22*dN13deta + invJ23*dN13dzeta;
			dN.entry[13][3] = invJ31*dN13drho + invJ32*dN13deta + invJ33*dN13dzeta;
			dN.entry[14][1] = invJ11*dN14drho + invJ12*dN14deta + invJ13*dN14dzeta;
			dN.entry[14][2] = invJ21*dN14drho + invJ22*dN14deta + invJ23*dN14dzeta;
			dN.entry[14][3] = invJ31*dN14drho + invJ32*dN14deta + invJ33*dN14dzeta;
			dN.entry[15][1] = invJ11*dN15drho + invJ12*dN15deta + invJ13*dN15dzeta;
			dN.entry[15][2] = invJ21*dN15drho + invJ22*dN15deta + invJ23*dN15dzeta;
			dN.entry[15][3] = invJ31*dN15drho + invJ32*dN15deta + invJ33*dN15dzeta;
			dN.entry[16][1] = invJ11*dN16drho + invJ12*dN16deta + invJ13*dN16dzeta;
			dN.entry[16][2] = invJ21*dN16drho + invJ22*dN16deta + invJ23*dN16dzeta;
			dN.entry[16][3] = invJ31*dN16drho + invJ32*dN16deta + invJ33*dN16dzeta;
			dN.entry[17][1] = invJ11*dN17drho + invJ12*dN17deta + invJ13*dN17dzeta;
			dN.entry[17][2] = invJ21*dN17drho + invJ22*dN17deta + invJ23*dN17dzeta;
			dN.entry[17][3] = invJ31*dN17drho + invJ32*dN17deta + invJ33*dN17dzeta;
			dN.entry[18][1] = invJ11*dN18drho + invJ12*dN18deta + invJ13*dN18dzeta;
			dN.entry[18][2] = invJ21*dN18drho + invJ22*dN18deta + invJ23*dN18dzeta;
			dN.entry[18][3] = invJ31*dN18drho + invJ32*dN18deta + invJ33*dN18dzeta;
			dN.entry[19][1] = invJ11*dN19drho + invJ12*dN19deta + invJ13*dN19dzeta;
			dN.entry[19][2] = invJ21*dN19drho + invJ22*dN19deta + invJ23*dN19dzeta;
			dN.entry[19][3] = invJ31*dN19drho + invJ32*dN19deta + invJ33*dN19dzeta;
			dN.entry[20][1] = invJ11*dN20drho + invJ12*dN20deta + invJ13*dN20dzeta;
			dN.entry[20][2] = invJ21*dN20drho + invJ22*dN20deta + invJ23*dN20dzeta;
			dN.entry[20][3] = invJ31*dN20drho + invJ32*dN20deta + invJ33*dN20dzeta;
		}


		void Hexahedron27(int element_id, const Vector<T>& point, Vector<T>& N, Matrix<T>& dN, T& det_J) const throw()
		{
			register T rho  = point.entry[1];
			register T eta  = point.entry[2];
			register T zeta = point.entry[3];

			// Shape functions
			N.entry[1]  =  ((eta - 1)*eta*(rho - 1)*rho*(zeta - 1)*zeta)*T(0.125);
			N.entry[2]  =  ((eta - 1)*eta*rho*(rho + 1)*(zeta - 1)*zeta)*T(0.125);
			N.entry[3]  =  (eta*(eta + 1)*rho*(rho + 1)*(zeta - 1)*zeta)*T(0.125);
			N.entry[4]  =  (eta*(eta + 1)*(rho - 1)*rho*(zeta - 1)*zeta)*T(0.125);
			N.entry[5]  =  ((eta - 1)*eta*(rho - 1)*rho*zeta*(zeta + 1))*T(0.125);
			N.entry[6]  =  ((eta - 1)*eta*rho*(rho + 1)*zeta*(zeta + 1))*T(0.125);
			N.entry[7]  =  (eta*(eta + 1)*rho*(rho + 1)*zeta*(zeta + 1))*T(0.125);
			N.entry[8]  =  (eta*(eta + 1)*(rho - 1)*rho*zeta*(zeta + 1))*T(0.125);
			N.entry[9]  = -((eta - 1)*eta*(rho - 1)*(rho + 1)*(zeta - 1)*zeta)*T(0.25);
			N.entry[10] = -((eta - 1)*(eta + 1)*rho*(rho + 1)*(zeta - 1)*zeta)*T(0.25);
			N.entry[11] = -(eta*(eta + 1)*(rho - 1)*(rho + 1)*(zeta - 1)*zeta)*T(0.25);
			N.entry[12] = -((eta - 1)*(eta + 1)*(rho - 1)*rho*(zeta - 1)*zeta)*T(0.25);
			N.entry[13] = -((eta - 1)*eta*(rho - 1)*rho*(zeta - 1)*(zeta + 1))*T(0.25);
			N.entry[14] = -((eta - 1)*eta*rho*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			N.entry[15] = -(eta*(eta + 1)*rho*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			N.entry[16] = -(eta*(eta + 1)*(rho - 1)*rho*(zeta - 1)*(zeta + 1))*T(0.25);
			N.entry[17] = -((eta - 1)*eta*(rho - 1)*(rho + 1)*zeta*(zeta + 1))*T(0.25);
			N.entry[18] = -((eta - 1)*(eta + 1)*rho*(rho + 1)*zeta*(zeta + 1))*T(0.25);
			N.entry[19] = -(eta*(eta + 1)*(rho - 1)*(rho + 1)*zeta*(zeta + 1))*T(0.25);
			N.entry[20] = -((eta - 1)*(eta + 1)*(rho - 1)*rho*zeta*(zeta + 1))*T(0.25);
			N.entry[21] =  ((eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*(zeta - 1)*zeta)*T(0.5);
			N.entry[22] =  ((eta - 1)*eta*(rho - 1)*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.5);
			N.entry[23] =  ((eta - 1)*(eta + 1)*rho*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.5);
			N.entry[24] =  (eta*(eta + 1)*(rho - 1)*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.5);
			N.entry[25] =  ((eta - 1)*(eta + 1)*(rho - 1)*rho*(zeta - 1)*(zeta + 1))*T(0.5);
			N.entry[26] =  ((eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*zeta*(zeta + 1))*T(0.5);
			N.entry[27] = -(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*(zeta - 1)*(zeta + 1);

			// Local coordinates derivative
			int node1  = mesh.connectivity.entry[element_id][1];
			int node2  = mesh.connectivity.entry[element_id][2];
			int node3  = mesh.connectivity.entry[element_id][3];
			int node4  = mesh.connectivity.entry[element_id][4];
			int node5  = mesh.connectivity.entry[element_id][5];
			int node6  = mesh.connectivity.entry[element_id][6];
			int node7  = mesh.connectivity.entry[element_id][7];
			int node8  = mesh.connectivity.entry[element_id][8];
			int node9  = mesh.connectivity.entry[element_id][9];
			int node10 = mesh.connectivity.entry[element_id][10];
			int node11 = mesh.connectivity.entry[element_id][11];
			int node12 = mesh.connectivity.entry[element_id][12];
			int node13 = mesh.connectivity.entry[element_id][13];
			int node14 = mesh.connectivity.entry[element_id][14];
			int node15 = mesh.connectivity.entry[element_id][15];
			int node16 = mesh.connectivity.entry[element_id][16];
			int node17 = mesh.connectivity.entry[element_id][17];
			int node18 = mesh.connectivity.entry[element_id][18];
			int node19 = mesh.connectivity.entry[element_id][19];
			int node20 = mesh.connectivity.entry[element_id][20];
			int node21 = mesh.connectivity.entry[element_id][21];
			int node22 = mesh.connectivity.entry[element_id][22];
			int node23 = mesh.connectivity.entry[element_id][23];
			int node24 = mesh.connectivity.entry[element_id][24];
			int node25 = mesh.connectivity.entry[element_id][25];
			int node26 = mesh.connectivity.entry[element_id][26];
			int node27 = mesh.connectivity.entry[element_id][27];
			T x1  = nodes.coordinate.entry[node1][1];
			T y1  = nodes.coordinate.entry[node1][2];
			T z1  = nodes.coordinate.entry[node1][3];
			T x2  = nodes.coordinate.entry[node2][1];
			T y2  = nodes.coordinate.entry[node2][2];
			T z2  = nodes.coordinate.entry[node2][3];
			T x3  = nodes.coordinate.entry[node3][1];
			T y3  = nodes.coordinate.entry[node3][2];
			T z3  = nodes.coordinate.entry[node3][3];
			T x4  = nodes.coordinate.entry[node4][1];
			T y4  = nodes.coordinate.entry[node4][2];
			T z4  = nodes.coordinate.entry[node4][3];
			T x5  = nodes.coordinate.entry[node5][1];
			T y5  = nodes.coordinate.entry[node5][2];
			T z5  = nodes.coordinate.entry[node5][3];
			T x6  = nodes.coordinate.entry[node6][1];
			T y6  = nodes.coordinate.entry[node6][2];
			T z6  = nodes.coordinate.entry[node6][3];
			T x7  = nodes.coordinate.entry[node7][1];
			T y7  = nodes.coordinate.entry[node7][2];
			T z7  = nodes.coordinate.entry[node7][3];
			T x8  = nodes.coordinate.entry[node8][1];
			T y8  = nodes.coordinate.entry[node8][2];
			T z8  = nodes.coordinate.entry[node8][3];
			T x9  = nodes.coordinate.entry[node9][1];
			T y9  = nodes.coordinate.entry[node9][2];
			T z9  = nodes.coordinate.entry[node9][3];
			T x10 = nodes.coordinate.entry[node10][1];
			T y10 = nodes.coordinate.entry[node10][2];
			T z10 = nodes.coordinate.entry[node10][3];
			T x11 = nodes.coordinate.entry[node11][1];
			T y11 = nodes.coordinate.entry[node11][2];
			T z11 = nodes.coordinate.entry[node11][3];
			T x12 = nodes.coordinate.entry[node12][1];
			T y12 = nodes.coordinate.entry[node12][2];
			T z12 = nodes.coordinate.entry[node12][3];
			T x13 = nodes.coordinate.entry[node13][1];
			T y13 = nodes.coordinate.entry[node13][2];
			T z13 = nodes.coordinate.entry[node13][3];
			T x14 = nodes.coordinate.entry[node14][1];
			T y14 = nodes.coordinate.entry[node14][2];
			T z14 = nodes.coordinate.entry[node14][3];
			T x15 = nodes.coordinate.entry[node15][1];
			T y15 = nodes.coordinate.entry[node15][2];
			T z15 = nodes.coordinate.entry[node15][3];
			T x16 = nodes.coordinate.entry[node16][1];
			T y16 = nodes.coordinate.entry[node16][2];
			T z16 = nodes.coordinate.entry[node16][3];
			T x17 = nodes.coordinate.entry[node17][1];
			T y17 = nodes.coordinate.entry[node17][2];
			T z17 = nodes.coordinate.entry[node17][3];
			T x18 = nodes.coordinate.entry[node18][1];
			T y18 = nodes.coordinate.entry[node18][2];
			T z18 = nodes.coordinate.entry[node18][3];
			T x19 = nodes.coordinate.entry[node19][1];
			T y19 = nodes.coordinate.entry[node19][2];
			T z19 = nodes.coordinate.entry[node19][3];
			T x20 = nodes.coordinate.entry[node20][1];
			T y20 = nodes.coordinate.entry[node20][2];
			T z20 = nodes.coordinate.entry[node20][3];
			T x21 = nodes.coordinate.entry[node21][1];
			T y21 = nodes.coordinate.entry[node21][2];
			T z21 = nodes.coordinate.entry[node21][3];
			T x22 = nodes.coordinate.entry[node22][1];
			T y22 = nodes.coordinate.entry[node22][2];
			T z22 = nodes.coordinate.entry[node22][3];
			T x23 = nodes.coordinate.entry[node23][1];
			T y23 = nodes.coordinate.entry[node23][2];
			T z23 = nodes.coordinate.entry[node23][3];
			T x24 = nodes.coordinate.entry[node24][1];
			T y24 = nodes.coordinate.entry[node24][2];
			T z24 = nodes.coordinate.entry[node24][3];
			T x25 = nodes.coordinate.entry[node25][1];
			T y25 = nodes.coordinate.entry[node25][2];
			T z25 = nodes.coordinate.entry[node25][3];
			T x26 = nodes.coordinate.entry[node26][1];
			T y26 = nodes.coordinate.entry[node26][2];
			T z26 = nodes.coordinate.entry[node26][3];
			T x27 = nodes.coordinate.entry[node27][1];
			T y27 = nodes.coordinate.entry[node27][2];
			T z27 = nodes.coordinate.entry[node27][3];

			T dxdrho  = (eta*(eta + 1)*(2*rho - 1)*x8*zeta*(zeta + 1) + eta*(eta + 1)*(2*rho + 1)*x7*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho + 1)*x6*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho - 1)*x5*zeta*(zeta + 1) + 8*(eta - 1)*(eta + 1)*rho*x26*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*x20*zeta*(zeta + 1) - 4*eta*(eta + 1)*rho*x19*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*x18*zeta*(zeta + 1) - 4*(eta - 1)*eta*rho*x17*zeta*(zeta + 1) - 16*(eta - 1)*(eta + 1)*rho*x27*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho - 1)*x25*(zeta - 1)*(zeta + 1) + 8*eta*(eta + 1)*rho*x24*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho + 1)*x23*(zeta - 1)*(zeta + 1) + 8*(eta - 1)*eta*rho*x22*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho - 1)*x16*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho + 1)*x15*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho + 1)*x14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho - 1)*x13*(zeta - 1)*(zeta + 1) - 4*(eta - 1)*eta*rho*x9*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho - 1)*x4*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho + 1)*x3*(zeta - 1)*zeta + 8*(eta - 1)*(eta + 1)*rho*x21*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho + 1)*x2*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*x12*(zeta - 1)*zeta - 4*eta*(eta + 1)*rho*x11*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*x10*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho - 1)*x1*(zeta - 1)*zeta)*T(0.125);
			T dxdeta  = ((2*eta + 1)*(rho - 1)*rho*x8*zeta*(zeta + 1) + (2*eta + 1)*rho*(rho + 1)*x7*zeta*(zeta + 1) + (2*eta - 1)*rho*(rho + 1)*x6*zeta*(zeta + 1) + (2*eta - 1)*(rho - 1)*rho*x5*zeta*(zeta + 1) + 8*eta*(rho - 1)*(rho + 1)*x26*zeta*(zeta + 1) - 4*eta*(rho - 1)*rho*x20*zeta*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*x19*zeta*(zeta + 1) - 4*eta*rho*(rho + 1)*x18*zeta*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*x17*zeta*(zeta + 1) - 16*eta*(rho - 1)*(rho + 1)*x27*(zeta - 1)*(zeta + 1) + 8*eta*(rho - 1)*rho*x25*(zeta - 1)*(zeta + 1) + 4*(2*eta + 1)*(rho - 1)*(rho + 1)*x24*(zeta - 1)*(zeta + 1) + 8*eta*rho*(rho + 1)*x23*(zeta - 1)*(zeta + 1) + 4*(2*eta - 1)*(rho - 1)*(rho + 1)*x22*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*rho*x16*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*rho*(rho + 1)*x15*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*rho*(rho + 1)*x14*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*rho*x13*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*x9*(zeta - 1)*zeta + (2*eta + 1)*(rho - 1)*rho*x4*(zeta - 1)*zeta + (2*eta + 1)*rho*(rho + 1)*x3*(zeta - 1)*zeta + 8*eta*(rho - 1)*(rho + 1)*x21*(zeta - 1)*zeta + (2*eta - 1)*rho*(rho + 1)*x2*(zeta - 1)*zeta - 4*eta*(rho - 1)*rho*x12*(zeta - 1)*zeta - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*x11*(zeta - 1)*zeta - 4*eta*rho*(rho + 1)*x10*(zeta - 1)*zeta + (2*eta - 1)*(rho - 1)*rho*x1*(zeta - 1)*zeta)*T(0.125);
			T dxdzeta = (eta*(eta + 1)*(rho - 1)*rho*x8*(2*zeta + 1) + eta*(eta + 1)*rho*(rho + 1)*x7*(2*zeta + 1) + (eta - 1)*eta*rho*(rho + 1)*x6*(2*zeta + 1) + (eta - 1)*eta*(rho - 1)*rho*x5*(2*zeta + 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*x26*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*x20*(2*zeta + 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*x19*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*x18*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*x17*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*x9*(2*zeta - 1) + eta*(eta + 1)*(rho - 1)*rho*x4*(2*zeta - 1) + eta*(eta + 1)*rho*(rho + 1)*x3*(2*zeta - 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*x21*(2*zeta - 1) + (eta - 1)*eta*rho*(rho + 1)*x2*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*x12*(2*zeta - 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*x11*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*x10*(2*zeta - 1) + (eta - 1)*eta*(rho - 1)*rho*x1*(2*zeta - 1) - 16*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*x27*zeta + 8*(eta - 1)*(eta + 1)*(rho - 1)*rho*x25*zeta + 8*eta*(eta + 1)*(rho - 1)*(rho + 1)*x24*zeta + 8*(eta - 1)*(eta + 1)*rho*(rho + 1)*x23*zeta + 8*(eta - 1)*eta*(rho - 1)*(rho + 1)*x22*zeta - 4*eta*(eta + 1)*(rho - 1)*rho*x16*zeta - 4*eta*(eta + 1)*rho*(rho + 1)*x15*zeta - 4*(eta - 1)*eta*rho*(rho + 1)*x14*zeta - 4*(eta - 1)*eta*(rho - 1)*rho*x13*zeta)*T(0.125);
			T dydrho  = (eta*(eta + 1)*(2*rho - 1)*y8*zeta*(zeta + 1) + eta*(eta + 1)*(2*rho + 1)*y7*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho + 1)*y6*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho - 1)*y5*zeta*(zeta + 1) + 8*(eta - 1)*(eta + 1)*rho*y26*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*y20*zeta*(zeta + 1) - 4*eta*(eta + 1)*rho*y19*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*y18*zeta*(zeta + 1) - 4*(eta - 1)*eta*rho*y17*zeta*(zeta + 1) - 16*(eta - 1)*(eta + 1)*rho*y27*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho - 1)*y25*(zeta - 1)*(zeta + 1) + 8*eta*(eta + 1)*rho*y24*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho + 1)*y23*(zeta - 1)*(zeta + 1) + 8*(eta - 1)*eta*rho*y22*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho - 1)*y16*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho + 1)*y15*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho + 1)*y14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho - 1)*y13*(zeta - 1)*(zeta + 1) - 4*(eta - 1)*eta*rho*y9*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho - 1)*y4*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho + 1)*y3*(zeta - 1)*zeta + 8*(eta - 1)*(eta + 1)*rho*y21*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho + 1)*y2*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*y12*(zeta - 1)*zeta - 4*eta*(eta + 1)*rho*y11*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*y10*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho - 1)*y1*(zeta - 1)*zeta)*T(0.125);
			T dydeta  = ((2*eta + 1)*(rho - 1)*rho*y8*zeta*(zeta + 1) + (2*eta + 1)*rho*(rho + 1)*y7*zeta*(zeta + 1) + (2*eta - 1)*rho*(rho + 1)*y6*zeta*(zeta + 1) + (2*eta - 1)*(rho - 1)*rho*y5*zeta*(zeta + 1) + 8*eta*(rho - 1)*(rho + 1)*y26*zeta*(zeta + 1) - 4*eta*(rho - 1)*rho*y20*zeta*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*y19*zeta*(zeta + 1) - 4*eta*rho*(rho + 1)*y18*zeta*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*y17*zeta*(zeta + 1) - 16*eta*(rho - 1)*(rho + 1)*y27*(zeta - 1)*(zeta + 1) + 8*eta*(rho - 1)*rho*y25*(zeta - 1)*(zeta + 1) + 4*(2*eta + 1)*(rho - 1)*(rho + 1)*y24*(zeta - 1)*(zeta + 1) + 8*eta*rho*(rho + 1)*y23*(zeta - 1)*(zeta + 1) + 4*(2*eta - 1)*(rho - 1)*(rho + 1)*y22*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*rho*y16*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*rho*(rho + 1)*y15*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*rho*(rho + 1)*y14*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*rho*y13*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*y9*(zeta - 1)*zeta + (2*eta + 1)*(rho - 1)*rho*y4*(zeta - 1)*zeta + (2*eta + 1)*rho*(rho + 1)*y3*(zeta - 1)*zeta + 8*eta*(rho - 1)*(rho + 1)*y21*(zeta - 1)*zeta + (2*eta - 1)*rho*(rho + 1)*y2*(zeta - 1)*zeta - 4*eta*(rho - 1)*rho*y12*(zeta - 1)*zeta - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*y11*(zeta - 1)*zeta - 4*eta*rho*(rho + 1)*y10*(zeta - 1)*zeta + (2*eta - 1)*(rho - 1)*rho*y1*(zeta - 1)*zeta)*T(0.125);
			T dydzeta = (eta*(eta + 1)*(rho - 1)*rho*y8*(2*zeta + 1) + eta*(eta + 1)*rho*(rho + 1)*y7*(2*zeta + 1) + (eta - 1)*eta*rho*(rho + 1)*y6*(2*zeta + 1) + (eta - 1)*eta*(rho - 1)*rho*y5*(2*zeta + 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*y26*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*y20*(2*zeta + 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*y19*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*y18*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*y17*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*y9*(2*zeta - 1) + eta*(eta + 1)*(rho - 1)*rho*y4*(2*zeta - 1) + eta*(eta + 1)*rho*(rho + 1)*y3*(2*zeta - 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*y21*(2*zeta - 1) + (eta - 1)*eta*rho*(rho + 1)*y2*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*y12*(2*zeta - 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*y11*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*y10*(2*zeta - 1) + (eta - 1)*eta*(rho - 1)*rho*y1*(2*zeta - 1) - 16*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*y27*zeta + 8*(eta - 1)*(eta + 1)*(rho - 1)*rho*y25*zeta + 8*eta*(eta + 1)*(rho - 1)*(rho + 1)*y24*zeta + 8*(eta - 1)*(eta + 1)*rho*(rho + 1)*y23*zeta + 8*(eta - 1)*eta*(rho - 1)*(rho + 1)*y22*zeta - 4*eta*(eta + 1)*(rho - 1)*rho*y16*zeta - 4*eta*(eta + 1)*rho*(rho + 1)*y15*zeta - 4*(eta - 1)*eta*rho*(rho + 1)*y14*zeta - 4*(eta - 1)*eta*(rho - 1)*rho*y13*zeta)*T(0.125);
			T dzdrho  = (eta*(eta + 1)*(2*rho - 1)*z8*zeta*(zeta + 1) + eta*(eta + 1)*(2*rho + 1)*z7*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho + 1)*z6*zeta*(zeta + 1) + (eta - 1)*eta*(2*rho - 1)*z5*zeta*(zeta + 1) + 8*(eta - 1)*(eta + 1)*rho*z26*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*z20*zeta*(zeta + 1) - 4*eta*(eta + 1)*rho*z19*zeta*(zeta + 1) - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*z18*zeta*(zeta + 1) - 4*(eta - 1)*eta*rho*z17*zeta*(zeta + 1) - 16*(eta - 1)*(eta + 1)*rho*z27*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho - 1)*z25*(zeta - 1)*(zeta + 1) + 8*eta*(eta + 1)*rho*z24*(zeta - 1)*(zeta + 1) + 4*(eta - 1)*(eta + 1)*(2*rho + 1)*z23*(zeta - 1)*(zeta + 1) + 8*(eta - 1)*eta*rho*z22*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho - 1)*z16*(zeta - 1)*(zeta + 1) - 2*eta*(eta + 1)*(2*rho + 1)*z15*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho + 1)*z14*(zeta - 1)*(zeta + 1) - 2*(eta - 1)*eta*(2*rho - 1)*z13*(zeta - 1)*(zeta + 1) - 4*(eta - 1)*eta*rho*z9*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho - 1)*z4*(zeta - 1)*zeta + eta*(eta + 1)*(2*rho + 1)*z3*(zeta - 1)*zeta + 8*(eta - 1)*(eta + 1)*rho*z21*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho + 1)*z2*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho - 1)*z12*(zeta - 1)*zeta - 4*eta*(eta + 1)*rho*z11*(zeta - 1)*zeta - 2*(eta - 1)*(eta + 1)*(2*rho + 1)*z10*(zeta - 1)*zeta + (eta - 1)*eta*(2*rho - 1)*z1*(zeta - 1)*zeta)*T(0.125);
			T dzdeta  = ((2*eta + 1)*(rho - 1)*rho*z8*zeta*(zeta + 1) + (2*eta + 1)*rho*(rho + 1)*z7*zeta*(zeta + 1) + (2*eta - 1)*rho*(rho + 1)*z6*zeta*(zeta + 1) + (2*eta - 1)*(rho - 1)*rho*z5*zeta*(zeta + 1) + 8*eta*(rho - 1)*(rho + 1)*z26*zeta*(zeta + 1) - 4*eta*(rho - 1)*rho*z20*zeta*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*z19*zeta*(zeta + 1) - 4*eta*rho*(rho + 1)*z18*zeta*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*z17*zeta*(zeta + 1) - 16*eta*(rho - 1)*(rho + 1)*z27*(zeta - 1)*(zeta + 1) + 8*eta*(rho - 1)*rho*z25*(zeta - 1)*(zeta + 1) + 4*(2*eta + 1)*(rho - 1)*(rho + 1)*z24*(zeta - 1)*(zeta + 1) + 8*eta*rho*(rho + 1)*z23*(zeta - 1)*(zeta + 1) + 4*(2*eta - 1)*(rho - 1)*(rho + 1)*z22*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*(rho - 1)*rho*z16*(zeta - 1)*(zeta + 1) - 2*(2*eta + 1)*rho*(rho + 1)*z15*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*rho*(rho + 1)*z14*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*rho*z13*(zeta - 1)*(zeta + 1) - 2*(2*eta - 1)*(rho - 1)*(rho + 1)*z9*(zeta - 1)*zeta + (2*eta + 1)*(rho - 1)*rho*z4*(zeta - 1)*zeta + (2*eta + 1)*rho*(rho + 1)*z3*(zeta - 1)*zeta + 8*eta*(rho - 1)*(rho + 1)*z21*(zeta - 1)*zeta + (2*eta - 1)*rho*(rho + 1)*z2*(zeta - 1)*zeta - 4*eta*(rho - 1)*rho*z12*(zeta - 1)*zeta - 2*(2*eta + 1)*(rho - 1)*(rho + 1)*z11*(zeta - 1)*zeta - 4*eta*rho*(rho + 1)*z10*(zeta - 1)*zeta + (2*eta - 1)*(rho - 1)*rho*z1*(zeta - 1)*zeta)*T(0.125);
			T dzdzeta = (eta*(eta + 1)*(rho - 1)*rho*z8*(2*zeta + 1) + eta*(eta + 1)*rho*(rho + 1)*z7*(2*zeta + 1) + (eta - 1)*eta*rho*(rho + 1)*z6*(2*zeta + 1) + (eta - 1)*eta*(rho - 1)*rho*z5*(2*zeta + 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*z26*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*z20*(2*zeta + 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*z19*(2*zeta + 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*z18*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*z17*(2*zeta + 1) - 2*(eta - 1)*eta*(rho - 1)*(rho + 1)*z9*(2*zeta - 1) + eta*(eta + 1)*(rho - 1)*rho*z4*(2*zeta - 1) + eta*(eta + 1)*rho*(rho + 1)*z3*(2*zeta - 1) + 4*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*z21*(2*zeta - 1) + (eta - 1)*eta*rho*(rho + 1)*z2*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*(rho - 1)*rho*z12*(2*zeta - 1) - 2*eta*(eta + 1)*(rho - 1)*(rho + 1)*z11*(2*zeta - 1) - 2*(eta - 1)*(eta + 1)*rho*(rho + 1)*z10*(2*zeta - 1) + (eta - 1)*eta*(rho - 1)*rho*z1*(2*zeta - 1) - 16*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*z27*zeta + 8*(eta - 1)*(eta + 1)*(rho - 1)*rho*z25*zeta + 8*eta*(eta + 1)*(rho - 1)*(rho + 1)*z24*zeta + 8*(eta - 1)*(eta + 1)*rho*(rho + 1)*z23*zeta + 8*(eta - 1)*eta*(rho - 1)*(rho + 1)*z22*zeta - 4*eta*(eta + 1)*(rho - 1)*rho*z16*zeta - 4*eta*(eta + 1)*rho*(rho + 1)*z15*zeta - 4*(eta - 1)*eta*rho*(rho + 1)*z14*zeta - 4*(eta - 1)*eta*(rho - 1)*rho*z13*zeta)*T(0.125);

			// Jacobian determinant
			det_J = dxdrho*(dydeta*dzdzeta - dydzeta*dzdeta) - dydrho*(dxdeta*dzdzeta - dxdzeta*dzdeta) + (dxdeta*dydzeta - dxdzeta*dydeta)*dzdrho;

			// Inverse of the Jacobian
			T invJ11 = (dydeta*dzdzeta - dydzeta*dzdeta)/det_J;
			T invJ12 = (dydzeta*dzdrho - dydrho*dzdzeta)/det_J;
			T invJ13 = (dydrho*dzdeta - dydeta*dzdrho)/det_J;
			T invJ21 = (dxdzeta*dzdeta - dxdeta*dzdzeta)/det_J;
			T invJ22 = (dxdrho*dzdzeta - dxdzeta*dzdrho)/det_J;
			T invJ23 = (dxdeta*dzdrho - dxdrho*dzdeta)/det_J;
			T invJ31 = (dxdeta*dydzeta - dxdzeta*dydeta)/det_J;
			T invJ32 = (dxdzeta*dydrho - dxdrho*dydzeta)/det_J;
			T invJ33 = (dxdrho*dydeta - dxdeta*dydrho)/det_J;

			// Global coordinates derivative
			T dN1drho   =  ((eta - 1)*eta*(2*rho - 1)*(zeta - 1)*zeta)*T(0.125);
			T dN1deta   =  ((2*eta - 1)*(rho - 1)*rho*(zeta - 1)*zeta)*T(0.125);
			T dN1dzeta  =  ((eta - 1)*eta*(rho - 1)*rho*(2*zeta - 1))*T(0.125);
			T dN2drho   =  ((eta - 1)*eta*(2*rho + 1)*(zeta - 1)*zeta)*T(0.125);
			T dN2deta   =  ((2*eta - 1)*rho*(rho + 1)*(zeta - 1)*zeta)*T(0.125);
			T dN2dzeta  =  ((eta - 1)*eta*rho*(rho + 1)*(2*zeta - 1))*T(0.125);
			T dN3drho   =  (eta*(eta + 1)*(2*rho + 1)*(zeta - 1)*zeta)*T(0.125);
			T dN3deta   =  ((2*eta + 1)*rho*(rho + 1)*(zeta - 1)*zeta)*T(0.125);
			T dN3dzeta  =  (eta*(eta + 1)*rho*(rho + 1)*(2*zeta - 1))*T(0.125);
			T dN4drho   =  (eta*(eta + 1)*(2*rho - 1)*(zeta - 1)*zeta)*T(0.125);
			T dN4deta   =  ((2*eta + 1)*(rho - 1)*rho*(zeta - 1)*zeta)*T(0.125);
			T dN4dzeta  =  (eta*(eta + 1)*(rho - 1)*rho*(2*zeta - 1))*T(0.125);
			T dN5drho   =  ((eta - 1)*eta*(2*rho - 1)*zeta*(zeta + 1))*T(0.125);
			T dN5deta   =  ((2*eta - 1)*(rho - 1)*rho*zeta*(zeta + 1))*T(0.125);
			T dN5dzeta  =  ((eta - 1)*eta*(rho - 1)*rho*(2*zeta + 1))*T(0.125);
			T dN6drho   =  ((eta - 1)*eta*(2*rho + 1)*zeta*(zeta + 1))*T(0.125);
			T dN6deta   =  ((2*eta - 1)*rho*(rho + 1)*zeta*(zeta + 1))*T(0.125);
			T dN6dzeta  =  ((eta - 1)*eta*rho*(rho + 1)*(2*zeta + 1))*T(0.125);
			T dN7drho   =  (eta*(eta + 1)*(2*rho + 1)*zeta*(zeta + 1))*T(0.125);
			T dN7deta   =  ((2*eta + 1)*rho*(rho + 1)*zeta*(zeta + 1))*T(0.125);
			T dN7dzeta  =  (eta*(eta + 1)*rho*(rho + 1)*(2*zeta + 1))*T(0.125);
			T dN8drho   =  (eta*(eta + 1)*(2*rho - 1)*zeta*(zeta + 1))*T(0.125);
			T dN8deta   =  ((2*eta + 1)*(rho - 1)*rho*zeta*(zeta + 1))*T(0.125);
			T dN8dzeta  =  (eta*(eta + 1)*(rho - 1)*rho*(2*zeta + 1))*T(0.125);
			T dN9drho   = -((eta - 1)*eta*rho*(zeta - 1)*zeta)*T(0.5);
			T dN9deta   = -((2*eta - 1)*(rho - 1)*(rho + 1)*(zeta - 1)*zeta)*T(0.25);
			T dN9dzeta  = -((eta - 1)*eta*(rho - 1)*(rho + 1)*(2*zeta - 1))*T(0.25);
			T dN10drho  = -((eta - 1)*(eta + 1)*(2*rho + 1)*(zeta - 1)*zeta)*T(0.25);
			T dN10deta  = -(eta*rho*(rho + 1)*(zeta - 1)*zeta)*T(0.5);
			T dN10dzeta = -((eta - 1)*(eta + 1)*rho*(rho + 1)*(2*zeta - 1))*T(0.25);
			T dN11drho  = -(eta*(eta + 1)*rho*(zeta - 1)*zeta)*T(0.5);
			T dN11deta  = -((2*eta + 1)*(rho - 1)*(rho + 1)*(zeta - 1)*zeta)*T(0.25);
			T dN11dzeta = -(eta*(eta + 1)*(rho - 1)*(rho + 1)*(2*zeta - 1))*T(0.25);
			T dN12drho  = -((eta - 1)*(eta + 1)*(2*rho - 1)*(zeta - 1)*zeta)*T(0.25);
			T dN12deta  = -(eta*(rho - 1)*rho*(zeta - 1)*zeta)*T(0.5);
			T dN12dzeta = -((eta - 1)*(eta + 1)*(rho - 1)*rho*(2*zeta - 1))*T(0.25);
			T dN13drho  = -((eta - 1)*eta*(2*rho - 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN13deta  = -((2*eta - 1)*(rho - 1)*rho*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN13dzeta = -((eta - 1)*eta*(rho - 1)*rho*zeta)*T(0.5);
			T dN14drho  = -((eta - 1)*eta*(2*rho + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN14deta  = -((2*eta - 1)*rho*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN14dzeta = -((eta - 1)*eta*rho*(rho + 1)*zeta)*T(0.5);
			T dN15drho  = -(eta*(eta + 1)*(2*rho + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN15deta  = -((2*eta + 1)*rho*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN15dzeta = -(eta*(eta + 1)*rho*(rho + 1)*zeta)*T(0.5);
			T dN16drho  = -(eta*(eta + 1)*(2*rho - 1)*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN16deta  = -((2*eta + 1)*(rho - 1)*rho*(zeta - 1)*(zeta + 1))*T(0.25);
			T dN16dzeta = -(eta*(eta + 1)*(rho - 1)*rho*zeta)*T(0.5);
			T dN17drho  = -((eta - 1)*eta*rho*zeta*(zeta + 1))*T(0.5);
			T dN17deta  = -((2*eta - 1)*(rho - 1)*(rho + 1)*zeta*(zeta + 1))*T(0.25);
			T dN17dzeta = -((eta - 1)*eta*(rho - 1)*(rho + 1)*(2*zeta + 1))*T(0.25);
			T dN18drho  = -((eta - 1)*(eta + 1)*(2*rho + 1)*zeta*(zeta + 1))*T(0.25);
			T dN18deta  = -(eta*rho*(rho + 1)*zeta*(zeta + 1))*T(0.5);
			T dN18dzeta = -((eta - 1)*(eta + 1)*rho*(rho + 1)*(2*zeta + 1))*T(0.25);
			T dN19drho  = -(eta*(eta + 1)*rho*zeta*(zeta + 1))*T(0.5);
			T dN19deta  = -((2*eta + 1)*(rho - 1)*(rho + 1)*zeta*(zeta + 1))*T(0.25);
			T dN19dzeta = -(eta*(eta + 1)*(rho - 1)*(rho + 1)*(2*zeta + 1))*T(0.25);
			T dN20drho  = -((eta - 1)*(eta + 1)*(2*rho - 1)*zeta*(zeta + 1))*T(0.25);
			T dN20deta  = -(eta*(rho - 1)*rho*zeta*(zeta + 1))*T(0.5);
			T dN20dzeta = -((eta - 1)*(eta + 1)*(rho - 1)*rho*(2*zeta + 1))*T(0.25);
			T dN21drho  =  (eta - 1)*(eta + 1)*rho*(zeta - 1)*zeta;
			T dN21deta  =  eta*(rho - 1)*(rho + 1)*(zeta - 1)*zeta;
			T dN21dzeta =  ((eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*(2*zeta - 1))*T(0.5);
			T dN22drho  =  (eta - 1)*eta*rho*(zeta - 1)*(zeta + 1);
			T dN22deta  =  ((2*eta - 1)*(rho - 1)*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.5);
			T dN22dzeta =  (eta - 1)*eta*(rho - 1)*(rho + 1)*zeta;
			T dN23drho  =  ((eta - 1)*(eta + 1)*(2*rho + 1)*(zeta - 1)*(zeta + 1))*T(0.5);
			T dN23deta  =  eta*rho*(rho + 1)*(zeta - 1)*(zeta + 1);
			T dN23dzeta =  (eta - 1)*(eta + 1)*rho*(rho + 1)*zeta;
			T dN24drho  =  eta*(eta + 1)*rho*(zeta - 1)*(zeta + 1);
			T dN24deta  =  ((2*eta + 1)*(rho - 1)*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.5);
			T dN24dzeta =  eta*(eta + 1)*(rho - 1)*(rho + 1)*zeta;
			T dN25drho  =  ((eta - 1)*(eta + 1)*(2*rho - 1)*(zeta - 1)*(zeta + 1))*T(0.5);
			T dN25deta  =  eta*(rho - 1)*rho*(zeta - 1)*(zeta + 1);
			T dN25dzeta =  (eta - 1)*(eta + 1)*(rho - 1)*rho*zeta;
			T dN26drho  =  (eta - 1)*(eta + 1)*rho*zeta*(zeta + 1);
			T dN26deta  =  eta*(rho - 1)*(rho + 1)*zeta*(zeta + 1);
			T dN26dzeta =  ((eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*(2*zeta + 1))*T(0.5);
			T dN27drho  = -2*(eta - 1)*(eta + 1)*rho*(zeta - 1)*(zeta + 1);
			T dN27deta  = -2*eta*(rho - 1)*(rho + 1)*(zeta - 1)*(zeta + 1);
			T dN27dzeta = -2*(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*zeta;

			dN.entry[1][1]  = invJ11*dN1drho + invJ12*dN1deta + invJ13*dN1dzeta;
			dN.entry[1][2]  = invJ21*dN1drho + invJ22*dN1deta + invJ23*dN1dzeta;
			dN.entry[1][3]  = invJ31*dN1drho + invJ32*dN1deta + invJ33*dN1dzeta;
			dN.entry[2][1]  = invJ11*dN2drho + invJ12*dN2deta + invJ13*dN2dzeta;
			dN.entry[2][2]  = invJ21*dN2drho + invJ22*dN2deta + invJ23*dN2dzeta;
			dN.entry[2][3]  = invJ31*dN2drho + invJ32*dN2deta + invJ33*dN2dzeta;
			dN.entry[3][1]  = invJ11*dN3drho + invJ12*dN3deta + invJ13*dN3dzeta;
			dN.entry[3][2]  = invJ21*dN3drho + invJ22*dN3deta + invJ23*dN3dzeta;
			dN.entry[3][3]  = invJ31*dN3drho + invJ32*dN3deta + invJ33*dN3dzeta;
			dN.entry[4][1]  = invJ11*dN4drho + invJ12*dN4deta + invJ13*dN4dzeta;
			dN.entry[4][2]  = invJ21*dN4drho + invJ22*dN4deta + invJ23*dN4dzeta;
			dN.entry[4][3]  = invJ31*dN4drho + invJ32*dN4deta + invJ33*dN4dzeta;
			dN.entry[5][1]  = invJ11*dN5drho + invJ12*dN5deta + invJ13*dN5dzeta;
			dN.entry[5][2]  = invJ21*dN5drho + invJ22*dN5deta + invJ23*dN5dzeta;
			dN.entry[5][3]  = invJ31*dN5drho + invJ32*dN5deta + invJ33*dN5dzeta;
			dN.entry[6][1]  = invJ11*dN6drho + invJ12*dN6deta + invJ13*dN6dzeta;
			dN.entry[6][2]  = invJ21*dN6drho + invJ22*dN6deta + invJ23*dN6dzeta;
			dN.entry[6][3]  = invJ31*dN6drho + invJ32*dN6deta + invJ33*dN6dzeta;
			dN.entry[7][1]  = invJ11*dN7drho + invJ12*dN7deta + invJ13*dN7dzeta;
			dN.entry[7][2]  = invJ21*dN7drho + invJ22*dN7deta + invJ23*dN7dzeta;
			dN.entry[7][3]  = invJ31*dN7drho + invJ32*dN7deta + invJ33*dN7dzeta;
			dN.entry[8][1]  = invJ11*dN8drho + invJ12*dN8deta + invJ13*dN8dzeta;
			dN.entry[8][2]  = invJ21*dN8drho + invJ22*dN8deta + invJ23*dN8dzeta;
			dN.entry[8][3]  = invJ31*dN8drho + invJ32*dN8deta + invJ33*dN8dzeta;
			dN.entry[9][1]  = invJ11*dN9drho + invJ12*dN9deta + invJ13*dN9dzeta;
			dN.entry[9][2]  = invJ21*dN9drho + invJ22*dN9deta + invJ23*dN9dzeta;
			dN.entry[9][3]  = invJ31*dN9drho + invJ32*dN9deta + invJ33*dN9dzeta;
			dN.entry[10][1] = invJ11*dN10drho + invJ12*dN10deta + invJ13*dN10dzeta;
			dN.entry[10][2] = invJ21*dN10drho + invJ22*dN10deta + invJ23*dN10dzeta;
			dN.entry[10][3] = invJ31*dN10drho + invJ32*dN10deta + invJ33*dN10dzeta;
			dN.entry[11][1] = invJ11*dN11drho + invJ12*dN11deta + invJ13*dN11dzeta;
			dN.entry[11][2] = invJ21*dN11drho + invJ22*dN11deta + invJ23*dN11dzeta;
			dN.entry[11][3] = invJ31*dN11drho + invJ32*dN11deta + invJ33*dN11dzeta;
			dN.entry[12][1] = invJ11*dN12drho + invJ12*dN12deta + invJ13*dN12dzeta;
			dN.entry[12][2] = invJ21*dN12drho + invJ22*dN12deta + invJ23*dN12dzeta;
			dN.entry[12][3] = invJ31*dN12drho + invJ32*dN12deta + invJ33*dN12dzeta;
			dN.entry[13][1] = invJ11*dN13drho + invJ12*dN13deta + invJ13*dN13dzeta;
			dN.entry[13][2] = invJ21*dN13drho + invJ22*dN13deta + invJ23*dN13dzeta;
			dN.entry[13][3] = invJ31*dN13drho + invJ32*dN13deta + invJ33*dN13dzeta;
			dN.entry[14][1] = invJ11*dN14drho + invJ12*dN14deta + invJ13*dN14dzeta;
			dN.entry[14][2] = invJ21*dN14drho + invJ22*dN14deta + invJ23*dN14dzeta;
			dN.entry[14][3] = invJ31*dN14drho + invJ32*dN14deta + invJ33*dN14dzeta;
			dN.entry[15][1] = invJ11*dN15drho + invJ12*dN15deta + invJ13*dN15dzeta;
			dN.entry[15][2] = invJ21*dN15drho + invJ22*dN15deta + invJ23*dN15dzeta;
			dN.entry[15][3] = invJ31*dN15drho + invJ32*dN15deta + invJ33*dN15dzeta;
			dN.entry[16][1] = invJ11*dN16drho + invJ12*dN16deta + invJ13*dN16dzeta;
			dN.entry[16][2] = invJ21*dN16drho + invJ22*dN16deta + invJ23*dN16dzeta;
			dN.entry[16][3] = invJ31*dN16drho + invJ32*dN16deta + invJ33*dN16dzeta;
			dN.entry[17][1] = invJ11*dN17drho + invJ12*dN17deta + invJ13*dN17dzeta;
			dN.entry[17][2] = invJ21*dN17drho + invJ22*dN17deta + invJ23*dN17dzeta;
			dN.entry[17][3] = invJ31*dN17drho + invJ32*dN17deta + invJ33*dN17dzeta;
			dN.entry[18][1] = invJ11*dN18drho + invJ12*dN18deta + invJ13*dN18dzeta;
			dN.entry[18][2] = invJ21*dN18drho + invJ22*dN18deta + invJ23*dN18dzeta;
			dN.entry[18][3] = invJ31*dN18drho + invJ32*dN18deta + invJ33*dN18dzeta;
			dN.entry[19][1] = invJ11*dN19drho + invJ12*dN19deta + invJ13*dN19dzeta;
			dN.entry[19][2] = invJ21*dN19drho + invJ22*dN19deta + invJ23*dN19dzeta;
			dN.entry[19][3] = invJ31*dN19drho + invJ32*dN19deta + invJ33*dN19dzeta;
			dN.entry[20][1] = invJ11*dN20drho + invJ12*dN20deta + invJ13*dN20dzeta;
			dN.entry[20][2] = invJ21*dN20drho + invJ22*dN20deta + invJ23*dN20dzeta;
			dN.entry[20][3] = invJ31*dN20drho + invJ32*dN20deta + invJ33*dN20dzeta;
			dN.entry[21][1] = invJ11*dN21drho + invJ12*dN21deta + invJ13*dN21dzeta;
			dN.entry[21][2] = invJ21*dN21drho + invJ22*dN21deta + invJ23*dN21dzeta;
			dN.entry[21][3] = invJ31*dN21drho + invJ32*dN21deta + invJ33*dN21dzeta;
			dN.entry[22][1] = invJ11*dN22drho + invJ12*dN22deta + invJ13*dN22dzeta;
			dN.entry[22][2] = invJ21*dN22drho + invJ22*dN22deta + invJ23*dN22dzeta;
			dN.entry[22][3] = invJ31*dN22drho + invJ32*dN22deta + invJ33*dN22dzeta;
			dN.entry[23][1] = invJ11*dN23drho + invJ12*dN23deta + invJ13*dN23dzeta;
			dN.entry[23][2] = invJ21*dN23drho + invJ22*dN23deta + invJ23*dN23dzeta;
			dN.entry[23][3] = invJ31*dN23drho + invJ32*dN23deta + invJ33*dN23dzeta;
			dN.entry[24][1] = invJ11*dN24drho + invJ12*dN24deta + invJ13*dN24dzeta;
			dN.entry[24][2] = invJ21*dN24drho + invJ22*dN24deta + invJ23*dN24dzeta;
			dN.entry[24][3] = invJ31*dN24drho + invJ32*dN24deta + invJ33*dN24dzeta;
			dN.entry[25][1] = invJ11*dN25drho + invJ12*dN25deta + invJ13*dN25dzeta;
			dN.entry[25][2] = invJ21*dN25drho + invJ22*dN25deta + invJ23*dN25dzeta;
			dN.entry[25][3] = invJ31*dN25drho + invJ32*dN25deta + invJ33*dN25dzeta;
			dN.entry[26][1] = invJ11*dN26drho + invJ12*dN26deta + invJ13*dN26dzeta;
			dN.entry[26][2] = invJ21*dN26drho + invJ22*dN26deta + invJ23*dN26dzeta;
			dN.entry[26][3] = invJ31*dN26drho + invJ32*dN26deta + invJ33*dN26dzeta;
			dN.entry[27][1] = invJ11*dN27drho + invJ12*dN27deta + invJ13*dN27dzeta;
			dN.entry[27][2] = invJ21*dN27drho + invJ22*dN27deta + invJ23*dN27dzeta;
			dN.entry[27][3] = invJ31*dN27drho + invJ32*dN27deta + invJ33*dN27dzeta;
		}


		bool TriangleGlobal(int element_id, const Vector<T>& point_global, T& rho, T& eta) const throw()
		{
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];

			// Affine transformation for triangle elements.wxm
			register T x = point_global.entry[1];
			register T y = point_global.entry[2];

			T A = x2 - x1;
			T B = x1 - x3;
			T C = 1.0/(A*y3 + B*y2 + (x3 - x2)*y1);
			T D = x1 - x;

			rho = (B*y - D*y3 - (x - x3)*y1)*C;
			eta = (D*y2 + (x - x2)*y1 + A*y)*C;

			if ((rho < 0) || (rho > 1) || (eta < 0) || (eta > 1) || (rho + eta > 1))
			{
				return false;
			}
			return true;
		}


		bool Triangle3Global(int element_id, const Vector<T>& point_global, Vector<T>& N) const throw()
		{
			T rho;
			T eta;
			if (TriangleGlobal(element_id, point_global, rho, eta))
			{
				N.entry[1] = 1 - rho - eta;
				N.entry[2] = rho;
				N.entry[3] = eta;
				return true;
			}
			return false;
		}


		bool Triangle6Global(int element_id, const Vector<T>& point_global, Vector<T>& N) const throw()
		{
			T rho;
			T eta;
			if (TriangleGlobal(element_id, point_global, rho, eta))
			{
				T xi = 1 - rho - eta;

				N.entry[1] = (2*xi - 1)*xi;;
				N.entry[2] = rho*(2*rho - 1);
				N.entry[3] = eta*(2*eta - 1);
				N.entry[4] = 4*rho*xi;
				N.entry[5] = 4*rho*eta;
				N.entry[6] = 4*eta*xi;
				return true;
			}
			return false;
		}


		bool QuadrilateralGlobal(int element_id, const Vector<T>& point_global, T& rho, T& eta) const throw()
		{
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];

			T x = point_global.entry[1];
			T y = point_global.entry[2];

			T A = (x1 + x2 + x3 + x4)*0.25;
			T B = (-x1 + x2 + x3 - x4)*0.25;
			T C = (-x1 - x2 + x3 + x4)*0.25;
			T D = (x1 - x2 + x3 - x4)*0.25;

			T E = (y1 + y2 + y3 + y4)*0.25;
			T F = (-y1 + y2 + y3 - y4)*0.25;
			T G = (-y1 - y2 + y3 + y4)*0.25;
			T H = (y1 - y2 + y3 - y4)*0.25;

			// Newton's method (http://en.wikipedia.org/wiki/Newton%27s_method#Nonlinear_systems_of_equations)
			register T r_rho = 0;
			register T r_eta = 0;
			for(int i = 0; i < 3; ++i)
			{
				T f1 = A + r_rho*B + r_eta*C + r_rho*r_eta*D - x;
				T f2 = E + r_rho*F + r_eta*G + r_rho*r_eta*H - y;

				if ((fabs(f1) < 1e-5) && (fabs(f2) < 1e-5))
				{
					if ((r_rho >= -1) && (r_rho <= 1) && (r_eta >= -1) && (r_eta <= 1))
					{
						rho = r_rho;
						eta = r_eta;
						return true;
					}
					return false;
				}

				T df1dp = B + r_eta*D;
				T df1dn = C + r_rho*D;
				T df2dp = F + r_eta*H;
				T df2dn = G + r_rho*H;

				T det = df1dp*df2dn - df1dn*df2dp;
				r_rho -= (df2dn*f1 - df1dn*f2)/det;
				r_eta -= (df1dp*f2 - df2dp*f1)/det;
			}

			return false;
		}


		bool Quadrilateral4Global(int element_id, const Vector<T>& point_global, Vector<T>& N) const throw()
		{
			T rho;
			T eta;
			if (QuadrilateralGlobal(element_id, point_global, rho, eta))
			{
				N.entry[1] = T(0.25)*(1 - rho)*(1 - eta);
				N.entry[2] = T(0.25)*(1 + rho)*(1 - eta);
				N.entry[3] = T(0.25)*(1 + rho)*(1 + eta);
				N.entry[4] = T(0.25)*(1 - rho)*(1 + eta);
				return true;
			}
			return false;
		}


		bool Quadrilateral8Global(int element_id, const Vector<T>& point_global, Vector<T>& N) const throw()
		{
			T rho;
			T eta;
			if (QuadrilateralGlobal(element_id, point_global, rho, eta))
			{
				N.entry[1] = T(0.25)*(1 - rho)*(1 - eta)*(-1 - rho - eta);
				N.entry[2] = T(0.25)*(1 + rho)*(1 - eta)*(-1 + rho - eta);
				N.entry[3] = T(0.25)*(1 + rho)*(1 + eta)*(-1 + rho + eta);
				N.entry[4] = T(0.25)*(1 - rho)*(1 + eta)*(-1 - rho + eta);
				N.entry[5] = T(0.5)*(1 - rho*rho)*(1 - eta);
				N.entry[6] = T(0.5)*(1 - eta*eta)*(1 + rho);
				N.entry[7] = T(0.5)*(1 - rho*rho)*(1 + eta);
				N.entry[8] = T(0.5)*(1 - eta*eta)*(1 - rho);
				return true;
			}
			return false;
		}


		bool Quadrilateral9Global(int element_id, const Vector<T>& point_global, Vector<T>& N) const throw()
		{
			T rho;
			T eta;
			if (QuadrilateralGlobal(element_id, point_global, rho, eta))
			{
				N.entry[1] = T( 0.25)*(1 - rho)*rho*(1 - eta)*eta;
				N.entry[2] = T(-0.25)*(1 + rho)*rho*(1 - eta)*eta;
				N.entry[3] = T( 0.25)*(1 + rho)*rho*(1 + eta)*eta;
				N.entry[4] = T(-0.25)*(1 - rho)*rho*(1 + eta)*eta;
				N.entry[5] = T(-0.5)*(1 - rho*rho)*(1 - eta)*eta;
				N.entry[6] = T( 0.5)*(1 - eta*eta)*(1 + rho)*rho;
				N.entry[7] = T( 0.5)*(1 - rho*rho)*(1 + eta)*eta;
				N.entry[8] = T(-0.5)*(1 - eta*eta)*(1 - rho)*rho;
				N.entry[9] = (1 - rho*rho)*(1 - eta*eta);
				return true;
			}
			return false;
		}


		bool TetrahedronGlobal(int element_id, const Vector<T>& point_global, T& rho, T& eta, T& zeta) const throw()
		{
			int node1 = mesh.connectivity.entry[element_id][1];
			int node2 = mesh.connectivity.entry[element_id][2];
			int node3 = mesh.connectivity.entry[element_id][3];
			int node4 = mesh.connectivity.entry[element_id][4];
			T x1 = nodes.coordinate.entry[node1][1];
			T y1 = nodes.coordinate.entry[node1][2];
			T z1 = nodes.coordinate.entry[node1][3];
			T x2 = nodes.coordinate.entry[node2][1];
			T y2 = nodes.coordinate.entry[node2][2];
			T z2 = nodes.coordinate.entry[node2][3];
			T x3 = nodes.coordinate.entry[node3][1];
			T y3 = nodes.coordinate.entry[node3][2];
			T z3 = nodes.coordinate.entry[node3][3];
			T x4 = nodes.coordinate.entry[node4][1];
			T y4 = nodes.coordinate.entry[node4][2];
			T z4 = nodes.coordinate.entry[node4][3];

			register T x = point_global.entry[1];
			register T y = point_global.entry[2];
			register T z = point_global.entry[3];

			T factor = 1.0/(((x2 - x1)*y3 + (x1 - x3)*y2 + (x3 - x2)*y1)*z4 + ((x1 - x2)*y4 + (x4 - x1)*y2 + (x2 - x4)*y1)*z3 + ((x3 - x1)*y4 + (x1 - x4)*y3 + (x4 - x3)*y1)*z2 + ((x2 - x3)*y4 + (x4 - x2)*y3 + (x3 - x4)*y2)*z1);
			rho = (x*((y3 - y1)*z4 + (y1 - y4)*z3 + (y4 - y3)*z1) + y*((x1 - x3)*z4 + (x4 - x1)*z3 + (x3 - x4)*z1) + (x3*y1 - x1*y3)*z4 + (x1*y4 - x4*y1)*z3 + (x4*y3 - x3*y4)*z1 + ((x3 - x1)*y4 + (x1 - x4)*y3 + (x4 - x3)*y1)*z)*factor;
			eta = -(x*((y2 - y1)*z4 + (y1 - y4)*z2 + (y4 - y2)*z1) + y*((x1 - x2)*z4 + (x4 - x1)*z2 + (x2 - x4)*z1) + (x2*y1 - x1*y2)*z4 + (x1*y4 - x4*y1)*z2 + (x4*y2 - x2*y4)*z1 + ((x2 - x1)*y4 + (x1 - x4)*y2 + (x4 - x2)*y1)*z)*factor;
			zeta = (x*((y2 - y1)*z3 + (y1 - y3)*z2 + (y3 - y2)*z1) + y*((x1 - x2)*z3 + (x3 - x1)*z2 + (x2 - x3)*z1) + (x2*y1 - x1*y2)*z3 + (x1*y3 - x3*y1)*z2 + (x3*y2 - x2*y3)*z1 + ((x2 - x1)*y3 + (x1 - x3)*y2 + (x3 - x2)*y1)*z)*factor;

			if ((rho < 0) || (rho > 1) || (eta < 0) || (eta > 1) || (zeta < 0) || (zeta > 1) || (rho + eta + zeta > 1))
			{
				return false;
			}
			return true;
		}


		bool Tetrahedron4Global(int element_id, const Vector<T>& point_global, Vector<T>& N) const throw()
		{
			T rho;
			T eta;
			T zeta;
			if (TetrahedronGlobal(element_id, point_global, rho, eta, zeta))
			{
				N.entry[1] = 1 - rho - eta - zeta;
				N.entry[2] = rho;
				N.entry[3] = eta;
				N.entry[4] = zeta;
				return true;
			}
			return false;
		}


		bool Tetrahedron10Global(int element_id, const Vector<T>& point_global, Vector<T>& N) const throw()
		{
			T rho;
			T eta;
			T zeta;
			if (TetrahedronGlobal(element_id, point_global, rho, eta, zeta))
			{
				N.entry[1]  = (1 - rho - eta - zeta)*(1 - 2*rho - 2*eta - 2*zeta);
				N.entry[2]  = rho*(2*rho - 1);
				N.entry[3]  = eta*(2*eta - 1);
				N.entry[4]  = zeta*(2*zeta - 1);
				N.entry[5]  = 4*rho*(1 - rho - eta - zeta);
				N.entry[6]  = 4*rho*eta;
				N.entry[7]  = 4*eta*(1 - rho - eta - zeta);
				N.entry[8]  = 4*zeta*(1 - rho - eta - zeta);
				N.entry[9]  = 4*rho*zeta;
				N.entry[10] = 4*eta*zeta;
				return true;
			}
			return false;
		}


		bool HexahedronGlobal(int element_id, const Vector<T>& point_global, T& rho, T& eta, T& zeta) const throw()
		{
			int node1  = mesh.connectivity.entry[element_id][1];
			int node2  = mesh.connectivity.entry[element_id][2];
			int node3  = mesh.connectivity.entry[element_id][3];
			int node4  = mesh.connectivity.entry[element_id][4];
			int node5  = mesh.connectivity.entry[element_id][5];
			int node6  = mesh.connectivity.entry[element_id][6];
			int node7  = mesh.connectivity.entry[element_id][7];
			int node8  = mesh.connectivity.entry[element_id][8];
			T x1  = nodes.coordinate.entry[node1][1];
			T y1  = nodes.coordinate.entry[node1][2];
			T z1  = nodes.coordinate.entry[node1][3];
			T x2  = nodes.coordinate.entry[node2][1];
			T y2  = nodes.coordinate.entry[node2][2];
			T z2  = nodes.coordinate.entry[node2][3];
			T x3  = nodes.coordinate.entry[node3][1];
			T y3  = nodes.coordinate.entry[node3][2];
			T z3  = nodes.coordinate.entry[node3][3];
			T x4  = nodes.coordinate.entry[node4][1];
			T y4  = nodes.coordinate.entry[node4][2];
			T z4  = nodes.coordinate.entry[node4][3];
			T x5  = nodes.coordinate.entry[node5][1];
			T y5  = nodes.coordinate.entry[node5][2];
			T z5  = nodes.coordinate.entry[node5][3];
			T x6  = nodes.coordinate.entry[node6][1];
			T y6  = nodes.coordinate.entry[node6][2];
			T z6  = nodes.coordinate.entry[node6][3];
			T x7  = nodes.coordinate.entry[node7][1];
			T y7  = nodes.coordinate.entry[node7][2];
			T z7  = nodes.coordinate.entry[node7][3];
			T x8  = nodes.coordinate.entry[node8][1];
			T y8  = nodes.coordinate.entry[node8][2];
			T z8  = nodes.coordinate.entry[node8][3];

			T x = point_global.entry[1];
			T y = point_global.entry[2];
			T z = point_global.entry[3];

			T A = (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)*0.125;
			T B = (-x1 + x2 + x3 - x4 - x5 + x6 + x7 - x8)*0.125;
			T C = (-x1 - x2 + x3 + x4 - x5 - x6 + x7 + x8)*0.125;
			T D = (-x1 - x2 - x3 - x4 + x5 + x6 + x7 + x8)*0.125;
			T E = (x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8)*0.125;
			T F = (x1 - x2 - x3 + x4 - x5 + x6 + x7 - x8)*0.125;
			T G = (x1 + x2 - x3 - x4 - x5 - x6 + x7 + x8)*0.125;
			T H = (-x1 + x2 - x3 + x4 + x5 - x6 + x7 - x8)*0.125;

			T I = (y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8)*0.125;
			T J = (-y1 + y2 + y3 - y4 - y5 + y6 + y7 - y8)*0.125;
			T K = (-y1 - y2 + y3 + y4 - y5 - y6 + y7 + y8)*0.125;
			T L = (-y1 - y2 - y3 - y4 + y5 + y6 + y7 + y8)*0.125;
			T M = (y1 - y2 + y3 - y4 + y5 - y6 + y7 - y8)*0.125;
			T N = (y1 - y2 - y3 + y4 - y5 + y6 + y7 - y8)*0.125;
			T O = (y1 + y2 - y3 - y4 - y5 - y6 + y7 + y8)*0.125;
			T P = (-y1 + y2 - y3 + y4 + y5 - y6 + y7 - y8)*0.125;

			T Q = (z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8)*0.125;
			T R = (-z1 + z2 + z3 - z4 - z5 + z6 + z7 - z8)*0.125;
			T S = (-z1 - z2 + z3 + z4 - z5 - z6 + z7 + z8)*0.125;
			T U = (-z1 - z2 - z3 - z4 + z5 + z6 + z7 + z8)*0.125;
			T V = (z1 - z2 + z3 - z4 + z5 - z6 + z7 - z8)*0.125;
			T W = (z1 - z2 - z3 + z4 - z5 + z6 + z7 - z8)*0.125;
			T X = (z1 + z2 - z3 - z4 - z5 - z6 + z7 + z8)*0.125;
			T Y = (-z1 + z2 - z3 + z4 + z5 - z6 + z7 - z8)*0.125;

			// Newton's method (http://en.wikipedia.org/wiki/Newton%27s_method#Nonlinear_systems_of_equations)
			register T r_rho = 0;
			register T r_eta = 0;
			register T r_zeta = 0;
			for(int i = 0; i < 5; ++i)
			{
				T f1 = A + B*r_rho + C*r_eta + D*r_zeta + E*r_rho*r_eta + F*r_rho*r_zeta + G*r_eta*r_zeta + H*r_rho*r_eta*r_zeta - x;
				T f2 = I + J*r_rho + K*r_eta + L*r_zeta + M*r_rho*r_eta + N*r_rho*r_zeta + O*r_eta*r_zeta + P*r_rho*r_eta*r_zeta - y;
				T f3 = Q + R*r_rho + S*r_eta + U*r_zeta + V*r_rho*r_eta + W*r_rho*r_zeta + X*r_eta*r_zeta + Y*r_rho*r_eta*r_zeta - z;

				if ((fabs(f1) < 1e-5) && (fabs(f2) < 1e-5) && (fabs(f3) < 1e-5))
				{
/**/// Relax this condition!! 1+eps
					if ((r_rho >= -1) && (r_rho <= 1) && (r_eta >= -1) && (r_eta <= 1) && (r_zeta >= -1) && (r_zeta <= 1))
					{
						rho = r_rho;
						eta = r_eta;
						zeta = r_zeta;
						return true;
					}
					return false;
				}

				T df1dp = B + E*r_eta + F*r_zeta + H*r_eta*r_zeta;
				T df1dn = C + E*r_rho + G*r_zeta + H*r_rho*r_zeta;
				T df1dq = D + F*r_rho + G*r_eta + H*r_rho*r_eta;
				T df2dp = J + M*r_eta + N*r_zeta + P*r_eta*r_zeta;
				T df2dn = K + M*r_rho + O*r_zeta + P*r_rho*r_zeta;
				T df2dq = L + N*r_rho + O*r_eta + P*r_rho*r_eta;
				T df3dp = R + V*r_eta + W*r_zeta + Y*r_eta*r_zeta;
				T df3dn = S + V*r_rho + X*r_zeta + Y*r_rho*r_zeta;
				T df3dq = U + W*r_rho + X*r_eta + Y*r_rho*r_eta;

				T det = df1dp*(df2dn*df3dq - df2dq*df3dn) - df1dn*(df3dq*df2dp - df2dq*df3dp) + df1dq*(df2dp*df3dn - df2dn*df3dp);

				r_rho -= ((df2dn*df3dq - df2dq*df3dn)*f1 + (df1dq*df3dn - df1dn*df3dq)*f2 + (df1dn*df2dq - df1dq*df2dn)*f3)/det;
				r_eta -= ((df2dq*df3dp - df2dp*df3dq)*f1 + (df1dp*df3dq - df1dq*df3dp)*f2 + (df1dq*df2dp - df1dp*df2dq)*f3)/det;
				r_zeta -= ((df2dp*df3dn - df2dn*df3dp)*f1 + (df3dp*df1dn - df1dp*df3dn)*f2 + (df1dp*df2dn - df1dn*df2dp)*f3)/det;
			}

			return false;
		}


		bool Hexahedron8Global(int element_id, const Vector<T>& point_global, Vector<T>& N) const throw()
		{
			T rho;
			T eta;
			T zeta;
			if (HexahedronGlobal(element_id, point_global, rho, eta, zeta))
			{
				N.entry[1] = T(0.125)*(1 - rho)*(1 - eta)*(1 - zeta);
				N.entry[2] = T(0.125)*(1 + rho)*(1 - eta)*(1 - zeta);
				N.entry[3] = T(0.125)*(1 + rho)*(1 + eta)*(1 - zeta);
				N.entry[4] = T(0.125)*(1 - rho)*(1 + eta)*(1 - zeta);
				N.entry[5] = T(0.125)*(1 - rho)*(1 - eta)*(1 + zeta);
				N.entry[6] = T(0.125)*(1 + rho)*(1 - eta)*(1 + zeta);
				N.entry[7] = T(0.125)*(1 + rho)*(1 + eta)*(1 + zeta);
				N.entry[8] = T(0.125)*(1 - rho)*(1 + eta)*(1 + zeta);
				return true;
			}
			return false;
		}


		bool Hexahedron20Global(int element_id, const Vector<T>& point_global, Vector<T>& N) const throw()
		{
			T rho;
			T eta;
			T zeta;
			if (HexahedronGlobal(element_id, point_global, rho, eta, zeta))
			{
				N.entry[1] = T(0.125)*(1 - rho)*(1 - eta)*(1 - zeta)*(-rho - eta - zeta - 2);
				N.entry[2] = T(0.125)*(1 + rho)*(1 - eta)*(1 - zeta)*(rho - eta - zeta - 2);
				N.entry[3] = T(0.125)*(1 + rho)*(1 + eta)*(1 - zeta)*(rho + eta - zeta - 2);
				N.entry[4] = T(0.125)*(1 - rho)*(1 + eta)*(1 - zeta)*(-rho + eta - zeta - 2);
				N.entry[5] = T(0.125)*(1 - rho)*(1 - eta)*(1 + zeta)*(-rho - eta + zeta - 2);
				N.entry[6] = T(0.125)*(1 + rho)*(1 - eta)*(1 + zeta)*(rho - eta + zeta - 2);
				N.entry[7] = T(0.125)*(1 + rho)*(1 + eta)*(1 + zeta)*(rho + eta + zeta - 2);
				N.entry[8] = T(0.125)*(1 - rho)*(1 + eta)*(1 + zeta)*(-rho + eta + zeta - 2);
				N.entry[9] = T(0.25)*(1 - rho*rho)*(1 - eta)*(1 - zeta);
				N.entry[10] = T(0.25)*(1 + rho)*(1 - eta*eta)*(1 - zeta);
				N.entry[11] = T(0.25)*(1 - rho*rho)*(1 + eta)*(1 - zeta);
				N.entry[12] = T(0.25)*(1 - rho)*(1 - eta*eta)*(1 - zeta);
				N.entry[13] = T(0.25)*(1 - rho)*(1 - eta)*(1 - zeta*zeta);
				N.entry[14] = T(0.25)*(1 + rho)*(1 - eta)*(1 - zeta*zeta);
				N.entry[15] = T(0.25)*(1 + rho)*(1 + eta)*(1 - zeta*zeta);
				N.entry[16] = T(0.25)*(1 - rho)*(1 + eta)*(1 - zeta*zeta);
				N.entry[17] = T(0.25)*(1 - rho*rho)*(1 - eta)*(1 + zeta);
				N.entry[18] = T(0.25)*(1 + rho)*(1 - eta*eta)*(1 + zeta);
				N.entry[19] = T(0.25)*(1 - rho*rho)*(1 + eta)*(1 + zeta);
				N.entry[20] = T(0.25)*(1 - rho)*(1 - eta*eta)*(1 + zeta);
				return true;
			}
			return false;
		}


		bool Hexahedron27Global(int element_id, const Vector<T>& point_global, Vector<T>& N) const throw()
		{
			T rho;
			T eta;
			T zeta;
			if (HexahedronGlobal(element_id, point_global, rho, eta, zeta))
			{
				N.entry[1]  =  ((eta - 1)*eta*(rho - 1)*rho*(zeta - 1)*zeta)*T(0.125);
				N.entry[2]  =  ((eta - 1)*eta*rho*(rho + 1)*(zeta - 1)*zeta)*T(0.125);
				N.entry[3]  =  (eta*(eta + 1)*rho*(rho + 1)*(zeta - 1)*zeta)*T(0.125);
				N.entry[4]  =  (eta*(eta + 1)*(rho - 1)*rho*(zeta - 1)*zeta)*T(0.125);
				N.entry[5]  =  ((eta - 1)*eta*(rho - 1)*rho*zeta*(zeta + 1))*T(0.125);
				N.entry[6]  =  ((eta - 1)*eta*rho*(rho + 1)*zeta*(zeta + 1))*T(0.125);
				N.entry[7]  =  (eta*(eta + 1)*rho*(rho + 1)*zeta*(zeta + 1))*T(0.125);
				N.entry[8]  =  (eta*(eta + 1)*(rho - 1)*rho*zeta*(zeta + 1))*T(0.125);
				N.entry[9]  = -((eta - 1)*eta*(rho - 1)*(rho + 1)*(zeta - 1)*zeta)*T(0.25);
				N.entry[10] = -((eta - 1)*(eta + 1)*rho*(rho + 1)*(zeta - 1)*zeta)*T(0.25);
				N.entry[11] = -(eta*(eta + 1)*(rho - 1)*(rho + 1)*(zeta - 1)*zeta)*T(0.25);
				N.entry[12] = -((eta - 1)*(eta + 1)*(rho - 1)*rho*(zeta - 1)*zeta)*T(0.25);
				N.entry[13] = -((eta - 1)*eta*(rho - 1)*rho*(zeta - 1)*(zeta + 1))*T(0.25);
				N.entry[14] = -((eta - 1)*eta*rho*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
				N.entry[15] = -(eta*(eta + 1)*rho*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.25);
				N.entry[16] = -(eta*(eta + 1)*(rho - 1)*rho*(zeta - 1)*(zeta + 1))*T(0.25);
				N.entry[17] = -((eta - 1)*eta*(rho - 1)*(rho + 1)*zeta*(zeta + 1))*T(0.25);
				N.entry[18] = -((eta - 1)*(eta + 1)*rho*(rho + 1)*zeta*(zeta + 1))*T(0.25);
				N.entry[19] = -(eta*(eta + 1)*(rho - 1)*(rho + 1)*zeta*(zeta + 1))*T(0.25);
				N.entry[20] = -((eta - 1)*(eta + 1)*(rho - 1)*rho*zeta*(zeta + 1))*T(0.25);
				N.entry[21] =  ((eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*(zeta - 1)*zeta)*T(0.5);
				N.entry[22] =  ((eta - 1)*eta*(rho - 1)*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.5);
				N.entry[23] =  ((eta - 1)*(eta + 1)*rho*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.5);
				N.entry[24] =  (eta*(eta + 1)*(rho - 1)*(rho + 1)*(zeta - 1)*(zeta + 1))*T(0.5);
				N.entry[25] =  ((eta - 1)*(eta + 1)*(rho - 1)*rho*(zeta - 1)*(zeta + 1))*T(0.5);
				N.entry[26] =  ((eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*zeta*(zeta + 1))*T(0.5);
				N.entry[27] = -(eta - 1)*(eta + 1)*(rho - 1)*(rho + 1)*(zeta - 1)*(zeta + 1);
				return true;
			}
			return false;
		}


		void FacetLinear2(const Vector<int>& facet_nodes, const Vector<T>& point, Vector<T>& N, T& det_J) const throw()
		{
			T rho = point.entry[1];

			// Shape functions
			N.entry[1] = (1 - rho)*0.5;
			N.entry[2] = (rho + 1)*0.5;

			// Global coordinates relative to (gx1, gy1)
			T gx1 = nodes.coordinate.entry[facet_nodes.entry[1]][1];
			T gy1 = nodes.coordinate.entry[facet_nodes.entry[1]][2];
			T gx2 = nodes.coordinate.entry[facet_nodes.entry[2]][1] - gx1;
			T gy2 = nodes.coordinate.entry[facet_nodes.entry[2]][2] - gy1;

			// Local coordinates
			T x1 = 0;
			T x2 = sqrt(gx2*gx2 + gy2*gy2);

			// Jacobian
			T dxdrho = (x2 - x1)*0.5;

			// Jacobian determinant
			det_J = dxdrho;
		}


		void FacetLinear3(const Vector<int>& facet_nodes, const Vector<T>& point, Vector<T>& N, T& det_J) const throw()
		{
			T rho = point.entry[1];

			// Shape functions
			N.entry[1] = rho*(rho - 1)*0.5;
			N.entry[2] = rho*(rho + 1)*0.5;
			N.entry[3] = (1 - rho)*(rho + 1);

			// Global coordinates relative to (gx1, gy1)
			T gx1 = nodes.coordinate.entry[facet_nodes.entry[1]][1];
			T gy1 = nodes.coordinate.entry[facet_nodes.entry[1]][2];
			T gx2 = nodes.coordinate.entry[facet_nodes.entry[2]][1] - gx1;
			T gy2 = nodes.coordinate.entry[facet_nodes.entry[2]][2] - gy1;
			T gx3 = nodes.coordinate.entry[facet_nodes.entry[3]][1] - gx1;
			T gy3 = nodes.coordinate.entry[facet_nodes.entry[3]][2] - gy1;
			T size_i = sqrt(gx2*gx2 + gy2*gy2);
			T ix = gx2/size_i;
			T iy = gy2/size_i;

			// Local coordinates
			T x1 = 0;
			T x2 = sqrt(gx2*gx2 + gy2*gy2);
			T x3 = gx3*ix + gy3*iy;

			// Jacobian
			T dxdrho = -(4*rho*x3 + (-2*rho - 1)*x2 + (1 - 2*rho)*x1)*0.5;

			// Jacobian determinant
			det_J = dxdrho;
		}


		void FacetTriangle3(const Vector<int>& facet_nodes, const Vector<T>& point, Vector<T>& N, T& det_J) const throw()
		{
			T rho = point.entry[1];
			T eta = point.entry[2];

			// Shape functions
			N.entry[1] = 1 - rho - eta;
			N.entry[2] = rho;
			N.entry[3] = eta;

			// Global coordinates relative to (ox, oy, oz)
			T gx1 = nodes.coordinate.entry[facet_nodes.entry[1]][1];
			T gy1 = nodes.coordinate.entry[facet_nodes.entry[1]][2];
			T gz1 = nodes.coordinate.entry[facet_nodes.entry[1]][3];
			T gx2 = nodes.coordinate.entry[facet_nodes.entry[2]][1] - gx1;
			T gy2 = nodes.coordinate.entry[facet_nodes.entry[2]][2] - gy1;
			T gz2 = nodes.coordinate.entry[facet_nodes.entry[2]][3] - gz1;
			T gx3 = nodes.coordinate.entry[facet_nodes.entry[3]][1] - gx1;
			T gy3 = nodes.coordinate.entry[facet_nodes.entry[3]][2] - gy1;
			T gz3 = nodes.coordinate.entry[facet_nodes.entry[3]][3] - gz1;
			T size_i = sqrt(gx2*gx2 + gy2*gy2 + gz2*gz2);
			T ix = gx2/size_i;
			T iy = gy2/size_i;
			T iz = gz2/size_i;
			T jx = gz2*(gx3*gz2 - gx2*gz3) - gy2*(gx2*gy3 - gx3*gy2);
			T jy = gx2*(gx2*gy3 - gx3*gy2) - gz2*(gy2*gz3 - gy3*gz2);
			T jz = gy2*(gy2*gz3 - gy3*gz2) - gx2*(gx3*gz2 - gx2*gz3);
			T size_j = sqrt(jx*jx + jy*jy + jz*jz);
			jx /= size_j;
			jy /= size_j;
			jz /= size_j;

			// Local coordinates 
			T x1 = 0;
			T y1 = 0;
			T x2 = sqrt(gx2*gx2 + gy2*gy2 + gz2*gz2);
			T y2 = 0;
			T x3 = gx3*ix + gy3*iy + gz3*iz;
			T y3 = gx3*jx + gy3*jy + gz3*jz;

			// Jacobian
			T dxdrho = x2 - x1;
			T dxdeta = x3 - x1;
			T dydrho = y2 - y1;
			T dydeta = y3 - y1;

			// Jacobian determinant
			det_J = dxdrho*dydeta - dxdeta*dydrho;
		}


		void FacetTriangle6(const Vector<int>& facet_nodes, const Vector<T>& point, Vector<T>& N, T& det_J) const throw()
		{
			T rho = point.entry[1];
			T eta = point.entry[2];
			T xi = 1 - rho - eta;

			// Shape functions
			N.entry[1] = (2*xi - 1)*xi;;
			N.entry[2] = rho*(2*rho - 1);
			N.entry[3] = eta*(2*eta - 1);
			N.entry[4] = 4*rho*xi;
			N.entry[5] = 4*rho*eta;
			N.entry[6] = 4*eta*xi;

			// Global coordinates relative to (gx1, gy1, gz1)
			T gx1 = nodes.coordinate.entry[facet_nodes.entry[1]][1];
			T gy1 = nodes.coordinate.entry[facet_nodes.entry[1]][2];
			T gz1 = nodes.coordinate.entry[facet_nodes.entry[1]][3];
			T gx2 = nodes.coordinate.entry[facet_nodes.entry[2]][1] - gx1;
			T gy2 = nodes.coordinate.entry[facet_nodes.entry[2]][2] - gy1;
			T gz2 = nodes.coordinate.entry[facet_nodes.entry[2]][3] - gz1;
			T gx3 = nodes.coordinate.entry[facet_nodes.entry[3]][1] - gx1;
			T gy3 = nodes.coordinate.entry[facet_nodes.entry[3]][2] - gy1;
			T gz3 = nodes.coordinate.entry[facet_nodes.entry[3]][3] - gz1;
			T gx4 = nodes.coordinate.entry[facet_nodes.entry[4]][1] - gx1;
			T gy4 = nodes.coordinate.entry[facet_nodes.entry[4]][2] - gy1;
			T gz4 = nodes.coordinate.entry[facet_nodes.entry[4]][3] - gz1;
			T gx5 = nodes.coordinate.entry[facet_nodes.entry[5]][1] - gx1;
			T gy5 = nodes.coordinate.entry[facet_nodes.entry[5]][2] - gy1;
			T gz5 = nodes.coordinate.entry[facet_nodes.entry[5]][3] - gz1;
			T gx6 = nodes.coordinate.entry[facet_nodes.entry[6]][1] - gx1;
			T gy6 = nodes.coordinate.entry[facet_nodes.entry[6]][2] - gy1;
			T gz6 = nodes.coordinate.entry[facet_nodes.entry[6]][3] - gz1;
			T size_i = sqrt(gx2*gx2 + gy2*gy2 + gz2*gz2);
			T ix = gx2/size_i;
			T iy = gy2/size_i;
			T iz = gz2/size_i;
			T jx = gz2*(gx3*gz2 - gx2*gz3) - gy2*(gx2*gy3 - gx3*gy2);
			T jy = gx2*(gx2*gy3 - gx3*gy2) - gz2*(gy2*gz3 - gy3*gz2);
			T jz = gy2*(gy2*gz3 - gy3*gz2) - gx2*(gx3*gz2 - gx2*gz3);
			T size_j = sqrt(jx*jx + jy*jy + jz*jz);
			jx /= size_j;
			jy /= size_j;
			jz /= size_j;

			// Local coordinates 
			T x1 = 0;
			T y1 = 0;
			T x2 = sqrt(gx2*gx2 + gy2*gy2 + gz2*gz2);
			T y2 = 0;
			T x3 = gx3*ix + gy3*iy + gz3*iz;
			T y3 = gx3*jx + gy3*jy + gz3*jz;
			T x4 = gx4*ix + gy4*iy + gz4*iz;
			T y4 = gx4*jx + gy4*jy + gz4*jz;
			T x5 = gx5*ix + gy5*iy + gz5*iz;
			T y5 = gx5*jx + gy5*jy + gz5*jz;
			T x6 = gx6*ix + gy6*iy + gz6*iz;
			T y6 = gx6*jx + gy6*jy + gz6*jz;

			// Jacobian
			T dxdrho = -4*eta*x6 + 4*eta*x5 + 4*(xi - rho)*x4 + (4*rho - 1)*x2 - (4*xi - 1)*x1;
			T dxdeta =  4*(xi - eta)*x6 + 4*rho*x5 - 4*rho*x4 + (4*eta - 1)*x3 - (4*xi - 1)*x1;
			T dydrho = -4*eta*y6 + 4*eta*y5 + 4*(xi - rho)*y4 + (4*rho - 1)*y2 - (4*xi - 1)*y1;
			T dydeta =  4*(xi - eta)*y6 + 4*rho*y5 - 4*rho*y4 + (4*eta - 1)*y3 - (4*xi - 1)*y1;

			// Jacobian determinant
			det_J = dxdrho*dydeta - dxdeta*dydrho;
		}


		void FacetQuadrilateral4(const Vector<int>& facet_nodes, const Vector<T>& point, Vector<T>& N, T& det_J) const throw()
		{
			T rho = point.entry[1];
			T eta = point.entry[2];

			// Shape functions
			N.entry[1] = T(0.25)*(1 - rho)*(1 - eta);
			N.entry[2] = T(0.25)*(1 + rho)*(1 - eta);
			N.entry[3] = T(0.25)*(1 + rho)*(1 + eta);
			N.entry[4] = T(0.25)*(1 - rho)*(1 + eta);

			// Global coordinates relative to (gx1, gy1, gz1)
			T gx1 = nodes.coordinate.entry[facet_nodes.entry[1]][1];
			T gy1 = nodes.coordinate.entry[facet_nodes.entry[1]][2];
			T gz1 = nodes.coordinate.entry[facet_nodes.entry[1]][3];
			T gx2 = nodes.coordinate.entry[facet_nodes.entry[2]][1] - gx1;
			T gy2 = nodes.coordinate.entry[facet_nodes.entry[2]][2] - gy1;
			T gz2 = nodes.coordinate.entry[facet_nodes.entry[2]][3] - gz1;
			T gx3 = nodes.coordinate.entry[facet_nodes.entry[3]][1] - gx1;
			T gy3 = nodes.coordinate.entry[facet_nodes.entry[3]][2] - gy1;
			T gz3 = nodes.coordinate.entry[facet_nodes.entry[3]][3] - gz1;
			T gx4 = nodes.coordinate.entry[facet_nodes.entry[4]][1] - gx1;
			T gy4 = nodes.coordinate.entry[facet_nodes.entry[4]][2] - gy1;
			T gz4 = nodes.coordinate.entry[facet_nodes.entry[4]][3] - gz1;
			T size_i = sqrt(gx2*gx2 + gy2*gy2 + gz2*gz2);
			T ix = gx2/size_i;
			T iy = gy2/size_i;
			T iz = gz2/size_i;
			T jx = gz2*(gx3*gz2 - gx2*gz3) - gy2*(gx2*gy3 - gx3*gy2);
			T jy = gx2*(gx2*gy3 - gx3*gy2) - gz2*(gy2*gz3 - gy3*gz2);
			T jz = gy2*(gy2*gz3 - gy3*gz2) - gx2*(gx3*gz2 - gx2*gz3);
			T size_j = sqrt(jx*jx + jy*jy + jz*jz);
			jx /= size_j;
			jy /= size_j;
			jz /= size_j;

			// Local coordinates 
			T x1 = 0;
			T y1 = 0;
			T x2 = sqrt(gx2*gx2 + gy2*gy2 + gz2*gz2);
			T y2 = 0;
			T x3 = gx3*ix + gy3*iy + gz3*iz;
			T y3 = gx3*jx + gy3*jy + gz3*jz;
			T x4 = gx4*ix + gy4*iy + gz4*iz;
			T y4 = gx4*jx + gy4*jy + gz4*jz;

			// Jacobian
			T dxdrho = T(-0.25)*((1 + eta)*x4 - (1 + eta)*x3 - (1 - eta)*x2 + (1 - eta)*x1);
			T dxdeta = T( 0.25)*((1 - rho)*x4 + (1 + rho)*x3 - (1 + rho)*x2 - (1 - rho)*x1);
			T dydrho = T(-0.25)*((1 + eta)*y4 - (1 + eta)*y3 - (1 - eta)*y2 + (1 - eta)*y1);
			T dydeta = T( 0.25)*((1 - rho)*y4 + (1 + rho)*y3 - (1 + rho)*y2 - (1 - rho)*y1);

			// Jacobian determinant
			det_J = dxdrho*dydeta - dxdeta*dydrho;
		}


		void FacetQuadrilateral8(const Vector<int>& facet_nodes, const Vector<T>& point, Vector<T>& N, T& det_J) const throw()
		{
			T rho = point.entry[1];
			T eta = point.entry[2];

			// Shape functions
			N.entry[1] = T(0.25)*(1 - rho)*(1 - eta)*(-1 - rho - eta);
			N.entry[2] = T(0.25)*(1 + rho)*(1 - eta)*(-1 + rho - eta);
			N.entry[3] = T(0.25)*(1 + rho)*(1 + eta)*(-1 + rho + eta);
			N.entry[4] = T(0.25)*(1 - rho)*(1 + eta)*(-1 - rho + eta);
			N.entry[5] = T(0.5)*(1 - rho*rho)*(1 - eta);
			N.entry[6] = T(0.5)*(1 - eta*eta)*(1 + rho);
			N.entry[7] = T(0.5)*(1 - rho*rho)*(1 + eta);
			N.entry[8] = T(0.5)*(1 - eta*eta)*(1 - rho);

			// Global coordinates relative to (gx1, gy1, gz1)
			T gx1 = nodes.coordinate.entry[facet_nodes.entry[1]][1];
			T gy1 = nodes.coordinate.entry[facet_nodes.entry[1]][2];
			T gz1 = nodes.coordinate.entry[facet_nodes.entry[1]][3];
			T gx2 = nodes.coordinate.entry[facet_nodes.entry[2]][1] - gx1;
			T gy2 = nodes.coordinate.entry[facet_nodes.entry[2]][2] - gy1;
			T gz2 = nodes.coordinate.entry[facet_nodes.entry[2]][3] - gz1;
			T gx3 = nodes.coordinate.entry[facet_nodes.entry[3]][1] - gx1;
			T gy3 = nodes.coordinate.entry[facet_nodes.entry[3]][2] - gy1;
			T gz3 = nodes.coordinate.entry[facet_nodes.entry[3]][3] - gz1;
			T gx4 = nodes.coordinate.entry[facet_nodes.entry[4]][1] - gx1;
			T gy4 = nodes.coordinate.entry[facet_nodes.entry[4]][2] - gy1;
			T gz4 = nodes.coordinate.entry[facet_nodes.entry[4]][3] - gz1;
			T gx5 = nodes.coordinate.entry[facet_nodes.entry[5]][1] - gx1;
			T gy5 = nodes.coordinate.entry[facet_nodes.entry[5]][2] - gy1;
			T gz5 = nodes.coordinate.entry[facet_nodes.entry[5]][3] - gz1;
			T gx6 = nodes.coordinate.entry[facet_nodes.entry[6]][1] - gx1;
			T gy6 = nodes.coordinate.entry[facet_nodes.entry[6]][2] - gy1;
			T gz6 = nodes.coordinate.entry[facet_nodes.entry[6]][3] - gz1;
			T gx7 = nodes.coordinate.entry[facet_nodes.entry[7]][1] - gx1;
			T gy7 = nodes.coordinate.entry[facet_nodes.entry[7]][2] - gy1;
			T gz7 = nodes.coordinate.entry[facet_nodes.entry[7]][3] - gz1;
			T gx8 = nodes.coordinate.entry[facet_nodes.entry[8]][1] - gx1;
			T gy8 = nodes.coordinate.entry[facet_nodes.entry[8]][2] - gy1;
			T gz8 = nodes.coordinate.entry[facet_nodes.entry[8]][3] - gz1;
			T size_i = sqrt(gx2*gx2 + gy2*gy2 + gz2*gz2);
			T ix = gx2/size_i;
			T iy = gy2/size_i;
			T iz = gz2/size_i;
			T jx = gz2*(gx3*gz2 - gx2*gz3) - gy2*(gx2*gy3 - gx3*gy2);
			T jy = gx2*(gx2*gy3 - gx3*gy2) - gz2*(gy2*gz3 - gy3*gz2);
			T jz = gy2*(gy2*gz3 - gy3*gz2) - gx2*(gx3*gz2 - gx2*gz3);
			T size_j = sqrt(jx*jx + jy*jy + jz*jz);
			jx /= size_j;
			jy /= size_j;
			jz /= size_j;

			// Local coordinates 
			T x1 = 0;
			T y1 = 0;
			T x2 = sqrt(gx2*gx2 + gy2*gy2 + gz2*gz2);
			T y2 = 0;
			T x3 = gx3*ix + gy3*iy + gz3*iz;
			T y3 = gx3*jx + gy3*jy + gz3*jz;
			T x4 = gx4*ix + gy4*iy + gz4*iz;
			T y4 = gx4*jx + gy4*jy + gz4*jz;
			T x5 = gx5*ix + gy5*iy + gz5*iz;
			T y5 = gx5*jx + gy5*jy + gz5*jz;
			T x6 = gx6*ix + gy6*iy + gz6*iz;
			T y6 = gx6*jx + gy6*jy + gz6*jz;
			T x7 = gx7*ix + gy7*iy + gz7*iz;
			T y7 = gx7*jx + gy7*jy + gz7*jz;
			T x8 = gx8*ix + gy8*iy + gz8*iz;
			T y8 = gx8*jx + gy8*jy + gz8*jz;

			// Jacobian
			T dxdrho = T(0.25)*((2*eta*eta - 2)*x8 + (-4*eta - 4)*rho*x7 + (2 - 2*eta*eta)*x6 + (4*eta - 4)*rho*x5 + ((2*eta + 2)*rho - eta*eta - eta)*x4 + ((2*eta + 2)*rho + eta*eta + eta)*x3 + ((2 - 2*eta)*rho + eta*eta - eta)*x2 + ((2 - 2*eta)*rho - eta*eta + eta)*x1);
			T dxdeta = T(0.25)*((4*eta*rho - 4*eta)*x8 + (2 - 2*rho*rho)*x7 + (-4*eta*rho - 4*eta)*x6 + (2*rho*rho - 2)*x5 + (rho*rho + (-2*eta - 1)*rho + 2*eta)*x4 + (rho*rho + (2*eta + 1)*rho + 2*eta)*x3 + (-rho*rho + (2*eta - 1)*rho + 2*eta)*x2 + (-rho*rho + (1 - 2*eta)*rho + 2*eta)*x1);
			T dydrho = T(0.25)*((2*eta*eta - 2)*y8 + (-4*eta - 4)*rho*y7 + (2 - 2*eta*eta)*y6 + (4*eta - 4)*rho*y5 + ((2*eta + 2)*rho - eta*eta - eta)*y4 + ((2*eta + 2)*rho + eta*eta + eta)*y3 + ((2 - 2*eta)*rho + eta*eta - eta)*y2 + ((2 - 2*eta)*rho - eta*eta + eta)*y1);
			T dydeta = T(0.25)*((4*eta*rho - 4*eta)*y8 + (2 - 2*rho*rho)*y7 + (-4*eta*rho - 4*eta)*y6 + (2*rho*rho - 2)*y5 + (rho*rho + (-2*eta - 1)*rho + 2*eta)*y4 + (rho*rho + (2*eta + 1)*rho + 2*eta)*y3 + (-rho*rho + (2*eta - 1)*rho + 2*eta)*y2 + (-rho*rho + (1 - 2*eta)*rho + 2*eta)*y1);

			// Jacobian determinant
			det_J = dxdrho*dydeta - dxdeta*dydrho;
		}


		void FacetQuadrilateral9(const Vector<int>& facet_nodes, const Vector<T>& point, Vector<T>& N, T& det_J) const throw()
		{
			register T rho = point.entry[1];
			register T eta = point.entry[2];

			// Shape functions
			N.entry[1] = T( 0.25)*(1 - rho)*rho*(1 - eta)*eta;
			N.entry[2] = T(-0.25)*(1 + rho)*rho*(1 - eta)*eta;
			N.entry[3] = T( 0.25)*(1 + rho)*rho*(1 + eta)*eta;
			N.entry[4] = T(-0.25)*(1 - rho)*rho*(1 + eta)*eta;
			N.entry[5] = T(-0.5)*(1 - rho*rho)*(1 - eta)*eta;
			N.entry[6] = T( 0.5)*(1 - eta*eta)*(1 + rho)*rho;
			N.entry[7] = T( 0.5)*(1 - rho*rho)*(1 + eta)*eta;
			N.entry[8] = T(-0.5)*(1 - eta*eta)*(1 - rho)*rho;
			N.entry[9] = (1 - rho*rho)*(1 - eta*eta);

			// Global coordinates relative to (gx1, gy1, gz1)
			T gx1 = nodes.coordinate.entry[facet_nodes.entry[1]][1];
			T gy1 = nodes.coordinate.entry[facet_nodes.entry[1]][2];
			T gz1 = nodes.coordinate.entry[facet_nodes.entry[1]][3];
			T gx2 = nodes.coordinate.entry[facet_nodes.entry[2]][1] - gx1;
			T gy2 = nodes.coordinate.entry[facet_nodes.entry[2]][2] - gy1;
			T gz2 = nodes.coordinate.entry[facet_nodes.entry[2]][3] - gz1;
			T gx3 = nodes.coordinate.entry[facet_nodes.entry[3]][1] - gx1;
			T gy3 = nodes.coordinate.entry[facet_nodes.entry[3]][2] - gy1;
			T gz3 = nodes.coordinate.entry[facet_nodes.entry[3]][3] - gz1;
			T gx4 = nodes.coordinate.entry[facet_nodes.entry[4]][1] - gx1;
			T gy4 = nodes.coordinate.entry[facet_nodes.entry[4]][2] - gy1;
			T gz4 = nodes.coordinate.entry[facet_nodes.entry[4]][3] - gz1;
			T gx5 = nodes.coordinate.entry[facet_nodes.entry[5]][1] - gx1;
			T gy5 = nodes.coordinate.entry[facet_nodes.entry[5]][2] - gy1;
			T gz5 = nodes.coordinate.entry[facet_nodes.entry[5]][3] - gz1;
			T gx6 = nodes.coordinate.entry[facet_nodes.entry[6]][1] - gx1;
			T gy6 = nodes.coordinate.entry[facet_nodes.entry[6]][2] - gy1;
			T gz6 = nodes.coordinate.entry[facet_nodes.entry[6]][3] - gz1;
			T gx7 = nodes.coordinate.entry[facet_nodes.entry[7]][1] - gx1;
			T gy7 = nodes.coordinate.entry[facet_nodes.entry[7]][2] - gy1;
			T gz7 = nodes.coordinate.entry[facet_nodes.entry[7]][3] - gz1;
			T gx8 = nodes.coordinate.entry[facet_nodes.entry[8]][1] - gx1;
			T gy8 = nodes.coordinate.entry[facet_nodes.entry[8]][2] - gy1;
			T gz8 = nodes.coordinate.entry[facet_nodes.entry[8]][3] - gz1;
			T gx9 = nodes.coordinate.entry[facet_nodes.entry[9]][1] - gx1;
			T gy9 = nodes.coordinate.entry[facet_nodes.entry[9]][2] - gy1;
			T gz9 = nodes.coordinate.entry[facet_nodes.entry[9]][3] - gz1;
			T size_i = sqrt(gx2*gx2 + gy2*gy2 + gz2*gz2);
			T ix = gx2/size_i;
			T iy = gy2/size_i;
			T iz = gz2/size_i;
			T jx = gz2*(gx3*gz2 - gx2*gz3) - gy2*(gx2*gy3 - gx3*gy2);
			T jy = gx2*(gx2*gy3 - gx3*gy2) - gz2*(gy2*gz3 - gy3*gz2);
			T jz = gy2*(gy2*gz3 - gy3*gz2) - gx2*(gx3*gz2 - gx2*gz3);
			T size_j = sqrt(jx*jx + jy*jy + jz*jz);
			jx /= size_j;
			jy /= size_j;
			jz /= size_j;

			// Local coordinates 
			T x1 = 0;
			T y1 = 0;
			T x2 = sqrt(gx2*gx2 + gy2*gy2 + gz2*gz2);
			T y2 = 0;
			T x3 = gx3*ix + gy3*iy + gz3*iz;
			T y3 = gx3*jx + gy3*jy + gz3*jz;
			T x4 = gx4*ix + gy4*iy + gz4*iz;
			T y4 = gx4*jx + gy4*jy + gz4*jz;
			T x5 = gx5*ix + gy5*iy + gz5*iz;
			T y5 = gx5*jx + gy5*jy + gz5*jz;
			T x6 = gx6*ix + gy6*iy + gz6*iz;
			T y6 = gx6*jx + gy6*jy + gz6*jz;
			T x7 = gx7*ix + gy7*iy + gz7*iz;
			T y7 = gx7*jx + gy7*jy + gz7*jz;
			T x8 = gx8*ix + gy8*iy + gz8*iz;
			T y8 = gx8*jx + gy8*jy + gz8*jz;
			T x9 = gx9*ix + gy9*iy + gz9*iz;
			T y9 = gx9*jx + gy9*jy + gz9*jz;

			// Jacobian
			T dxdrho = T(0.25)*(8*(eta*eta - 1)*rho*x9 + (4*(1 - eta*eta)*rho + 2*eta*eta - 2)*x8 - 4*( eta + eta*eta)*rho*x7 + (4*(1 - eta*eta)*rho - 2*eta*eta + 2)*x6 + 4*(eta - eta*eta)*rho*x5 + (2*(eta*eta + eta)*rho - eta*eta - eta)*x4 + (2*(eta*eta + eta)*rho + eta*eta + eta)*x3 + (2*(eta*eta - eta)*rho + eta*eta - eta)*x2 + (2*(eta*eta - eta)*rho - eta*eta + eta)*x1);
			T dxdeta = T(0.25)*(8*(eta*rho*rho - eta)*x9 + 4*(eta*rho - eta*rho*rho)*x8 + (( - 4*eta - 2)*rho*rho + 4*eta + 2)*x7 - 4*(eta*rho*rho + eta*rho)*x6 + ((2 - 4*eta)*rho*rho + 4*eta - 2)*x5 + ((2*eta + 1)*rho*rho + ( - 2*eta - 1)*rho)*x4 + ((2*eta + 1)*rho*rho + (2*eta + 1)*rho)*x3 + ((2*eta - 1)*rho*rho + (2*eta - 1)*rho)*x2 + ((2*eta - 1)*rho*rho + (1 - 2*eta)*rho)*x1);
			T dydrho = T(0.25)*(8*(eta*eta - 1)*rho*y9 + (4*(1 - eta*eta)*rho + 2*eta*eta - 2)*y8 - 4*( eta + eta*eta)*rho*y7 + (4*(1 - eta*eta)*rho - 2*eta*eta + 2)*y6 + 4*(eta - eta*eta)*rho*y5 + (2*(eta*eta + eta)*rho - eta*eta - eta)*y4 + (2*(eta*eta + eta)*rho + eta*eta + eta)*y3 + (2*(eta*eta - eta)*rho + eta*eta - eta)*y2 + (2*(eta*eta - eta)*rho - eta*eta + eta)*y1);
			T dydeta = T(0.25)*(8*(eta*rho*rho - eta)*y9 + 4*(eta*rho - eta*rho*rho)*y8 + (( - 4*eta - 2)*rho*rho + 4*eta + 2)*y7 - 4*(eta*rho*rho + eta*rho)*y6 + ((2 - 4*eta)*rho*rho + 4*eta - 2)*y5 + ((2*eta + 1)*rho*rho + ( - 2*eta - 1)*rho)*y4 + ((2*eta + 1)*rho*rho + (2*eta + 1)*rho)*y3 + ((2*eta - 1)*rho*rho + (2*eta - 1)*rho)*y2 + ((2*eta - 1)*rho*rho + (1 - 2*eta)*rho)*y1);

			// Jacobian determinant
			det_J = dxdrho*dydeta - dxdeta*dydrho;
		}


	private:

		ShapeFunctions& operator = (const ShapeFunctions&) throw()
		{
			return *this;
		}
};

#endif
