// ShapeIntegrationRule.h
// Copyright (C) 2010 Miguel Vargas (miguel.vargas@gmail.com)
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

#ifndef _ShapeIntegrationRule_h_
#define _ShapeIntegrationRule_h_

#include <Basic/Memory.h>
#include <Container/Vector.h>
#include <FiniteElement/Shape.h>


template <typename T>
struct IntegrationRule
{
	Vector<T> point;
	T weight;
};


template <typename T>
void ShapeIntegrationRule(const ShapeType shape_type, const int points_count, Vector<IntegrationRule<T> >& integration_rule) throw(Memory::Exception)
{
	Assert(shape_type != shape_undefined);

	try
	{
		if ((shape_type == shape_linear) || (shape_type == shape_quadrilateral) || (shape_type == shape_hexahedron))
		{
			Assert((points_count >= 1) && (points_count <= 7));

			// Abscissae and weights of Gauss quadrature of with 1 point
			static const T a1[] =
			{
				T( 0.0e0L)
			};
			static const T w1[] =
			{
				T( 2.0e0L)
			};

			// Abscissae and weights of Gauss quadrature of with 2 points
			static const T a2[] =
			{
				T(-0.57735026918962576450914878050196e0L),
				T( 0.57735026918962576450914878050196e0L)
			};
			static const T w2[] =
			{
				T( 1.0e0L),
				T( 1.0e0L)
			};

			// Abscissae and weights of Gauss quadrature of with 3 points
			static const T a3[] =
			{
				T(-0.77459666924148337703585307995648e0L),
				T( 0.0e0L),
				T( 0.77459666924148337703585307995648e0L)
			};
			static const T w3[] =
			{
				T( 0.55555555555555555555555555555556e0L),
				T( 0.88888888888888888888888888888889e0L),
				T( 0.55555555555555555555555555555556e0L)
			};

			// Abscissae and weights of Gauss quadrature of width 4 points
			static const T a4[] =
			{
				T(-0.86113631159405257522394648889281e0L),
				T(-0.33998104358485626480266575910324e0L),
				T( 0.33998104358485626480266575910324e0L),
				T( 0.86113631159405257522394648889281e0L)
			};
			static const T w4[] =
			{
				T( 0.34785484513745385737306394922200e0L),
				T( 0.65214515486254614262693605077800e0L),
				T( 0.65214515486254614262693605077800e0L),
				T( 0.34785484513745385737306394922200e0L)
			};

			// Abscissae and weights of Gauss quadrature of width 5 points
			static const T a5[] =
			{
				T(-0.90617984593866399279762687829939e0L),
				T(-0.53846931010568309103631442070021e0L),
				T( 0.0e0L),
				T( 0.53846931010568309103631442070021e0L),
				T( 0.90617984593866399279762687829939e0L)
			};
			static const T w5[] =
			{
				T( 0.23692688505618908751426404071991e0L),
				T( 0.47862867049936646804129151483563e0L),
				T( 0.56888888888888888888888888888888e0L),
				T( 0.47862867049936646804129151483563e0L),
				T( 0.23692688505618908751426404071991e0L)
			};

			// Abscissae and weights of Gauss quadrature of width 6 points
			static const T a6[] =
			{
				T(-0.93246951420315202781230155449399e0L),
				T(-0.66120938646626451366139959501990e0L),
				T(-0.23861918608319690863050172168071e0L),
				T( 0.23861918608319690863050172168071e0L),
				T( 0.66120938646626451366139959501990e0L),
				T( 0.93246951420315202781230155449399e0L)
			};
			static const T w6[] =
			{
				T( 0.17132449237917034504029614217273e0L),
				T( 0.36076157304813860756983351383771e0L),
				T( 0.46791393457269104738987034398955e0L),
				T( 0.46791393457269104738987034398955e0L),
				T( 0.36076157304813860756983351383771e0L),
				T( 0.17132449237917034504029614217273e0L)
			};

			// Abscissae and weights of Gauss quadrature of width 7 points
			static const T a7[] =
			{
				T(-0.94910791234275852452618968404785e0L),
				T(-0.74153118559939443986386477328078e0L),
				T(-0.40584515137739716690660641207696e0L),
				T( 0.e0L),
				T( 0.40584515137739716690660641207696e0L),
				T( 0.74153118559939443986386477328078e0L),
				T( 0.94910791234275852452618968404785e0L)
			};
			static const T w7[] =
			{
				T( 0.12948496616886969327061143267908e0L),
				T( 0.27970539148927666790146777142377e0L),
				T( 0.38183005050511894495036977548897e0L),
				T( 0.41795918367346938775510204081632e0L),
				T( 0.38183005050511894495036977548897e0L),
				T( 0.27970539148927666790146777142377e0L),
				T( 0.12948496616886969327061143267908e0L)
			};

			static const T* abscissae[] =
			{
				a1, a2, a3, a4, a5, a6, a7
			};

			static const T* weights[] =
			{
				w1, w2, w3, w4, w5, w6, w7
			};

			if (shape_type == shape_linear)
			{
				integration_rule.Resize(points_count);
				for (int i = 1; i <= points_count; ++i)
				{
					integration_rule.entry[i].point.Resize(1);
					integration_rule.entry[i].point.entry[1] = abscissae[points_count - 1][i - 1];
					integration_rule.entry[i].weight = weights[points_count - 1][i - 1];
				}
			}
			else if (shape_type == shape_quadrilateral)
			{
				integration_rule.Resize(points_count*points_count);
				int o = points_count - 1;
				int i = 0;
				for (int d1 = 0; d1 < points_count; ++d1)
				{
					for (int d2 = 0; d2 < points_count; ++d2)
					{
						++i;
						integration_rule.entry[i].point.Resize(2);
						integration_rule.entry[i].point.entry[1] = abscissae[o][d1];
						integration_rule.entry[i].point.entry[2] = abscissae[o][d2];
						integration_rule.entry[i].weight = weights[o][d1]*weights[o][d2];
					}
				}
			}
			else // shape_hexahedron
			{
				integration_rule.Resize(points_count*points_count*points_count);
				int o = points_count - 1;
				int i = 0;
				for (int d1 = 0; d1 < points_count; ++d1)
				{
					for (int d2 = 0; d2 < points_count; ++d2)
					{
						for (int d3 = 0; d3 < points_count; ++d3)
						{
							++i;
							integration_rule.entry[i].point.Resize(3);
							integration_rule.entry[i].point.entry[1] = abscissae[o][d1];
							integration_rule.entry[i].point.entry[2] = abscissae[o][d2];
							integration_rule.entry[i].point.entry[3] = abscissae[o][d3];
							integration_rule.entry[i].weight = weights[o][d1]*weights[o][d2]*weights[o][d3];
						}
					}
				}
			}
		}
		else if (shape_type == shape_triangle)
		{
			Assert((points_count == 1) || (points_count == 3) || (points_count == 4) || (points_count == 7));

			// J. E. Akin. Finite Element Analysis with Error Estimators. Elsevier Butterworth-Heinemann. 2005. p 271.

			// Abscissae and weights with 1 point
			static const T a1[] =
			{
				T( 0.33333333333333333333333333333333e0L), T( 0.33333333333333333333333333333333e0L)
			};
			static const T w1[] =
			{
				T( 0.5e0L)
			};

			// Abscissae and weights with 3 points
			static const T a3[] =
			{
				T( 0.16666666666666666666666666666667e0L), T( 0.16666666666666666666666666666667e0L),
				T( 0.66666666666666666666666666666667e0L), T( 0.16666666666666666666666666666667e0L),
				T( 0.16666666666666666666666666666667e0L), T( 0.66666666666666666666666666666667e0L)
			};
			static const T w3[] =
			{
				T( 0.16666666666666666666666666666667e0L),
				T( 0.16666666666666666666666666666667e0L),
				T( 0.16666666666666666666666666666667e0L)
			};

			// Abscissae and weights with 4 points
			static const T a4[] =
			{
				T( 0.33333333333333333333333333333333e0L), T( 0.33333333333333333333333333333333e0L),
				T( 0.60000000000000000000000000000000e0L), T( 0.20000000000000000000000000000000e0L),
				T( 0.20000000000000000000000000000000e0L), T( 0.60000000000000000000000000000000e0L),
				T( 0.20000000000000000000000000000000e0L), T( 0.20000000000000000000000000000000e0L)
			};
			static const T w4[] =
			{
				T(-0.28125000000000000000000000000000e0L),
				T( 0.26041666666666666666666666666667e0L),
				T( 0.26041666666666666666666666666667e0L),
				T( 0.26041666666666666666666666666667e0L)
			};

			// Abscissae and weights with 7 points
			static const T a7[] =
			{
				T( 0.00000000000000000000000000000000e0L), T( 0.00000000000000000000000000000000e0L),
				T( 0.50000000000000000000000000000000e0L), T( 0.00000000000000000000000000000000e0L),
				T( 1.00000000000000000000000000000000e0L), T( 0.00000000000000000000000000000000e0L),
				T( 0.50000000000000000000000000000000e0L), T( 0.50000000000000000000000000000000e0L),
				T( 0.00000000000000000000000000000000e0L), T( 1.00000000000000000000000000000000e0L),
				T( 0.00000000000000000000000000000000e0L), T( 0.50000000000000000000000000000000e0L),
				T( 0.33333333333333333333333333333333e0L), T( 0.33333333333333333333333333333333e0L)
			};
			static const T w7[] =
			{
				T( 0.02500000000000000000000000000000e0L),
				T( 0.06666666666666666666666666666667e0L),
				T( 0.02500000000000000000000000000000e0L),
				T( 0.06666666666666666666666666666667e0L),
				T( 0.02500000000000000000000000000000e0L),
				T( 0.06666666666666666666666666666667e0L),
				T( 0.22500000000000000000000000000000e0L)
			};

			static const T* abscissae[] =
			{
				a1, (const T*)0, a3, a4, (const T*)0, (const T*)0, a7
			};

			static const T* weights[] =
			{
				w1, (const T*)0, w3, w4, (const T*)0, (const T*)0, w7
			};

			int o = points_count - 1;
			integration_rule.Resize(points_count);
			for (int i = 1; i <= points_count; ++i)
			{
				int j = (i - 1)*2;
				integration_rule.entry[i].point.Resize(2);
				integration_rule.entry[i].point.entry[1] = abscissae[o][j];
				integration_rule.entry[i].point.entry[2] = abscissae[o][j + 1];
				integration_rule.entry[i].weight = weights[o][i - 1];
			}
		}
		else if (shape_type == shape_tetrahedron)
		{
			Assert((points_count == 1) || (points_count == 4) || (points_count == 5) || (points_count == 11));

			// J. E. Akin. Finite Element Analysis with Error Estimators. Elsevier Butterworth-Heinemann. 2005. p 272.

			// Abscissae and weights with 1 point
			static const T a1[] =
			{
				T( 0.25000000000000000000000000000000e0L), T( 0.25000000000000000000000000000000e0L), T( 0.25000000000000000000000000000000e0L)
				
			};
			static const T w1[] =
			{
				T( 0.16666666666666666666666666666667e0L)
			};

			// Abscissae and weights with 4 points
			static const T a4[] =
			{
				T( 0.58541019662496845446137605030969e0L), T( 0.13819660112501051517954131656344e0L), T( 0.13819660112501051517954131656344e0L),
				T( 0.13819660112501051517954131656344e0L), T( 0.58541019662496845446137605030969e0L), T( 0.13819660112501051517954131656344e0L),
				T( 0.13819660112501051517954131656344e0L), T( 0.13819660112501051517954131656344e0L), T( 0.58541019662496845446137605030969e0L),
				T( 0.13819660112501051517954131656344e0L), T( 0.13819660112501051517954131656344e0L), T( 0.13819660112501051517954131656344e0L)
			};
			static const T w4[] =
			{
				T( 0.04166666666666666666666666666667e0L),
				T( 0.04166666666666666666666666666667e0L),
				T( 0.04166666666666666666666666666667e0L),
				T( 0.04166666666666666666666666666667e0L)
			};

			// Abscissae and weights with 5 points
			static const T a5[] =
			{
				T( 0.25000000000000000000000000000000e0L), T( 0.25000000000000000000000000000000e0L), T( 0.25000000000000000000000000000000e0L),
				T( 0.50000000000000000000000000000000e0L), T( 0.16666666666666666666666666666667e0L), T( 0.16666666666666666666666666666667e0L),
				T( 0.16666666666666666666666666666667e0L), T( 0.50000000000000000000000000000000e0L), T( 0.16666666666666666666666666666667e0L),
				T( 0.16666666666666666666666666666667e0L), T( 0.16666666666666666666666666666667e0L), T( 0.50000000000000000000000000000000e0L),
				T( 0.16666666666666666666666666666667e0L), T( 0.16666666666666666666666666666667e0L), T( 0.16666666666666666666666666666667e0L)
			};
			static const T w5[] =
			{
				T(-0.13333333333333333333333333333333e0L),
				T( 0.07500000000000000000000000000000e0L),
				T( 0.07500000000000000000000000000000e0L),
				T( 0.07500000000000000000000000000000e0L),
				T( 0.07500000000000000000000000000000e0L)
			};

			// Abscissae and weights with 11 points
			static const T a11[] =
			{
				T( 0.25000000000000000000000000000000e0L), T( 0.25000000000000000000000000000000e0L), T( 0.25000000000000000000000000000000e0L),
				T( 0.78571428571428571428571428571429e0L), T( 0.07142857142857142857142857142857e0L), T( 0.07142857142857142857142857142857e0L),
				T( 0.07142857142857142857142857142857e0L), T( 0.78571428571428571428571428571429e0L), T( 0.07142857142857142857142857142857e0L),
				T( 0.07142857142857142857142857142857e0L), T( 0.07142857142857142857142857142857e0L), T( 0.78571428571428571428571428571429e0L),
				T( 0.07142857142857142857142857142857e0L), T( 0.07142857142857142857142857142857e0L), T( 0.07142857142857142857142857142857e0L),
				T( 0.39940357616679920499610214746164e0L), T( 0.39940357616679920499610214746164e0L), T( 0.10059642383320079500389785253836e0L),
				T( 0.39940357616679920499610214746164e0L), T( 0.10059642383320079500389785253836e0L), T( 0.39940357616679920499610214746164e0L),
				T( 0.39940357616679920499610214746164e0L), T( 0.10059642383320079500389785253836e0L), T( 0.10059642383320079500389785253836e0L),
				T( 0.10059642383320079500389785253836e0L), T( 0.39940357616679920499610214746164e0L), T( 0.39940357616679920499610214746164e0L),
				T( 0.10059642383320079500389785253836e0L), T( 0.39940357616679920499610214746164e0L), T( 0.10059642383320079500389785253836e0L),
				T( 0.10059642383320079500389785253836e0L), T( 0.10059642383320079500389785253836e0L), T( 0.39940357616679920499610214746164e0L)
			};
			static const T w11[] =
			{
				T(-0.01315555555555555555555555555556e0L),
				T( 0.00762222222222222222222222222222e0L),
				T( 0.00762222222222222222222222222222e0L),
				T( 0.00762222222222222222222222222222e0L),
				T( 0.00762222222222222222222222222222e0L),
				T( 0.02488888888888888888888888888889e0L),
				T( 0.02488888888888888888888888888889e0L),
				T( 0.02488888888888888888888888888889e0L),
				T( 0.02488888888888888888888888888889e0L),
				T( 0.02488888888888888888888888888889e0L),
				T( 0.02488888888888888888888888888889e0L)
			};

			static const T* abscissae[] =
			{
				a1, (const T*)0, (const T*)0, a4, a5, (const T*)0, (const T*)0, (const T*)0, (const T*)0, (const T*)0, a11
			};

			static const T* weights[] =
			{
				w1, (const T*)0, (const T*)0, w4, w5, (const T*)0, (const T*)0, (const T*)0, (const T*)0, (const T*)0, w11
			};

			int o = points_count - 1;
			integration_rule.Resize(points_count);
			for (int i = 1; i <= points_count; ++i)
			{
				int j = (i - 1)*3;
				integration_rule.entry[i].point.Resize(3);
				integration_rule.entry[i].point.entry[1] = abscissae[o][j];
				integration_rule.entry[i].point.entry[2] = abscissae[o][j + 1];
				integration_rule.entry[i].point.entry[3] = abscissae[o][j + 2];
				integration_rule.entry[i].weight = weights[o][i - 1];
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}

#endif
