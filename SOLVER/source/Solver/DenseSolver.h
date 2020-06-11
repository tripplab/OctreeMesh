// DenseSolver.h
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

#ifndef _DenseSolver_h_
#define _DenseSolver_h_

#include <Basic/Memory.h>
#include <Basic/Float.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>


template <typename T>
class DenseSolver
{
	public:

		Matrix<T>& A;
		Vector<T>& x;
		Vector<T>& b;
		Vector<bool>& fixed;
		int threads;


		DenseSolver(Matrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, int threads) throw(Memory::Exception)
		:	A(A),
			x(x),
			b(b),
			fixed(fixed),
			threads(threads),
			b_adjust(A.rows)
		{
			// Keeping fixed conditions in b_adjust
			b_adjust.Fill(T(0));
			#pragma omp parallel for default(shared) num_threads(threads)
			for (int i = 1; i <= A.rows; ++i)
			{
				T* __restrict A_entry_i = A.entry[i];
				if (!fixed.entry[i])
				{
					for (int j = 1; j <= A.columns; ++j)
					{
						if (fixed.entry[j])
						{
							b_adjust.entry[i] += A_entry_i[j]*x.entry[j];
						}
					}
				}
			}

			// Remove fixed entries from A
			for (int i = 1; i <= A.rows; ++i)
			{
				if (fixed.entry[i])
				{
					for (int j = 1; j <= A.columns; ++j)
					{
						A.entry[i][j] = (T)0;
					}
				}
			}
			for (int j = 1; j <= A.columns; ++j)
			{
				if (fixed.entry[j])
				{
					for (int i = 1; i <= A.rows; ++i)
					{
						A.entry[i][j] = (T)0;
					}
					A.entry[j][j] = (T)1;
				}
			}
		}


		virtual ~DenseSolver()
		{
		}


		void CompensateFixed() throw()
		{
			// Compensate vector b considering fixed entries
			#pragma omp parallel for default(shared) num_threads(threads)
			for (int i = 1; i <= b.size; ++i)
			{
				if (!fixed.entry[i])
				{
					b.entry[i] -= b_adjust.entry[i];
				}
				else
				{
					b.entry[i] = x.entry[i];
				}
			}
		}


		virtual bool Calculate() throw()
		{
			return false;
		}


		bool CheckSolution() throw()
		{
			// Check for a valid solution
			int invalid = 0;
			#pragma omp parallel for default(shared) reduction(+:invalid) num_threads(threads)
			for (int i = 1; i <= x.size; ++i)
			{
				if ((x.entry[i] == Float<T>::infinite) || (x.entry[i] == -Float<T>::infinite) || Float<T>::IsNaN(x.entry[i]))
				{
					++invalid;
				}
			}
			return (invalid > 0) ? false : true;
		}


	private:

		bool matrix_adjust;
		Vector<T> b_adjust; // Will store the fixed enties needed to initilize the b vector (important when calling the solver several times)

		DenseSolver& operator = (const DenseSolver&)
		{
			return *this;
		}
};

#endif
