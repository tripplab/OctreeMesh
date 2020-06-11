// Solver.h
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

#ifndef _Solver_h_
#define _Solver_h_

#include <Basic/Memory.h>
#include <Basic/Float.h>
#include <Container/CSRMatrix.h>
#include <Container/Vector.h>


template <typename T>
class Solver
{
	public:

		CSRMatrix<T>& A;
		Vector<T>& x;
		Vector<T>& b;
		Vector<bool>& fixed;
		int threads;


		Solver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	A(A),
			x(x),
			b(b),
			fixed(fixed),
			threads(threads),
			matrix_adjust(matrix_adjust),
			b_adjust(matrix_adjust ? 0 : A.rows),
			A_adjust(matrix_adjust ? A.rows : 0, matrix_adjust ? A.rows : 0)
		{
			try
			{
				if (matrix_adjust)
				{
					// Keeping A_adjust
					for (int i = 1; i <= A.rows; ++i)
					{
						if (!fixed.entry[i])
						{
							int* __restrict A_index_i = A.index[i];
							T* __restrict A_entry_i = A.entry[i];

							int* __restrict A_adjust_index_i = A_adjust.index[i];
							T* __restrict A_adjust_entry_i = A_adjust.entry[i];

							int fixed_count = 0;
							int k_max = A.Count(i);
							for (int k = 1; k <= k_max; ++k)
							{
								register int j = A_index_i[k];
								if (fixed.entry[j])
								{
									++fixed_count;
								}
							}
							A_adjust.AllocateRow(i, fixed_count);
							A_adjust_index_i = A_adjust.index[i];
							A_adjust_entry_i = A_adjust.entry[i];
							for (int k = 1, bk = 1; k <= k_max; ++k)
							{
								register int j = A_index_i[k];
								if (fixed.entry[j])
								{
									A_adjust_index_i[bk] = j;
									A_adjust_entry_i[bk] = A_entry_i[k];
									++bk;
								}
							}
						}
					}
				}
				else
				{
					// Keeping fixed conditions in b_adjust
					b_adjust.Fill(T(0));
					for (int i = 1; i <= A.rows; ++i)
					{
						if (!fixed.entry[i])
						{
							int* __restrict A_index_i = A.index[i];
							T* __restrict A_entry_i = A.entry[i];

							int k_max = A.Count(i);
							for (int k = 1; k <= k_max; ++k)
							{
								register int j = A_index_i[k];
								if (fixed.entry[j])
								{
									b_adjust.entry[i] += A_entry_i[k]*x.entry[j];
								}
							}
						}
					}
				}

				// Remove fixed entries from A
				for (int i = 1; i <= A.rows; ++i)
				{
					if (fixed.entry[i])
					{
						int* __restrict A_index_i = A.index[i];

						int k_max = A.Count(i);
						for (int k = 1; k <= k_max; ++k)
						{
							int j = A_index_i[k];
							if (j != i)
							{
								A.RemoveEntry(j, i);
							}
						}
						A.AllocateRow(i, 1);
						A.index[i][1] = i;
						A.entry[i][1] = 1.0;
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		virtual ~Solver()
		{
		}


		void CompensateFixed() throw()
		{
			// Compensate vector b considering fixed entries
			if (matrix_adjust)
			{
				#pragma omp parallel for default(shared) num_threads(threads)
				for (int i = 1; i <= A_adjust.rows; ++i)
				{
					if (!fixed.entry[i])
					{
						int* __restrict A_adjust_index_i = A_adjust.index[i];
						T* __restrict A_adjust_entry_i = A_adjust.entry[i];

						int k_max = A_adjust.Count(i);
						T sum = 0;
						for (register int k = 1; k <= k_max; ++k)
						{
							register int j = A_adjust_index_i[k];
							sum += A_adjust_entry_i[k]*x.entry[j];
						}
						b.entry[i] -= sum;
					}
					else
					{
						b.entry[i] = x.entry[i];
					}
				}
			}
			else
			{
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
		}


		virtual bool Calculate() throw(Memory::Exception)
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
		CSRMatrix<T> A_adjust; // Will store the fixed enties needed to initilize the b vector (this method is needed for domain decomposition iterations or when values on fixed conditions change)

		Solver& operator = (const Solver&)
		{
			return *this;
		}
};

#endif
