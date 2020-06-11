// LU.h
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

#ifndef _LU_h_
#define _LU_h_

#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Container/CSRMatrix.h>
#include <Solver/Cholesky.h>
#include <Solver/Solver.h>
#include <Solver/TriangularSystem.h>


// This function only works on symmetric structure LU decomposition
template <typename T>
void FillLUDecomposition(const CSRMatrix<T>& A, CSRMatrix<T>& L, CSRMatrix<T>& U, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "FillLUDecomposition:");

		// Generate the structure of U'
		CSRMatrix<T> Ut(L.rows, L.columns);
		for (int i = 1; i <= Ut.rows; ++i)
		{
			int count_i = L.Count(i);
			Ut.AllocateRow(i, count_i);

			int* __restrict Ut_index_i = Ut.index[i];
			int* __restrict L_index_i = L.index[i];
			for (register int k = 1; k <= count_i; ++k)
			{
				Ut_index_i[k] = L_index_i[k];
			}
		}
		Log(1, "-U' structure created");

		// Generate LU decomposition
		for (int j = 1; j <= L.columns; ++j)
		{
			T* __restrict L_entry_j = L.entry[j];
			T* __restrict U_entry_j = U.entry[j];
			T* __restrict Ut_entry_j = Ut.entry[j];

			T U_jj = A(j, j);
			int L_count_j = L.Count(j);
			for (register int q = 1; q < L_count_j; ++q)
			{
				U_jj -= L_entry_j[q]*Ut_entry_j[q];
			}
			L_entry_j[L_count_j] = 1; // L(j)(j) = 1
			U_entry_j[1] = U_jj;     // U(j)(j) = (A(j)(j) - sum(k = 1, j - 1, L(j)(k)*U(k)(j)))/L(j)(j)
			Ut_entry_j[L_count_j] = U_jj;

			int U_count_j = U.Count(j);
			#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
			for (int q = 2; q <= U_count_j; ++q)
			{
				int i = U.index[j][q];

				T* __restrict L_entry_i = L.entry[i];
				T* __restrict Ut_entry_i = Ut.entry[i];

				T L_ij = A(i, j);
				T U_ji = A(j, i);

				const register int* __restrict L_index_j = L.index[j];
				const register int* __restrict L_index_i = L.index[i];

				register int qi = 1;
				register int qj = 1;
				register int ki = L_index_i[qi];
				register int kj = L_index_j[qj];
				for (bool next = true; next; )
				{
					while (ki < kj)
					{
						++qi;
						ki = L_index_i[qi];
					}
					while (ki > kj)
					{
						++qj;
						kj = L_index_j[qj];
					}
					while (ki == kj)
					{
						if (ki == j)
						{
							next = false;
							break;
						}
						L_ij -= L_entry_i[qi]*Ut_entry_j[qj];
						U_ji -= L_entry_j[qj]*Ut_entry_i[qi];
						++qi;
						++qj;
						ki = L_index_i[qi];
						kj = L_index_j[qj];
					}
				}
				L_ij /= U_jj;
				L_entry_i[qi] = L_ij; // L(i)(j) = (A(i)(j) - sum(k = 1, j - 1, L(i)(k)*U(k)(j)))/U(i)(i)
				U_entry_j[q] = U_ji;  // U(j)(i) = (A(j)(i) - sum(k = 1, i - 1, L(j)(k)*U(k)(i)))/L(i)(i)
				Ut_entry_i[qi] = U_ji;
			}

			if (log_level >= 2)
			{
				if (!(j % 1000))
				{
					Log(2, "-Column: %7i", j);
				}
			}
		}
		Log(1, "-L, U and U' filled");
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
class LUSolver : public Solver<T>
{
	public:

		LUSolver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			L(A.rows, A.columns),
			U(A.rows, A.columns),
			c(A.rows)
		{
			try
			{
				SymbolicCholeskyDecomposition(A, L, U);
				FillLUDecomposition(A, L, U, this->threads);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

		virtual bool Calculate() throw()
		{
			LowerTriangularSystem(L, c, this->b);
			UpperTriangularSystem(U, this->x, c);
			return this->CheckSolution();
		}


	private:

		CSRMatrix<T> L;
		CSRMatrix<T> U;
		Vector<T> c;
};

#endif
