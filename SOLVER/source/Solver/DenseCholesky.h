// DenseCholesky.h
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

#ifndef _DenseCholesky_h_
#define _DenseCholesky_h_

#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Container/Matrix.h>
#include <Container/TriangularMatrix.h>
#include <Solver/DiagonalSystem.h>
#include <Solver/DenseSolver.h>
#include <Solver/DenseTriangularSystem.h>
#include <math.h>


template <typename T>
class DenseCholeskySolver : public DenseSolver<T>
{
	public:

		DenseCholeskySolver(Matrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, int threads) throw(Memory::Exception)
		:	DenseSolver<T>(A, x, b, fixed, threads),
			L(A.rows, A.columns),
			Lt(A.rows, A.columns),
			c(A.rows)
		{
			DenseCholesky(A, L, Lt, this->threads);
		}

		virtual bool Calculate() throw()
		{
			LowerTriangularSystem(L, c, this->b);
			UpperTriangularSystem(Lt, this->x, c);
			return this->CheckSolution();
		}

	private:

		LowerTriangularMatrix<T> L;
		UpperTriangularMatrix<T> Lt;
		Vector<T> c;
};


template <typename T>
void DenseCholesky(const Matrix<T>& A, LowerTriangularMatrix<T>& L, UpperTriangularMatrix<T>& Lt, int threads) throw()
{
	Assert(A.rows == A.columns);
	Assert(L.rows == A.rows);
	Assert(L.columns == A.columns);
	Assert(Lt.rows == L.rows);
	Assert(Lt.columns == L.columns);

	Log(1, "DenseCholeskyDecomposition:");

	// Generate LL' decomposition
	for (int j = 1; j <= L.columns; ++j)
	{
		T* __restrict L_entry_j = L.entry[j];
		T* __restrict Lt_entry_j = Lt.entry[j];

		T L_jj = A.entry[j][j];
		for (register int k = 1; k < j; ++k)
		{
			L_jj -= L_entry_j[k]*L_entry_j[k];
		}
		L_jj = sqrt(L_jj);
		L_entry_j[j] = L_jj; // L(j)(j) = sqrt(A(j)(j) - sum(k = 1, j - 1, L(j)(k)*L(j)(k)))
		Lt_entry_j[j] = L_jj;

		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int i = j + 1; i <= L.rows; ++i)
		{
			T* __restrict L_entry_i = L.entry[i];

			register T L_ij = A.entry[i][j];
			for (register int k = 1; k < j; ++k)
			{
				L_ij -= L_entry_i[k]*L_entry_j[k];
			}
			L_ij /= L_jj;
			L_entry_i[j] = L_ij;  // L(i)(j) = (A(i)(j) - sum(k = 1, j - 1, L(i)(k)*L(j)(k)))/L(j)(j)
			Lt_entry_j[i] = L_ij;
		}

		if (log_level >= 2)
		{
			if (!(j % 100))
			{
				Log(2, "-Column: %7i", j);
			}
		}
	}

	Log(1, "-L and L' filled");
}

#endif
