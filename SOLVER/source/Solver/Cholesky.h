// Cholesky.h
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

#ifndef _Cholesky_h_
#define _Cholesky_h_

#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Container/CSRMatrix.h>
#include <Container/Set.h>
#include <Container/TriangularMatrix.h>
#include <Solver/DiagonalSystem.h>
#include <Solver/Solver.h>
#include <Solver/TriangularSystem.h>
#include <math.h>


extern "C"
{
	#include <metis.h>
}


template <typename T>
void FindingAnOrdering(const CSRMatrix<T>& A, Vector<int>& index, Vector<int>& inverse_index) throw(Memory::Exception)
{
	try
	{
		// Create adjacency structures for METIS
		int N = A.rows;
		int total_edges = A.NonZero() - N;
		Vector<int> adjacency_start(N + 1);
		Vector<int> adjacency_nodes(total_edges);
		for (int n = 1, a = 1; n <= N; ++n)
		{
			adjacency_start.entry[n] = a;
			int k_max = A.Count(n);
			for (int k = 1; k <= k_max; ++k)
			{
				int ne = A.index[n][k];
				if (ne != n)
				{			
					adjacency_nodes.entry[a] = ne;
					++a;
				}
			}
		}
		adjacency_start.entry[N + 1] = total_edges + 1;

		// Call reordering routine
		int numflag = 1; // Used to indicate which numbering scheme is used for the adjacency structure of the graph: 0 C-style, 1 Fortran-style
		int options[8] = {0, 0, 0, 0, 0, 0, 0, 0};
		index.Resize(N);
		inverse_index.Resize(N);
		METIS_NodeND(&N, adjacency_start.data, adjacency_nodes.data, &numflag, options, index.data, inverse_index.data);
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void Reorder(const Vector<int>& index, const Vector<int>& inverse_index, const CSRMatrix<T>& A, CSRMatrix<T>& Ar, int threads) throw(Memory::Exception)
{
	Assert(Ar.rows == A.rows);
	Assert(Ar.columns == A.columns);

	try
	{
		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int n = 1; n <= Ar.rows; ++n)
		{
			int g = index.entry[n];
			Ar.AllocateRow(n, A.Count(g));
			int k_max = Ar.Count(n);
			for (int k = 1; k <= k_max; ++k)
			{
				int rg = inverse_index.entry[A.index[g][k]];
				Ar.index[n][k] = rg;
				Ar.entry[n][k] = A.entry[g][k];
			}
			Ar.SortByIndexes(n);
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void Reorder(const Vector<int>& index, const Vector<T>& v, Vector<T>& vr, int threads) throw()
{
	#pragma omp parallel for default(shared) num_threads(threads)
	for (int n = 1; n <= v.size; ++n)
	{
		int g = index.entry[n];
		vr.entry[n] = v.entry[g];
	}
}


template <typename T>
void SymbolicCholeskyDecomposition(const CSRMatrix<T>& A, CSRMatrix<T>& L, CSRMatrix<T>& Lt, int ick = -1) throw(Memory::Exception)
{
	// Symbolic Cholesky decomposition (Parallel algorithms for matrix computations. Gallivan)

	Assert(L.rows == A.rows);
	Assert(L.columns == A.columns);
	Assert(Lt.rows == A.columns);
	Assert(Lt.columns == A.rows);

	try
	{
		Log(1, "SymbolicCholeskyDecomposition:");
		Log(1, "-Level:   %i", ick);

		if (ick == -1)
		{
			ick = A.rows;
		}

		Vector<int> L_count(A.rows);
		L_count.Fill(0);
		{
			Vector<int> index(A.rows);
			Vector<bool> used(A.rows);
			used.Fill(false);
			Vector<Set<int, 16> > R(A.rows);
			for (int j = 1; j <= A.rows; ++j)
			{
				int* __restrict A_index_j = A.index[j];

				int N = 0;
				int p = A.rows + 1;
				int k_max = A.Count(j);
				for (int k = 1; k <= k_max; ++k)
				{
					register int i = A_index_j[k];
					if (i > j)
					{
						if (!used.entry[i])
						{
							++N;
							index.entry[N] = i;
							used.entry[i] = true;
							if (i < p)
							{
								p = i;
							}
						}
					}
				}
				if (ick > 0)
				{
					for (SetItem<int>* r = R.entry[j].first; r; r = r->next)
					{
						int i = r->value;
						int* __restrict Lt_index_i = Lt.index[i];

						int Lt_count_i = Lt.Count(i);
						for (int k = 2; k <= Lt_count_i; ++k)
						{
							register int l = Lt_index_i[k];
							if ((l - j) > ick)
							{
								r = R.entry[j].last;
								break;
							}
							if (l != j)
							{
								if (!used.entry[l])
								{
									++N;
									index.entry[N] = l;
									used.entry[l] = true;
									if (l < p)
									{
										p = l;
									}
								}
							}
						}
					}
				}
				R.entry[(p == A.rows + 1) ? j : p].Append(j);

				// Generate L'(j)
				++N;
				index.entry[N] = j;
				Lt.AllocateRow(j, N);
				for (register int c = 1; c <= N; ++c)
				{
					int i = index.entry[c];
					Lt.index[j][c] = i;
					used.entry[i] = false;
					++L_count.entry[i];
				}
				Lt.SortIndexes(j);

				if (log_level >= 2)
				{
					if (!(j % 1000))
					{
						Log(2, "-Column: %7i", j);
					}
				}
			}
		}
		Log(1, "-nnz(L'): %i", Lt.NonZero());

		// Generate L
		for (int i = 1; i <= L.rows; ++i)
		{
			L.AllocateRow(i, L_count.entry[i]);
		}
		for (int j = 1; j <= Lt.columns; ++j)
		{
			int* __restrict Lt_index_j = Lt.index[j];
			int Lt_count_j = Lt.Count(j);
			for (register int c = 1; c <= Lt_count_j; ++c)
			{
				register int i = Lt_index_j[c];
				L.index[i][L_count.entry[i]--] = j;
			}
		}
		L.SortIndexes();
		Log(1, "-nnz(L):  %i", Lt.NonZero());
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void FillCholeskyDecomposition(const CSRMatrix<T>& A, CSRMatrix<T>& L, CSRMatrix<T>& Lt, int threads) throw()
{
	Log(1, "FillCholeskyDecomposition:");

	// Generate LL' decomposition
	for (int j = 1; j <= L.columns; ++j)
	{
		T* __restrict L_entry_j = L.entry[j];
		T* __restrict Lt_entry_j = Lt.entry[j];

		T L_jj = A(j, j);
		int L_count_j = L.Count(j);
		for (register int q = 1; q < L_count_j; ++q)
		{
			L_jj -= L_entry_j[q]*L_entry_j[q];
		}
		L_jj = sqrt(L_jj);
		L_entry_j[L_count_j] = L_jj; // L(j)(j) = sqrt(A(j)(j) - sum(k = 1, j - 1, L(j)(k)*L(j)(k)))
		Lt_entry_j[1] = L_jj;

		int Lt_count_j = Lt.Count(j);
		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int q = 2; q <= Lt_count_j; ++q)
		{
			int i = Lt.index[j][q];

			T* __restrict L_entry_i = L.entry[i];

			T L_ij = A(i, j);

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
					L_ij -= L_entry_i[qi]*L_entry_j[qj];
					++qi;
					++qj;
					ki = L_index_i[qi];
					kj = L_index_j[qj];
				}
			}
			L_ij /= L_jj;
			L_entry_i[qi] = L_ij;  // L(i)(j) = (A(i)(j) - sum(k = 1, j - 1, L(i)(k)*L(j)(k)))/L(j)(j)
			Lt_entry_j[q] = L_ij;
		}

		if (log_level >= 2)
		{
			if (!(j % 1000))
			{
				Log(2, "-Column: %7i", j);
			}
		}
	}
	Log(1, "-L and L' filled");
}


template <typename T>
class CholeskySolver : public Solver<T>
{
	public:

		CholeskySolver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			L(A.rows, A.columns),
			Lt(A.rows, A.columns),
			c(A.rows)
		{
			try
			{
				SymbolicCholeskyDecomposition(A, L, Lt);
				FillCholeskyDecomposition(A, L, Lt, this->threads);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

		virtual bool Calculate() throw()
		{
			LowerTriangularSystem(L, c, this->b);
			UpperTriangularSystem(Lt, this->x, c);
			return this->CheckSolution();
		}

	private:

		CSRMatrix<T> L;
		CSRMatrix<T> Lt;
		Vector<T> c;
};


template <typename T>
void FillCholesky2Decomposition(const CSRMatrix<T>& A, CSRMatrix<T>& L, Vector<T>& D, CSRMatrix<T>& Lt, int threads) throw()
{
	Log(1, "FillCholesky2Decomposition");

	// Generate LDL' decomposition
	for (int j = 1; j <= L.columns; ++j)
	{
		T* __restrict L_entry_j = L.entry[j];
		T* __restrict Lt_entry_j = Lt.entry[j];
		const register int* __restrict L_index_j = L.index[j];

		T D_j = A(j, j);
		int L_count_j = L.Count(j);
		for (register int q = 1; q < L_count_j; ++q)
		{
			D_j -= L_entry_j[q]*L_entry_j[q]*D.entry[L_index_j[q]];
		}
		D.entry[j] = D_j; // D(j) = A(j)(j) - sum(k = 1, j - 1, L(j)(k)*L(j)(k)*D(k))
		L_entry_j[L_count_j] = 1; // L(j)(j) = 1
		Lt_entry_j[1] = 1;

		int Lt_count_j = Lt.Count(j);
		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int q = 2; q <= Lt_count_j; ++q)
		{
			int i = Lt.index[j][q];

			T* __restrict L_entry_i = L.entry[i];

			T L_ij = A(i, j);

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
					L_ij -= L_entry_i[qi]*L_entry_j[qj]*D.entry[ki];
					++qi;
					++qj;
					ki = L_index_i[qi];
					kj = L_index_j[qj];
				}
			}
			L_ij /= D_j;
			L_entry_i[qi] = L_ij;  // L(i)(j) = (A(i)(j) - sum(k = 1, j - 1, L(i)(k)*L(j)(k)))/L(j)(j)
			Lt_entry_j[q] = L_ij;
		}

		if (log_level >= 2)
		{
			if (!(j % 1000))
			{
				Log(2, "-Column: %7i", j);
			}
		}
	}
	Log(1, "-L, D and L' filled");
}


template <typename T>
void FillCholesky2Decomposition(const CSRMatrix<T>& A, CSRMatrix<T>& L, Vector<T>& D, CSRMatrix<T>& Lt, int threads, T alpha, T u) throw()
{
	Log(1, "FillCholesky2Decomposition, alpha = %g, u = %g", alpha, u);

	// Generate LDL' decomposition
	for (int j = 1; j <= L.columns; ++j)
	{
		T* __restrict L_entry_j = L.entry[j];
		T* __restrict Lt_entry_j = Lt.entry[j];
		const register int* __restrict L_index_j = L.index[j];

		T D_j = alpha*A(j, j);
		int L_count_j = L.Count(j);
		for (register int q = 1; q < L_count_j; ++q)
		{
			D_j -= L_entry_j[q]*L_entry_j[q]*D.entry[L_index_j[q]];
		}

		// N. Munksgaard. Solving sparse symmetric sets of linear equations by preconditioned conjugate gradients. ACM Transactions on Mathematical Software, Vol 6-2, pp. 206-219. 1980
		{
			T* __restrict A_entry_j = A.entry[j];
			int* __restrict A_index_j = A.index[j];
			register T Sum_A_j = 0;
			int A_count_j = A.Count(j);
			for (register int q = 1; q <= A_count_j; ++q)
			{
				if (A_index_j[q] != j)
				{
					Sum_A_j += A_entry_j[q];
				}
			}
			Sum_A_j = fabs(Sum_A_j);
			if (D_j < u*Sum_A_j)
			{
				D_j = (Sum_A_j == 0) ? 1 : Sum_A_j;
			}
		}

		D.entry[j] = D_j; // D(j) = A(j)(j) - sum(k = 1, j - 1, L(j)(k)*L(j)(k)*D(k))
		L_entry_j[L_count_j] = 1; // L(j)(j) = 1
		Lt_entry_j[1] = 1;

		int Lt_count_j = Lt.Count(j);
		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int q = 2; q <= Lt_count_j; ++q)
		{
			int i = Lt.index[j][q];

			T* __restrict L_entry_i = L.entry[i];

			T L_ij = A(i, j);

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
					L_ij -= L_entry_i[qi]*L_entry_j[qj]*D.entry[ki];
					++qi;
					++qj;
					ki = L_index_i[qi];
					kj = L_index_j[qj];
				}
			}
			L_ij /= D_j;
			L_entry_i[qi] = L_ij;  // L(i)(j) = (A(i)(j) - sum(k = 1, j - 1, L(i)(k)*L(j)(k)))/L(j)(j)
			Lt_entry_j[q] = L_ij;
		}

		if (log_level >= 2)
		{
			if (!(j % 1000))
			{
				Log(2, "-Column: %7i", j);
			}
		}
	}
	Log(1, "-L, D and L' filled");
}


template <typename T>
class Cholesky2Solver : public Solver<T>
{
	public:

		Cholesky2Solver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			L(A.rows, A.columns),
			D(A.rows),
			Lt(A.rows, A.columns),
			c(A.rows),
			d(A.rows)
		{
			try
			{
				SymbolicCholeskyDecomposition(A, L, Lt);
				FillCholesky2Decomposition(A, L, D, Lt, this->threads);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

		virtual bool Calculate() throw()
		{
			LowerTriangularSystem(L, c, this->b);
			DiagonalSystem(D, d, c);
			UpperTriangularSystem(Lt, this->x, d);
			return this->CheckSolution();
		}

	private:

		CSRMatrix<T> L;
		Vector<T> D;
		CSRMatrix<T> Lt;
		Vector<T> c;
		Vector<T> d;
};

#endif
