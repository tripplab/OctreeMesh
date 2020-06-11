// SparseApproximateInverse.h
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

#ifndef _SparseApproximateInverse_h_
#define _SparseApproximateInverse_h_

#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Container/CSRMatrix.h>
#include <Container/CSVector.h>
#include <Container/List.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>

#include <math.h>
#include <omp.h>


template <typename T>
void PowerPreconditionerPattern(const CSRMatrix<T>& A, CSRMatrix<T>& M, T threshold, int preconditioner_level, bool factorized, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "PowerPreconditionerPattern, threshold = %e, level = %i, factorized = %i", threshold, preconditioner_level, factorized ? 1 : 0);
		// A must have a symetric structure

		int N = A.rows;

		// Scalling diagonal
		Vector<T> Ds(N);
		for (int i = 1; i <= N; ++i)
		{
			T v = fabs(A(i, i));
			Ds.entry[i] = v > 0 ? 1/sqrt(v) : 1;
		}

		// Create thrA (threshoding A)
		int nnz_thrA = 0;
		Vector<CSVector<int> > thrA(N);
		#pragma omp parallel for default(shared) reduction(+:nnz_thrA) schedule(guided) num_threads(threads)
		for (int i = 1; i <= N; ++i)
		{
			List<int> connectivity;
			int A_count_i = A.Count(i);
			for (int c = 1; c <= A_count_i; ++c)
			{
				int j = A.index[i][c];
				if (i != j)
				{
					T v = fabs(Ds.entry[i]*A.entry[i][c]*Ds.entry[j]);
					if (v > threshold)
					{
						connectivity.AppendLast(j);
					}
				}
			}

			thrA.entry[i].Resize(N);
			thrA.entry[i].Allocate(connectivity.size);
			int c = 0;
			for (register ListItem<int>* conection = connectivity.first; conection; conection = conection->next)
			{
				++c;
				register int j = conection->value;
				thrA.entry[i].entry[c] = j;
			}
			nnz_thrA += connectivity.size;
		}
		Log(1, "Matrix threshoding, nnz(thrA) = %i", nnz_thrA);

		// Create pattern of powers sparse matrix
		Vector<Vector<int> > thread_index(threads);
		Vector<Vector<bool> > thread_used(threads);
		for (int t = 1; t <= threads; ++t)
		{
			thread_index.entry[t].Resize(N);
			thread_used.entry[t].Resize(N);
			thread_used.entry[t].Fill(false);
		}
		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int i = 1; i <= N; ++i)
		{
			int thread = omp_get_thread_num() + 1;
			Vector<int>& index = thread_index.entry[thread];
			Vector<bool>& used = thread_used.entry[thread];

			int index_count = thrA.entry[i].count;
			for (int b = 1; b <= index_count; ++b)
			{
				int j = thrA.entry[i].entry[b];
				index.entry[b] = j;
				used.entry[j] = true;
			}
			used.entry[i] = true;

			int new_index_count = index_count;
			for (int l = 1; l <= preconditioner_level; ++l)
			{
				for (int b = 1; b <= index_count; ++b)
				{
					int j = index.entry[b];
					int thrA_count = thrA.entry[j].count;
					for (register int c = 1; c <= thrA_count; ++c)
					{
						register int k = thrA.entry[j].entry[c];
						if (!used.entry[k])
						{
							++new_index_count;
							index.entry[new_index_count] = k;
							used.entry[k] = true;
						}
					}
				}
				index_count = new_index_count;
			}
			++index_count;
			index.entry[index_count] = i;

			// Allocate preconditioner
			if (factorized)
			{
				int M_count_i = 0;
				for (register int b = 1; b <= index_count; ++b)
				{
					int j = index.entry[b];
					used.entry[j] = false;
					if (j <= i)
					{
						++M_count_i;
					}
				}
				M.AllocateRow(i, M_count_i);
				for (register int b = 1, c = 0; b <= index_count; ++b)
				{
					if (index.entry[b] <= i)
					{
						++c;
						M.index[i][c] = index.entry[b];
					}
				}
			}
			else
			{
				M.AllocateRow(i, index_count);
				for (register int b = 1; b <= index_count; ++b)
				{
					int j = index.entry[b];
					used.entry[j] = false;
					M.index[i][b] = index.entry[b];
				}
			}
			used.entry[i] = false;
		}
		M.SortIndexes();
		Log(1, "Pattern calculated, nnz(M) = %i", M.NonZero());
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void SparseApproximateInverseSymmetric(const CSRMatrix<T>& A, CSRMatrix<T>& G, CSRMatrix<T>& Gt, T threshold, int preconditioner_level, int threads) throw(Memory::Exception)
{
	// E. Chow
	// Parallel implementation and practical use of sparse approximate inverse preconditioners with a priori sparsity patterns
	// International Journal of High Performance Computing, Vol 15. pp 56-74. 2001.

	try
	{
		Log(1, "SparseApproximateInverseSymmetric, threshold = %e, level = %i", threshold, preconditioner_level);

		PowerPreconditionerPattern(A, G, threshold, preconditioner_level, true, threads);

		// Allocate memory for Cholesky factorizations
		Vector<Vector<T> > thread_L(threads);
		Vector<Vector<T> > thread_z(threads);
		Vector<Vector<T> > thread_x(threads);
		{
			int N_max = G.Count(1);
			for (register int s = 2; s <= G.rows; ++s)
			{
				if (N_max < G.Count(s))
				{
					N_max = G.Count(s);
				}
			}

			int triangular_size_max = (N_max*(N_max + 1)) >> 1;
			for (register int t = 1; t <= threads; ++t)
			{
				thread_L.entry[t].Resize(triangular_size_max);
				thread_z.entry[t].Resize(triangular_size_max);
				thread_x.entry[t].Resize(triangular_size_max);
			}
		}

		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int s = 1; s <= G.rows; ++s)
		{
			// Solve the small SPD system, then fill the triangular superior part of M
			int N = G.Count(s);
			if (N == 1)
			{
				register int i = G.index[s][1];
				G.entry[s][1] = 1/A(i, i);
			}
			else if (N == 2)
			{
				register int i1 = G.index[s][1];
				register int i2 = G.index[s][2];
				register T a = A(i1, i1);
				register T b = A(i1, i2), c = A(i2, i2);
				register T det = a*c - b*b;
				G.entry[s][1] = -b/det;
				G.entry[s][2] = a/det;
			}
			else if (N == 3)
			{
				register int i1 = G.index[s][1];
				register int i2 = G.index[s][2];
				register int i3 = G.index[s][3];
				T a = A(i1, i1);
				T b = A(i2, i1), c = A(i2, i2);
				T d = A(i3, i1), e = A(i3, i2), f = A(i3, i3);
				T becd = b*e - c*d;
				T dbae = d*b - a*e;
				T acbb = a*c - b*b;
				T det = d*becd + e*dbae + f*acbb;
				G.entry[s][1] = becd/det;
				G.entry[s][2] = dbae/det;
				G.entry[s][3] = acbb/det;
			}
			else
			{
				int thread = omp_get_thread_num() + 1;
				Vector<T>& L = thread_L.entry[thread];
				Vector<T>& z = thread_z.entry[thread];
				Vector<T>& x = thread_x.entry[thread];

				// Cholesky factorization
				for (register int i = 1; i <= N; ++i)
				{
					register int ii = (i*(i - 1)) >> 1;
					T Lii = A(G.index[s][i], G.index[s][i]);
					for (register int j = 1; j < i; ++j)
					{
						register int jj = (j*(j - 1)) >> 1;
						T Lij = A(G.index[s][i], G.index[s][j]);
						for (register int k = 1; k < j; ++k)
						{
							Lij -= L.entry[ii + k]*L.entry[jj + k];
						}
						Lij /= L.entry[jj + j];
						L.entry[ii + j] = Lij;
						Lii -= Lij*Lij;
					}
					L.entry[ii + i] = sqrt(Lii);
				}
				// Forward sustitution
				for (register int i = 1; i <= N; ++i)
				{
					register int ii = (i*(i - 1)) >> 1;
					register T zi = (i == N) ? 1.0 : 0.0;
					for (register int k = 1; k < i; ++k)
					{
						zi -= L.entry[ii + k]*z.entry[k];
					}
					z.entry[i] = zi/L.entry[ii + i];
				}
				// Backward sustitution
				for (register int i = N; i; --i)
				{
					register int ii = (i*(i - 1)) >> 1;
					register T xi = z.entry[i];
					for (register int k = i + 1; k <= N; ++k)
					{
						register int kk = (k*(k - 1)) >> 1;
						xi -= L.entry[kk + i]*x.entry[k];
					}
					x.entry[i] = xi/L.entry[ii + i];
				}
				for (register int i = 1; i <= N; ++i)
				{
					G.entry[s][i] = x.entry[i];
				}
			}
		}

		// Calculate D and calculate G
		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int i = 1; i <= G.rows; ++i)
		{
			int count = G.Count(i);
			T Dii = 1/sqrt(G.entry[i][count]);
			for (register int c = 1; c <= count; ++c)
			{
				G.entry[i][c] *= Dii;
			}
		}
		Log(1, "G filled");

		Transpose(G, Gt);
		Log(1, "G' filled");
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void SparseApproximateInverseUnsymmetric(const CSRMatrix<T>& A, CSRMatrix<T>& G, CSRMatrix<T>& Gt, T threshold, int preconditioner_level, int threads) throw(Memory::Exception)
{
	// A. Y. Yeremin, A. A. Nikishin
	// Factorized-Sparse-Approximate-Inverse Preconditionings of Linear Systems with Unsymmetric Matrices
	// Journal of Mathematical Sciences, Vol. 121-4. 2004

	try
	{
		Log(1, "SparseApproximateInverseUnsymmetric, threshold = %e, level = %i", threshold, preconditioner_level);

		PowerPreconditionerPattern(A, G, threshold, preconditioner_level, true, threads);

		// Allocate memory for LU factorizations
		Vector<Vector<T> > thread_L(threads);
		Vector<Vector<T> > thread_Ut(threads);
		Vector<Vector<T> > thread_z(threads);
		Vector<Vector<T> > thread_x(threads);
		{
			int N_max = G.Count(1);
			for (register int s = 2; s <= G.rows; ++s)
			{
				if (N_max < G.Count(s))
				{
					N_max = G.Count(s);
				}
			}

			int triangular_size_max = (N_max*(N_max + 1)) >> 1;
			for (register int t = 1; t <= threads; ++t)
			{
				thread_L.entry[t].Resize(triangular_size_max);
				thread_Ut.entry[t].Resize(triangular_size_max);
				thread_z.entry[t].Resize(triangular_size_max);
				thread_x.entry[t].Resize(triangular_size_max);
			}
		}

		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int s = 1; s <= G.rows; ++s)
		{
			// Solve the small definite system, then fill the triangular superior part of M
			int N = G.Count(s);
			if (N == 1)
			{
				register int i = G.index[s][1];
				G.entry[s][1] = 1/A(i, i);
			}
			else if (N == 2)
			{
				register int i1 = G.index[s][1];
				register int i2 = G.index[s][2];
				register T a = A(i1, i1);
				register T b = A(i1, i2);
				register T c = A(i2, i1);
				register T d = A(i2, i2);
				register T det = a*d - b*c;
				G.entry[s][1] = -b/det;
				G.entry[s][2] = a/det;
			}
			else if (N == 3)
			{
				register int i1 = G.index[s][1];
				register int i2 = G.index[s][2];
				register int i3 = G.index[s][3];
				register T a = A(i1, i1);
				register T b = A(i1, i2);
				register T c = A(i1, i3);
				register T d = A(i2, i1);
				register T e = A(i2, i2);
				register T f = A(i2, i3);
				register T g = A(i3, i1);
				register T h = A(i3, i2);
				register T k = A(i3, i3);
				T det = a*(e*k - f*h) - b*(k*d - f*g) + c*(d*h - e*g);
				G.entry[s][1] = (b*f - c*e)/det;
				G.entry[s][2] = (c*d - a*f)/det;
				G.entry[s][3] = (a*e - b*d)/det;
			}
			else
			{
				int thread = omp_get_thread_num() + 1;
				Vector<T>& L = thread_L.entry[thread];
				Vector<T>& Ut = thread_Ut.entry[thread];
				Vector<T>& z = thread_z.entry[thread];
				Vector<T>& x = thread_x.entry[thread];

				// LU factorization
				for (int i = 1; i <= N; ++i)
				{
					int ii = (i*(i - 1)) >> 1;

					register T Utii = A(G.index[s][i], G.index[s][i]);
					for (register int j = 1; j < i; ++j)
					{
						register int jj = (j*(j - 1)) >> 1;

						register T Lij = A(G.index[s][i], G.index[s][j]);
						register T Utij = A(G.index[s][j], G.index[s][i]);
						for (register int k = 1; k < j; ++k)
						{
							Lij -= L.entry[ii + k]*Ut.entry[jj + k];
							Utij -= L.entry[jj + k]*Ut.entry[ii + k];
						}
						Lij /= Ut.entry[jj + j];
						L.entry[ii + j] = Lij;
						Ut.entry[ii + j] = Utij;
						Utii -= Lij*Ut.entry[ii + j];
					}
					Ut.entry[ii + i] = Utii;
					L.entry[ii + i] = 1;
				}

				// Forward sustitution
				for (register int i = 1; i <= N; ++i)
				{
					register int ii = (i*(i - 1)) >> 1;
					register T zi = (i == N) ? 1.0 : 0.0;
					for (register int k = 1; k < i; ++k)
					{
						zi -= L.entry[ii + k]*z.entry[k];
					}
					z.entry[i] = zi/L.entry[ii + i];
				}
				// Backward sustitution
				for (register int i = N; i; --i)
				{
					register int ii = (i*(i - 1)) >> 1;
					register T xi = z.entry[i];
					for (register int k = i + 1; k <= N; ++k)
					{
						register int kk = (k*(k - 1)) >> 1;
						xi -= Ut.entry[kk + i]*x.entry[k];
					}
					x.entry[i] = xi/Ut.entry[ii + i];
				}
				for (register int i = 1; i <= N; ++i)
				{
					G.entry[s][i] = x.entry[i];
				}
			}
		}

		// Calculate D and calculate G
		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int i = 1; i <= G.rows; ++i)
		{
			int count = G.Count(i);
			T Dii = 1/sqrt(G.entry[i][count]);
			for (register int c = 1; c <= count; ++c)
			{
				G.entry[i][c] *= Dii;
			}
		}
		Log(1, "G filled");

		Transpose(G, Gt);
		Log(1, "G' filled");
	}
	catch (Exception&)
	{
		ReThrow();
	}
}

#endif
