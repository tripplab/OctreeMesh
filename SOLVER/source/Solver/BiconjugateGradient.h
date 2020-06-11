// BiconjugateGradient.h
// Copyright (C) 2012 Miguel Vargas (miguel.vargas@gmail.com)
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

#ifndef _BiconjugateGradient_h_
#define _BiconjugateGradient_h_

#include <Basic/Float.h>
#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Container/CSRMatrix.h>
#include <Container/Vector.h>
#include <Solver/LU.h>
#include <Solver/Solver.h>
#include <Solver/SparseApproximateInverse.h>
#include <Solver/TriangularSystem.h>


template <typename T>
int BiconjugateGradient(const CSRMatrix<T>& A, const CSRMatrix<T>& At, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, int threads) throw(Memory::Exception)
{
	// U. Meier-Yang
	// Preconditioned conjugate gradient-like methods for nonsymmetric linear systems
	// University of Illinois
	// 1992
	// pp. 6-7

	try
	{
		Log(1, "BiconjugateGradient:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		int n = x.size;

		Vector<T> r1(n); // Residual 1
		Vector<T> r2(n); // Residual 2
		Vector<T> p1(n); // Descent direcction 1
		Vector<T> p2(n); // Descent direcction 2
		Vector<T> w1(n); // w1 = A*p
		Vector<T> w2(n); // w2 = A'*p2

		T r1r1 = 0;
		T r2r1 = 0;
		#pragma omp parallel for default(shared) reduction(+:r2r1,r1r1) schedule(guided) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			int* __restrict A_index_i = A.index[i];
			T* __restrict A_entry_i = A.entry[i];

			T sum = 0.0;
			int k_max = A.Count(i);
			for (register int k = 1; k <= k_max; ++k)
			{
				sum += A_entry_i[k]*x.entry[A_index_i[k]];
			}
			r1.entry[i] = b.entry[i] - sum;  // r1 = b - Ax;
			r2.entry[i] = r1.entry[i];       // r2 = r1
			p1.entry[i] = r1.entry[i];       // p1 = r1
			p2.entry[i] = r1.entry[i];       // p2 = r1
			r2r1 += r2.entry[i]*r1.entry[i]; // r2r1 = r2'*r1
			r1r1 += r1.entry[i]*r1.entry[i]; // r1r1 = r1'*r1
		}

		Log(2, "-Step  r'*r");
		T epsilon = tolerance*tolerance;
		int step = 0;
		while (step < max_steps)
		{
			// Test for invalid value
			if (Float<T>::IsNaN(r1r1))
			{
				Log(1, "-[Error] BiconjugateGradientJacobi got NaN. at step %i", step);
				return -1;
			};

			// Test termination condition
			if (r1r1 <= epsilon) // r1'*r1 <= tolerance
			{
				break;
			}

			T p2w1 = 0;
			#pragma omp parallel for default(shared) reduction(+:p2w1) schedule(guided) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict A_index_i = A.index[i];
				T* __restrict A_entry_i = A.entry[i];
				int* __restrict At_index_i = At.index[i];
				T* __restrict At_entry_i = At.entry[i];

				T sum1 = 0.0;
				int k_max1 = A.Count(i);
				for (register int k = 1; k <= k_max1; ++k)
				{
					sum1 += A_entry_i[k]*p1.entry[A_index_i[k]];
				}
				T sum2 = 0.0;
				int k_max2 = At.Count(i);
				for (register int k = 1; k <= k_max2; ++k)
				{
					sum2 += At_entry_i[k]*p2.entry[At_index_i[k]];
				}
				w1.entry[i] = sum1;        // w1 = A*p1
				w2.entry[i] = sum2;        // w2 = A'*p2
				p2w1 += p2.entry[i]*w1.entry[i]; // p2w1 = p2'*w1
			}

			T alpha = r2r1/p2w1; // alpha = (r2'*r1)/(p2'*w1)

			r1r1 = 0;
			T r2nr1n = 0;
			#pragma omp parallel for default(shared) reduction(+:r2nr1n, r1r1) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p1.entry[i];   // Xn = x + alpha*p
				r1.entry[i] -= alpha*w1.entry[i];  // r1n = r1 - alpha*w1
				r2.entry[i] -= alpha*w2.entry[i];  // r2n = r2 - alpha*w2
				r2nr1n += r2.entry[i]*r1.entry[i]; // r2nr1n = r2n'*r1n
				r1r1 += r1.entry[i]*r1.entry[i]; // r1r1 = r1n'*r1n
			}

			T beta = r2nr1n/r2r1; // beta = (r2n'*r1n)/(r2'*r1)

			#pragma omp parallel for default(shared) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				p1.entry[i] = r1.entry[i] + beta*p1.entry[i]; // p1n = r1 + beta*p1
				p2.entry[i] = r2.entry[i] + beta*p2.entry[i]; // p2n = r2 + beta*p2
			}

			r2r1 = r2nr1n;

			Log(2, "%5i  %.5e", step, r1r1);
			++step;
		}
		Log(1, "-Total steps: %i", step);

		if (step >= max_steps)
		{
			Log(1, "-[Error] BiconjugateGradient did not converge in %i steps", max_steps);
			return -1;
		}

		return step;
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
class BiconjugateGradientSolver : public Solver<T>
{
	public:

		const int steps;

		BiconjugateGradientSolver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps),
			At(A.columns, A.rows)
		{
			try
			{
				Transpose(A, At);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

		virtual bool Calculate() throw(Memory::Exception)
		{
			try
			{
				*(int*)&steps = BiconjugateGradient(this->A, At, this->x, this->b, tolerance, max_steps, this->threads);
				if (steps == -1)
				{
					return false;
				}
				return this->CheckSolution();
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

	private:

		T tolerance;
		int max_steps;
		CSRMatrix<T> At;
};


template <typename T>
int BiconjugateGradientJacobi(const CSRMatrix<T>& A, const CSRMatrix<T>& At, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, const Vector<T>& M, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "BiconjugateGradientJacobi:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		// http://en.wikipedia.org/wiki/Biconjugate_gradient_method

		int n = x.size;

		Vector<T> r1(n); // Residual 1
		Vector<T> r2(n); // Residual 2
		Vector<T> p1(n); // Descent direcction 1
		Vector<T> p2(n); // Descent direcction 2
		Vector<T> z1(n); // Preconditioner system solution 1
		Vector<T> z2(n); // Preconditioner system solution 2
		Vector<T> w1(n); // w1 = A*p
		Vector<T> w2(n); // w2 = A'*p2

		T r1r1 = 0;
		T r2z1 = 0;
		#pragma omp parallel for default(shared) reduction(+:r2z1,r1r1) schedule(guided) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			int* __restrict A_index_i = A.index[i];
			T* __restrict A_entry_i = A.entry[i];

			T sum = 0.0;
			int k_max = A.Count(i);
			for (register int k = 1; k <= k_max; ++k)
			{
				sum += A_entry_i[k]*x.entry[A_index_i[k]];
			}
			r1.entry[i] = b.entry[i] - sum;  // r1 = b - Ax;
			r2.entry[i] = r1.entry[i];       // r2 = r1
			z1.entry[i] = r1.entry[i]/M.entry[i];  // Solve for q: M*z1 = r1
			z2.entry[i] = r2.entry[i]/M.entry[i];  // Solve for q: M'*z2 = r2
			p1.entry[i] = z1.entry[i];       // p1 = z1
			p2.entry[i] = z2.entry[i];       // p2 = z2
			r2z1 += r2.entry[i]*z1.entry[i]; // r2z1 = r2'*z1
			r1r1 += r1.entry[i]*r1.entry[i]; // r1r1 = r1'*r1
		}

		Log(2, "-Step  r'*r");
		T epsilon = tolerance*tolerance;
		int step = 0;
		while (step < max_steps)
		{
			// Test for invalid value
			if (Float<T>::IsNaN(r1r1))
			{
				Log(1, "-[Error] BiconjugateGradientJacobi got NaN. at step %i", step);
				return -1;
			};

			// Test termination condition
			if (r1r1 <= epsilon) // r1'*r1 <= tolerance
			{
				break;
			}

			T p2w1 = 0.0;
			#pragma omp parallel for default(shared) reduction(+:p2w1) schedule(guided) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict A_index_i = A.index[i];
				T* __restrict A_entry_i = A.entry[i];
				int* __restrict At_index_i = At.index[i];
				T* __restrict At_entry_i = At.entry[i];

				T sum1 = 0.0;
				int k_max1 = A.Count(i);
				for (register int k = 1; k <= k_max1; ++k)
				{
					sum1 += A_entry_i[k]*p1.entry[A_index_i[k]];
				}
				T sum2 = 0.0;
				int k_max2 = At.Count(i);
				for (register int k = 1; k <= k_max2; ++k)
				{
					sum2 += At_entry_i[k]*p2.entry[At_index_i[k]];
				}
				w1.entry[i] = sum1;        // w1 = A*p1
				w2.entry[i] = sum2;        // w2 = A'*p2
				p2w1 += p2.entry[i]*w1.entry[i]; // p2w1 = p2'*w1
			}

			T alpha = r2z1/p2w1; // alpha = (r2'*r1)/(p2'*w1)

			r1r1 = 0;
			T r2nz1n = 0;
			#pragma omp parallel for default(shared) reduction(+:r2nz1n,r1r1) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p1.entry[i];   // Xn = x + alpha*p
				r1.entry[i] -= alpha*w1.entry[i];  // r1n = r1 - alpha*w1
				r2.entry[i] -= alpha*w2.entry[i];  // r2n = r2 - alpha*w2
				z1.entry[i] = r1.entry[i]/M.entry[i];    // Solve for z1: M*z1 = r1
				z2.entry[i] = r2.entry[i]/M.entry[i];    // Solve for z2: M'*z2 = r2
				r2nz1n += r2.entry[i]*z1.entry[i]; // r2nz1n = r2n'*z1n
				r1r1 += r1.entry[i]*r1.entry[i]; // r1r1 = r1n'*r1n
			}

			T beta = r2nz1n/r2z1; // beta = (r2n'*r1n)/(r2'*r1)

			#pragma omp parallel for default(shared) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				p1.entry[i] = z1.entry[i] + beta*p1.entry[i]; // p1n = r1 + beta*p1
				p2.entry[i] = z2.entry[i] + beta*p2.entry[i]; // p2n = r2 + beta*p2
			}

			r2z1 = r2nz1n;

			Log(2, "%5i  %.5e", step, r1r1);
			++step;
		}
		Log(1, "-Total steps: %i", step);

		if (step >= max_steps)
		{
			Log(1, "-[Error] BiconjugateGradientJacobi did not converge in %i steps", max_steps);
			return -1;
		}

		return step;
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
class BiconjugateGradientJacobiSolver : public Solver<T>
{
	public:

		const int steps;

		BiconjugateGradientJacobiSolver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps),
			At(A.columns, A.rows),
			M(A.columns)
		{
			try
			{
				Transpose(A, At);

				// Initialize preconditioner
				for (int i = 1; i <= A.rows; ++i)
				{
					int* __restrict A_index_i = A.index[i];
					T* __restrict A_entry_i = A.entry[i];

					int k_max = A.Count(i);
					for (register int k = 1; k <= k_max; ++k)
					{
						register int j = A_index_i[k];
						if (i == j)
						{
							M.entry[i] = A_entry_i[k]; // M(i, i) = A(i, i)
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

		virtual bool Calculate() throw(Memory::Exception)
		{
			try
			{
				*(int*)&steps = BiconjugateGradientJacobi(this->A, At, this->x, this->b, tolerance, max_steps, M, this->threads);
				if (steps == -1)
				{
					return false;
				}
				return this->CheckSolution();
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

	private:

		T tolerance;
		int max_steps;
		CSRMatrix<T> At;
		Vector<T> M; // Jacobi preconditioner
};


template <typename T>
int BiconjugateGradientIncompleteLU(const CSRMatrix<T>& A, const CSRMatrix<T>& At, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, const CSRMatrix<T>& L, const CSRMatrix<T>& U, const CSRMatrix<T>& Lt, const CSRMatrix<T>& Ut, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "BiconjugateGradientIncompleteLU:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		// http://en.wikipedia.org/wiki/Biconjugate_gradient_method

		int n = x.size;

		Vector<T> r1(n); // Residual 1
		Vector<T> r2(n); // Residual 2
		Vector<T> p1(n); // Descent direcction 1
		Vector<T> p2(n); // Descent direcction 2
		Vector<T> z1(n); // Preconditioner system solution 1
		Vector<T> z2(n); // Preconditioner system solution 2
		Vector<T> y1(n); // Preconditioner system intermediate solution 1
		Vector<T> y2(n); // Preconditioner system intermediate solution 2
		Vector<T> w1(n); // w1 = A*p
		Vector<T> w2(n); // w2 = A'*p2

		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			int* __restrict A_index_i = A.index[i];
			T* __restrict A_entry_i = A.entry[i];

			T sum = 0.0;
			int k_max = A.Count(i);
			for (register int k = 1; k <= k_max; ++k)
			{
				sum += A_entry_i[k]*x.entry[A_index_i[k]];
			}
			r1.entry[i] = b.entry[i] - sum;  // r1 = b - Ax;
			r2.entry[i] = r1.entry[i];       // r2 = r1
		}

		// Solve for z1: M*z1 = r1
		LowerTriangularSystem(L, y1, r1);
		UpperTriangularSystem(U, z1, y1);

		// Solve for z2: M'*z2 = r2
		LowerTriangularSystem(Ut, y2, r2);
		UpperTriangularSystem(Lt, z2, y2);

		T r1r1 = 0;
		T r2z1 = 0;
		#pragma omp parallel for default(shared) reduction(+:r2z1,r1r1) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			p1.entry[i] = z1.entry[i];       // p1 = z1
			p2.entry[i] = z2.entry[i];       // p2 = z2
			r2z1 += r2.entry[i]*z1.entry[i]; // r2z1 = r2'*z1
			r1r1 += r1.entry[i]*r1.entry[i]; // r1r1 = r1'*r1
		}

		Log(2, "-Step  r'*r");
		T epsilon = tolerance*tolerance;
		int step = 0;
		while (step < max_steps)
		{
			// Test for invalid value
			if (Float<T>::IsNaN(r1r1))
			{
				Log(1, "-[Error] BiconjugateGradientIncompleteLU got NaN.");
				return -1;
			};

			// Test termination condition
			if (r1r1 <= epsilon) // r1'*r1 <= tolerance
			{
				break;
			}

			T p2w1 = 0.0;
			#pragma omp parallel for default(shared) reduction(+:p2w1) schedule(guided) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict A_index_i = A.index[i];
				T* __restrict A_entry_i = A.entry[i];
				int* __restrict At_index_i = At.index[i];
				T* __restrict At_entry_i = At.entry[i];

				T sum1 = 0.0;
				int k_max1 = A.Count(i);
				for (register int k = 1; k <= k_max1; ++k)
				{
					sum1 += A_entry_i[k]*p1.entry[A_index_i[k]];
				}
				T sum2 = 0.0;
				int k_max2 = At.Count(i);
				for (register int k = 1; k <= k_max2; ++k)
				{
					sum2 += At_entry_i[k]*p2.entry[At_index_i[k]];
				}
				w1.entry[i] = sum1;        // w1 = A*p1
				w2.entry[i] = sum2;        // w2 = A'*p2
				p2w1 += p2.entry[i]*w1.entry[i]; // p2w1 = p2'*w1
			}

			T alpha = r2z1/p2w1; // alpha = (r2'*r1)/(p2'*w1)

			#pragma omp parallel for default(shared) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p1.entry[i];   // Xn = x + alpha*p
				r1.entry[i] -= alpha*w1.entry[i];  // r1n = r1 - alpha*w1
				r2.entry[i] -= alpha*w2.entry[i];  // r2n = r2 - alpha*w2
			}

			// Solve for z1: M*z1 = r1
			LowerTriangularSystem(L, y1, r1);
			UpperTriangularSystem(U, z1, y1);

			// Solve for z2: M'*z2 = r2
			LowerTriangularSystem(Ut, y2, r2);
			UpperTriangularSystem(Lt, z2, y2);

			r1r1 = 0;
			T r2nz1n = 0;
			#pragma omp parallel for default(shared) reduction(+:r2nz1n,r1r1) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				r2nz1n += r2.entry[i]*z1.entry[i]; // r2nz1n = r2n'*z1n
				r1r1 += r1.entry[i]*r1.entry[i]; // r1r1 = r1n'*r1n
			}

			T beta = r2nz1n/r2z1; // beta = (r2n'*r1n)/(r2'*r1)

			#pragma omp parallel for default(shared) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				p1.entry[i] = z1.entry[i] + beta*p1.entry[i]; // p1n = r1 + beta*p1
				p2.entry[i] = z2.entry[i] + beta*p2.entry[i]; // p2n = r2 + beta*p2
			}

			r2z1 = r2nz1n;

			Log(2, "%5i  %.5e", step, r1r1);
			++step;
		}
		Log(1, "-Total steps: %i", step);

		if (step >= max_steps)
		{
			Log(1, "-[Error] BiconjugateGradientIncompleteLU did not converge in %i steps", max_steps);
			return -1;
		}

		return step;
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
class BiconjugateGradientIncompleteLUSolver : public Solver<T>
{
	public:

		const int steps;

		BiconjugateGradientIncompleteLUSolver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, int k, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps),
			At(A.columns, A.rows),
			L(A.rows, A.columns),
			U(A.rows, A.columns),
			Lt(A.columns, A.rows),
			Ut(A.columns, A.rows)
		{
			try
			{
				Transpose(A, At);

				// Initialize preconditioners
				SymbolicCholeskyDecomposition(A, L, U, k);
				FillLUDecomposition(A, L, U, this->threads);
				Transpose(L, Lt);
				Transpose(U, Ut);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

		virtual bool Calculate() throw(Memory::Exception)
		{
			try
			{
				*(int*)&steps = BiconjugateGradientIncompleteLU(this->A, At, this->x, this->b, tolerance, max_steps, L, U, Lt, Ut, this->threads);
				if (steps == -1)
				{
					return false;
				}
				return this->CheckSolution();
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

	private:

		T tolerance;
		int max_steps;
		CSRMatrix<T> At;
		CSRMatrix<T> L; // Incomplete LU preconditioner factor L
		CSRMatrix<T> U; // Incomplete LU preconditioner factor U
		CSRMatrix<T> Lt; // Incomplete (LU)' preconditioner factor L'
		CSRMatrix<T> Ut; // Incomplete (LU)' preconditioner factor U'
};


template <typename T>
int BiconjugateGradientSparseApproximateInverse(const CSRMatrix<T>& A, const CSCMatrix<T>& Ac, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, const CSRMatrix<T>& G, const CSRMatrix<T>& Gt, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "BiconjugateGradientSparseApproximateInverse:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		// http://en.wikipedia.org/wiki/Biconjugate_gradient_method

		int n = x.size;

		Vector<T> r1(n); // Residual 1
		Vector<T> r2(n); // Residual 2
		Vector<T> p1(n); // Descent direcction 1
		Vector<T> p2(n); // Descent direcction 2
		Vector<T> Ap1(n); // Ap1 = A*p1
		Vector<T> p2A(n); // p2A = p2*A
		Vector<T> Mr1(n); // Mr1 = M*r1
		Vector<T> r2M(n); // r2M = r2*M
		Vector<T> t1(n); // Middle preconditioner multiplication
		Vector<T> t2(n); // Middle preconditioner multiplication

		#pragma omp parallel for default(shared)schedule(guided) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			int* __restrict A_index_i = A.index[i];
			T* __restrict A_entry_i = A.entry[i];
			T Ax_sum = 0;
			int A_count_i = A.Count(i);
			for (register int k = 1; k <= A_count_i; ++k)
			{
				Ax_sum += A_entry_i[k]*x.entry[A_index_i[k]];
			}
			r1.entry[i] = b.entry[i] - Ax_sum;  // r1 = b - A*x;
			r2.entry[i] = r1.entry[i]; // r2 = r1
		}

		#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			int* __restrict G_index_i = G.index[i];
			T* __restrict G_entry_i = G.entry[i];
			T t1_sum = 0;
			T t2_sum = 0;
			int G_count_i = G.Count(i);
			for (register int k = 1; k <= G_count_i; ++k)
			{
				register int j = G_index_i[k];
				register T v = G_entry_i[k];
				t1_sum += v*r1.entry[j];
				t2_sum += v*r2.entry[j];
			}
			t1.entry[i] = t1_sum; // Mr1 = M*r1
			t2.entry[i] = t2_sum; // r2M = r2*M
		}

		T r2Mr1 = 0;
		T r1r1 = 0;
		#pragma omp parallel for default(shared) reduction(+:r2Mr1,r1r1) schedule(guided) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			int* __restrict Gt_index_i = Gt.index[i];
			T* __restrict Gt_entry_i = Gt.entry[i];
			T Mr1_sum = 0;
			T r2M_sum = 0;
			int Gt_count_i = Gt.Count(i);
			for (register int k = 1; k <= Gt_count_i; ++k)
			{
				register int j = Gt_index_i[k];
				register T v = Gt_entry_i[k];
				Mr1_sum += v*t1.entry[j];
				r2M_sum += v*t2.entry[j];
			}
			Mr1.entry[i] = Mr1_sum; // Mr1 = M*r1
			r2M.entry[i] = r2M_sum; // r2M = r2*M

			p1.entry[i] = Mr1_sum; // p1 = M*r1
			p2.entry[i] = r2M_sum; // p2 = r2*M

			T r1_i = r1.entry[i];
			r2Mr1 += r2M_sum*r1_i; // r2Mr1 = r2*M*r1
			r1r1 += r1_i*r1_i; // r1r1 = r1'*r1
		}

		Log(2, "-Step  r'*r");
		T epsilon = tolerance*tolerance;
		int step = 0;
		while (step < max_steps)
		{
			// Test for invalid value
			if (Float<T>::IsNaN(r1r1))
			{
				Log(1, "-[Error] BiconjugateGradientSparseApproximateInverse got NaN.");
				return -1;
			};

			// Test termination condition
			if (r1r1 <= epsilon) // r1'*r1 <= tolerance
			{
				break;
			}

			T p2Ap1 = 0;
			#pragma omp parallel for default(shared) reduction(+:p2Ap1) schedule(guided) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict A_index_i = A.index[i];
				T* __restrict A_entry_i = A.entry[i];
				T Ap1_sum = 0;
				int A_count_i = A.Count(i);
				for (register int k = 1; k <= A_count_i; ++k)
				{
					Ap1_sum += A_entry_i[k]*p1.entry[A_index_i[k]];
				}
				Ap1.entry[i] = Ap1_sum; // Ap1 = A*p1

				int* __restrict Ac_index_i = Ac.index[i];
				T* __restrict Ac_entry_i = Ac.entry[i];
				T p2A_sum = 0;
				int Ac_count_i = Ac.Count(i);
				for (register int k = 1; k <= Ac_count_i; ++k)
				{
					p2A_sum += p2.entry[Ac_index_i[k]]*Ac_entry_i[k];
				}
				p2A.entry[i] = p2A_sum; // p2A = p2*A

				p2Ap1 += p2.entry[i]*Ap1_sum;
			}

			T alpha = r2Mr1/p2Ap1; // alpha = (r2*M*r1)/(p2*A*p1)

			r1r1 = 0;
			#pragma omp parallel for default(shared) reduction(+:r1r1) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p1.entry[i];   // Xn = x + alpha*p
				r1.entry[i] -= alpha*Ap1.entry[i];  // r1n = r1 - alpha*Ap1
				r2.entry[i] -= alpha*p2A.entry[i];  // r2n = r2 - alpha*p2A

				r1r1 += r1.entry[i]*r1.entry[i]; // r1r1 = r1n'*r1n
			}

			#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict G_index_i = G.index[i];
				T* __restrict G_entry_i = G.entry[i];
				T t1_sum = 0;
				T t2_sum = 0;
				int G_count_i = G.Count(i);
				for (register int k = 1; k <= G_count_i; ++k)
				{
					register int j = G_index_i[k];
					register T v = G_entry_i[k];
					t1_sum += v*r1.entry[j];
					t2_sum += v*r2.entry[j];
				}
				t1.entry[i] = t1_sum; // Mr1 = M*r1
				t2.entry[i] = t2_sum; // r2M = r2*M
			}

			T r2nMr1n = 0;
			#pragma omp parallel for default(shared) reduction(+:r2nMr1n) schedule(guided) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict Gt_index_i = Gt.index[i];
				T* __restrict Gt_entry_i = Gt.entry[i];
				T Mr1_sum = 0;
				T r2M_sum = 0;
				int Gt_count_i = Gt.Count(i);
				for (register int k = 1; k <= Gt_count_i; ++k)
				{
					register int j = Gt_index_i[k];
					register T v = Gt_entry_i[k];
					Mr1_sum += v*t1.entry[j];
					r2M_sum += v*t2.entry[j];
				}
				Mr1.entry[i] = Mr1_sum; // Mr1 = M*r1
				r2M.entry[i] = r2M_sum; // r2M = r2*M

				r2nMr1n += r2M_sum*r1.entry[i];
			}

			T beta = r2nMr1n/r2Mr1; // beta = (r2n*M*r1n)/(r2*M*r1)

			#pragma omp parallel for default(shared) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				p1.entry[i] = Mr1.entry[i] + beta*p1.entry[i]; // p1n = r1 + beta*p1
				p2.entry[i] = r2M.entry[i] + beta*p2.entry[i]; // p2n = r2 + beta*p2
			}

			r2Mr1 = r2nMr1n;

			Log(2, "%5i  %.5e", step, r1r1);
			++step;
		}
		Log(1, "-Total steps: %i", step);

		if (step >= max_steps)
		{
			Log(1, "-[Error] BiconjugateGradientSparseApproximateInverse did not converge in %i steps", max_steps);
			return -1;
		}
		return step;
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
class BiconjugateGradientSparseApproximateInverseSolver : public Solver<T>
{
	public:

		const int steps;


		BiconjugateGradientSparseApproximateInverseSolver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, T preconditioner_threshold, int preconditioner_level, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps),
			Ac(A.rows, A.columns),
			G(A.rows, A.columns),
			Gt(A.columns, A.rows)
		{
			try
			{
				Transpose(A, Ac);

				// Initialize preconditioner
				SparseApproximateInverseUnsymmetric(A, G, Gt, preconditioner_threshold, preconditioner_level, this->threads);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

		virtual bool Calculate() throw(Memory::Exception)
		{
			try
			{
				*(int*)&steps = BiconjugateGradientSparseApproximateInverse(this->A, Ac, this->x, this->b, tolerance, max_steps, G, Gt, this->threads);
				if (steps == -1)
				{
					return false;
				}
				return this->CheckSolution();
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}

	private:

		T tolerance;
		int max_steps;
		CSCMatrix<T> Ac;
		CSRMatrix<T> G; // Sparse approximate inverse preconditioner
		CSRMatrix<T> Gt; // Sparse approximate inverse preconditioner (transpose)
};

#endif
