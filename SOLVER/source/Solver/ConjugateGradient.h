// ConjugateGradient.h
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

#ifndef _ConjugateGradient_h_
#define _ConjugateGradient_h_

#include <Basic/Float.h>
#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Container/CSRMatrix.h>
#include <Container/Vector.h>
#include <Solver/Cholesky.h>
#include <Solver/Solver.h>
#include <Solver/SparseApproximateInverse.h>
#include <Solver/TriangularSystem.h>


template <typename T>
int ConjugateGradient(const CSRMatrix<T>& A, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "ConjugateGradient:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		int n = x.size;

		Vector<T> g(n); // Gradient
		Vector<T> p(n); // Descent direcction
		Vector<T> w(n); // w = A*p

		T gg = 0.0;
		#pragma omp parallel for default(shared) reduction(+:gg) schedule(guided, 64) num_threads(threads)
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
			g.entry[i] = sum - b.entry[i]; // g = AX - b;
			p.entry[i] = -g.entry[i]; // p = -g
			gg += g.entry[i]*g.entry[i]; // gg = g'*g
		}

		Log(2, "-Step  r'*r");
		T epsilon = tolerance*tolerance;
		int step = 0;
		while (step < max_steps)
		{
			// Test for invalid value
			if (Float<T>::IsNaN(gg))
			{
				Log(1, "-[Error] ConjugateGradient got NaN.");
				return -1;
			};

			// Test termination condition
			if (gg <= epsilon) // Norm(Gn) <= tolerance
			{
				break;
			}

			T pw = 0.0;
			#pragma omp parallel for default(shared) reduction(+:pw) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict A_index_i = A.index[i];
				T* __restrict A_entry_i = A.entry[i];

				T sum = 0.0;
				int k_max = A.Count(i);
				for (register int k = 1; k <= k_max; ++k)
				{
					sum += A_entry_i[k]*p.entry[A_index_i[k]];
				}
				w.entry[i] = sum; // w = AP
				pw += p.entry[i]*w.entry[i]; // pw = p'*w
			}

			T alpha = gg/pw; // alpha = (g'*g)/(p'*w)

			T gngn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:gngn) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p.entry[i]; // Xn = x + alpha*p
				g.entry[i] += alpha*w.entry[i]; // Gn = g + alpha*w
				gngn += g.entry[i]*g.entry[i]; // gngn = Gn'*Gn
			}

			T beta = gngn/gg; // beta = (Gn'*Gn)/(g'*g)

			#pragma omp parallel for default(shared) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				p.entry[i] = beta*p.entry[i] - g.entry[i]; // Pn = -g + beta*p
			}

			gg = gngn;

			Log(2, "%5i  %.5e", step, gg);
			++step;
		}
		Log(1, "-Total steps: %i", step);

		if (step >= max_steps)
		{
			Log(1, "-[Error] ConjugateGradient did not converge in %i steps", max_steps);
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
class ConjugateGradientSolver : public Solver<T>
{
	public:

		const int steps;

		ConjugateGradientSolver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps)
		{
		}


		virtual bool Calculate() throw(Memory::Exception)
		{
			try
			{
				*(int*)&steps = ConjugateGradient(this->A, this->x, this->b, tolerance, max_steps, this->threads);
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
};


template <typename T>
int ConjugateGradientJacobi(const CSRMatrix<T>& A, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, const Vector<T>& M, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "ConjugateGradientJacobi:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		int n = x.size;

		Vector<T> g(n); // Gradient
		Vector<T> p(n); // Descent direcction
		Vector<T> w(n); // w = A*p
		Vector<T> q(n); // Preconditioner system solution

		T gg = 0.0;
		T gq = 0.0;
		#pragma omp parallel for default(shared) reduction(+:gg,gq) schedule(guided, 64) num_threads(threads)
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
			g.entry[i] = sum - b.entry[i]; // g = AX - b
			q.entry[i] = g.entry[i]/M.entry[i]; // Solve for q: M*q = g
			p.entry[i] = -q.entry[i]; // p = -q
			gg += g.entry[i]*g.entry[i]; // gg = g'*g
			gq += g.entry[i]*q.entry[i]; // gq = g'*q
		}

		Log(2, "-Step  r'*r");
		T epsilon = tolerance*tolerance;
		int step = 0;
		while (step < max_steps)
		{
			// Test for invalid value
			if (Float<T>::IsNaN(gg))
			{
				Log(1, "-[Error] ConjugateGradientJacobi got NaN.");
				return -1;
			};

			// Test termination condition
			if (gg <= epsilon) // Norm(Gn) <= tolerance
			{
				break;
			}

			T pw = 0.0;
			#pragma omp parallel for default(shared) reduction(+:pw) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict A_index_i = A.index[i];
				T* __restrict A_entry_i = A.entry[i];

				T sum = 0.0;
				int k_max = A.Count(i);
				for (register int k = 1; k <= k_max; ++k)
				{
					sum += A_entry_i[k]*p.entry[A_index_i[k]];
				}
				w.entry[i] = sum; // w = AP
				pw += p.entry[i]*w.entry[i]; // pw = p'*w
			}

			T alpha = gq/pw; // alpha = (g'*q)/(p'*w)

			T gngn = 0.0;
			T gnqn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:gngn,gnqn) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p.entry[i]; // Xn = x + alpha*p
				g.entry[i] += alpha*w.entry[i]; // Gn = g + alpha*w
				q.entry[i] = g.entry[i]/M.entry[i]; // Solve for q: MQ = g
				gngn += g.entry[i]*g.entry[i]; // gngn = Gn'*Gn
				gnqn += g.entry[i]*q.entry[i]; // gnqn = Gn'*Qn
			}

			T beta = gnqn/gq; // beta = (Gn'*Qn)/(g'*Q)

			#pragma omp parallel for default(shared) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				p.entry[i] = beta*p.entry[i] - q.entry[i]; // Pn = -q + beta*p
			}

			gg = gngn;
			gq = gnqn;

			Log(2, "%5i  %.5e", step, gg);
			++step;
		}
		Log(1, "-Total steps: %i", step);

		if (step >= max_steps)
		{
			Log(1, "-[Error] ConjugateGradientJacobi did not converge in %i steps", max_steps);
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
class ConjugateGradientJacobiSolver : public Solver<T>
{
	public:

		const int steps;

		ConjugateGradientJacobiSolver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps),
			M(A.rows)
		{
			// Initialize preconditioner
			for (int i = 1; i <= A.rows; ++i)
			{
				M.entry[i] = A(i, i); // M(i, i) = A(i, i)
			}
		}


		virtual bool Calculate() throw(Memory::Exception)
		{
			try
			{
				*(int*)&steps = ConjugateGradientJacobi(this->A, this->x, this->b, tolerance, max_steps, M, this->threads);
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
		Vector<T> M; // Jacobi preconditioner
};


template <typename T>
int ConjugateGradientIncompleteCholesky(const CSRMatrix<T>& A, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, const CSRMatrix<T>& L, const CSRMatrix<T>& Lt, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "ConjugateGradientIncompleteCholesky:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		int n = x.size;

		Vector<T> g(n); // Gradient
		Vector<T> p(n); // Descent direcction
		Vector<T> w(n); // w = A*p
		Vector<T> q(n); // Preconditioner system solution
		Vector<T> r(n); // Preconditioner system intermediate solution

		T gg = 0.0;
		#pragma omp parallel for default(shared) reduction(+:gg) schedule(guided, 64) num_threads(threads)
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
			g.entry[i] = sum - b.entry[i]; // g = AX - b
			gg += g.entry[i]*g.entry[i]; // gg = g'*g
		}

		// Solve for q: MQ = g
		LowerTriangularSystem(L, r, g);
		UpperTriangularSystem(Lt, q, r);

		T gq = 0.0;
		#pragma omp parallel for default(shared) reduction(+:gq) schedule(guided, 64) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			p.entry[i] = -q.entry[i]; // p = -q
			gq += g.entry[i]*q.entry[i]; // gq = g'*q
		}

		Log(2, "-Step  r'*r");
		T epsilon = tolerance*tolerance;
		int step = 0;
		while (step < max_steps)
		{
			// Test for invalid value
			if (Float<T>::IsNaN(gg))
			{
				Log(1, "-[Error] ConjugateGradientIncompleteCholesky got NaN.");
				return -1;
			};

			// Test termination condition
			if (gg <= epsilon) // Norm(Gn) <= tolerance
			{
				break;
			}

			T pw = 0.0;
			#pragma omp parallel for default(shared) reduction(+:pw) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict A_index_i = A.index[i];
				T* __restrict A_entry_i = A.entry[i];

				T sum = 0.0;
				int k_max = A.Count(i);
				for (register int k = 1; k <= k_max; ++k)
				{
					sum += A_entry_i[k]*p.entry[A_index_i[k]];
				}
				w.entry[i] = sum; // w = AP
				pw += p.entry[i]*w.entry[i]; // pw = p'*w
			}

			T alpha = gq/pw; // alpha = (g'*q)/(p'*w)

			T gngn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:gngn) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p.entry[i]; // Xn = x + alpha*p
				g.entry[i] += alpha*w.entry[i]; // Gn = g + alpha*w
				gngn += g.entry[i]*g.entry[i]; // gngn = Gn'*Gn
			}

			// Solve for q: MQ = g
			LowerTriangularSystem(L, r, g);
			UpperTriangularSystem(Lt, q, r);

			T gnqn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:gnqn) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				gnqn += g.entry[i]*q.entry[i]; // gnqn = Gn'*Qn
			}

			T beta = gnqn/gq; // beta = (Gn'*Gn)/(g'*g)

			#pragma omp parallel for default(shared) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				p.entry[i] = beta*p.entry[i] - q.entry[i]; // Pn = -q + beta*p
			}

			gg = gngn;
			gq = gnqn;

			Log(2, "%5i  %.5e", step, gg);
			++step;
		}
		Log(1, "-Total steps: %i", step);

		if (step >= max_steps)
		{
			Log(1, "-[Error] ConjugateGradientIncompleteCholesky did not converge in %i steps", max_steps);
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
class ConjugateGradientIncompleteCholeskySolver : public Solver<T>
{
	public:

		const int steps;


		ConjugateGradientIncompleteCholeskySolver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, int preconditioner_level, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps),
			L(A.rows, A.columns),
			Lt(A.rows, A.columns)
		{
			try
			{
				// Initialize preconditioner
				SymbolicCholeskyDecomposition(A, L, Lt, preconditioner_level);
				FillCholeskyDecomposition(A, L, Lt, this->threads);
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
				*(int*)&steps = ConjugateGradientIncompleteCholesky(this->A, this->x, this->b, tolerance, max_steps, L, Lt, this->threads);
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
		CSRMatrix<T> L;  // Incomplete Cholesky preconditioner factor L
		CSRMatrix<T> Lt; // Incomplete Cholesky preconditioner factor L'
};


template <typename T>
int ConjugateGradientIncompleteCholesky2(const CSRMatrix<T>& A, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, const CSRMatrix<T>& L, const Vector<T>& D, const CSRMatrix<T>& Lt, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "ConjugateGradientIncompleteCholesky2:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		int n = x.size;

		Vector<T> g(n); // Gradient
		Vector<T> p(n); // Descent direcction
		Vector<T> w(n); // w = A*p
		Vector<T> q(n); // Preconditioner system solution
		Vector<T> r(n); // Preconditioner system intermediate solution
		Vector<T> s(n); // Preconditioner system intermediate solution

		T gg = 0.0;
		#pragma omp parallel for default(shared) reduction(+:gg) schedule(guided, 64) num_threads(threads)
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
			g.entry[i] = sum - b.entry[i]; // g = AX - b
			gg += g.entry[i]*g.entry[i]; // gg = g'*g
		}

		// Solve for q: LDL'q = g
		LowerTriangularSystem(L, r, g);
		DiagonalSystem(D, s, r);
		UpperTriangularSystem(Lt, q, s);

		T gq = 0.0;
		#pragma omp parallel for default(shared) reduction(+:gq) schedule(guided, 64) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			p.entry[i] = -q.entry[i]; // p = -q
			gq += g.entry[i]*q.entry[i]; // gq = g'*q
		}

		Log(2, "-Step  r'*r");
		T epsilon = tolerance*tolerance;
		int step = 0;
		while (step < max_steps)
		{
			// Test for invalid value
			if (Float<T>::IsNaN(gg))
			{
				Log(1, "-[Error] ConjugateGradientIncompleteCholesky2 got NaN.");
				return -1;
			};

			// Test termination condition
			if (gg <= epsilon) // Norm(Gn) <= tolerance
			{
				break;
			}

			T pw = 0.0;
			#pragma omp parallel for default(shared) reduction(+:pw) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict A_index_i = A.index[i];
				T* __restrict A_entry_i = A.entry[i];

				T sum = 0.0;
				int k_max = A.Count(i);
				for (register int k = 1; k <= k_max; ++k)
				{
					sum += A_entry_i[k]*p.entry[A_index_i[k]];
				}
				w.entry[i] = sum; // w = AP
				pw += p.entry[i]*w.entry[i]; // pw = p'*w
			}

			T alpha = gq/pw; // alpha = (g'*q)/(p'*w)

			T gngn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:gngn) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p.entry[i]; // Xn = x + alpha*p
				g.entry[i] += alpha*w.entry[i]; // Gn = g + alpha*w
				gngn += g.entry[i]*g.entry[i]; // gngn = Gn'*Gn
			}

			// Solve for q: LDL'q = g
			LowerTriangularSystem(L, r, g);
			DiagonalSystem(D, s, r);
			UpperTriangularSystem(Lt, q, s);

			T gnqn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:gnqn) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				gnqn += g.entry[i]*q.entry[i]; // gnqn = Gn'*Qn
			}

			T beta = gnqn/gq; // beta = (Gn'*Gn)/(g'*g)

			#pragma omp parallel for default(shared) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				p.entry[i] = beta*p.entry[i] - q.entry[i]; // Pn = -q + beta*p
			}

			gg = gngn;
			gq = gnqn;

			Log(2, "%5i  %.5e", step, gg);
			++step;
		}
		Log(1, "-Total steps: %i", step);

		if (step >= max_steps)
		{
			Log(1, "-[Error] ConjugateGradientIncompleteCholesky2 did not converge in %i steps", max_steps);
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
class ConjugateGradientIncompleteCholesky2Solver : public Solver<T>
{
	public:

		const int steps;

		ConjugateGradientIncompleteCholesky2Solver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, T preconditioner_threshold, int preconditioner_level, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps),
			L(A.rows, A.columns),
			D(A.rows),
			Lt(A.rows, A.columns)
		{
			try
			{
				// Initialize preconditioner
				SymbolicCholeskyDecomposition(A, L, Lt, preconditioner_level);
				FillCholesky2Decomposition(A, L, D, Lt, this->threads, preconditioner_threshold, 0.01);

				// Check for a valid preconditioner
				for (int i = 1; i < D.size; ++i)
				{
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
				*(int*)&steps = ConjugateGradientIncompleteCholesky2(this->A, this->x, this->b, tolerance, max_steps, L, D, Lt, this->threads);
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
		CSRMatrix<T> L;  // Incomplete Cholesky preconditioner factor L
		Vector<T> D;     // Incomplete Cholesky preconditioner factor D
		CSRMatrix<T> Lt; // Incomplete Cholesky preconditioner factor L'
};


template <typename T>
int ConjugateGradientSparseApproximateInverse(const CSRMatrix<T>& A, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, const CSRMatrix<T>& G, const CSRMatrix<T>& Gt, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "ConjugateGradientSparseApproximateInverse:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		int n = x.size;

		Vector<T> r(n); // Gradient
		Vector<T> p(n); // Descent direcction
		Vector<T> Ap(n); // Ap = A*p
		Vector<T> Mr(n); // Preconditioner system solution
		Vector<T> t(n); // Middle preconditioner multiplication

		T rr = 0;
		#pragma omp parallel for default(shared) reduction(+:rr) schedule(guided, 64) num_threads(threads)
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
			r.entry[i] = b.entry[i] - Ax_sum; // r = b - A*x

			rr += r.entry[i]*r.entry[i]; // rr = r'*r
		}

		#pragma omp parallel for default(shared) schedule(guided, 64) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			int* __restrict G_index_i = G.index[i];
			T* __restrict G_entry_i = G.entry[i];
			T ti_sum = 0;
			int G_count_i = G.Count(i);
			for (register int k = 1; k <= G_count_i; ++k)
			{
				ti_sum += G_entry_i[k]*r.entry[G_index_i[k]];
			}
			t.entry[i] = ti_sum; // Mr = M*r
		}

		T rMr = 0;
		#pragma omp parallel for default(shared) reduction(+:rMr) schedule(guided, 64) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			int* __restrict Gt_index_i = Gt.index[i];
			T* __restrict Gt_entry_i = Gt.entry[i];
			T Mr_sum = 0;
			int Gt_count_i = Gt.Count(i);
			for (register int k = 1; k <= Gt_count_i; ++k)
			{
				Mr_sum += Gt_entry_i[k]*t.entry[Gt_index_i[k]];
			}
			Mr.entry[i] = Mr_sum; // Mr = M*r

			p.entry[i] = Mr_sum;
			rMr += r.entry[i]*Mr.entry[i]; // rMr = r'*M*r
		}

		Log(2, "-Step  r'*r");
		T epsilon = tolerance*tolerance;
		int step = 0;
		while (step < max_steps)
		{
			// Test for invalid value
			if (Float<T>::IsNaN(rr))
			{
				Log(1, "-[Error] ConjugateGradientSparseApproximateInverse got NaN.");
				return -1;
			};

			// Test termination condition
			if (rr <= epsilon) // Norm(Gn) <= tolerance
			{
				break;
			}

			T pAp = 0.0;
			#pragma omp parallel for default(shared) reduction(+:pAp) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict A_index_i = A.index[i];
				T* __restrict A_entry_i = A.entry[i];
				T Ap_sum = 0.0;
				int A_count_i = A.Count(i);
				for (register int k = 1; k <= A_count_i; ++k)
				{
					Ap_sum += A_entry_i[k]*p.entry[A_index_i[k]];
				}
				Ap.entry[i] = Ap_sum; // Ap = A*p

				pAp += p.entry[i]*Ap.entry[i]; // pAp = p'*A*p
			}

			T alpha = rMr/pAp; // alpha = (r'*M*r)/(p'*A*p)

			T rnrn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:rnrn) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p.entry[i]; // xn = x + alpha*p
				r.entry[i] -= alpha*Ap.entry[i]; // rn = r - alpha*A*p

				rnrn += r.entry[i]*r.entry[i]; // rnrn = rn'*rn
			}

			#pragma omp parallel for default(shared) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict G_index_i = G.index[i];
				T* __restrict G_entry_i = G.entry[i];
				T ti_sum = 0;
				int G_count_i = G.Count(i);
				for (register int k = 1; k <= G_count_i; ++k)
				{
					ti_sum += G_entry_i[k]*r.entry[G_index_i[k]];
				}
				t.entry[i] = ti_sum; // Mr = M*r
			}

			T rnMrn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:rnMrn) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict Gt_index_i = Gt.index[i];
				T* __restrict Gt_entry_i = Gt.entry[i];
				T Mr_sum = 0;
				int Gt_count_i = Gt.Count(i);
				for (register int k = 1; k <= Gt_count_i; ++k)
				{
					Mr_sum += Gt_entry_i[k]*t.entry[Gt_index_i[k]];
				}
				Mr.entry[i] = Mr_sum; // Mr = M*r

				rnMrn += r.entry[i]*Mr.entry[i]; // rnMrn = rn'*M*rn
			}

			T beta = rnMrn/rMr; // beta = (rn'*M*rn)/(r'*M*r)

			#pragma omp parallel for default(shared) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				p.entry[i] = Mr.entry[i] + beta*p.entry[i]; // Pn = M*r + beta*p
			}

			rr = rnrn;
			rMr = rnMrn;

			Log(2, "%5i  %.5e", step, rr);
			++step;
		}
		Log(1, "-Total steps: %i", step);

		if (step >= max_steps)
		{
			Log(1, "-[Error] ConjugateGradientSparseApproximateInverse did not converge in %i steps", max_steps);
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
class ConjugateGradientSparseApproximateInverseSolver : public Solver<T>
{
	public:

		const int steps;


		ConjugateGradientSparseApproximateInverseSolver(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, T preconditioner_threshold, int preconditioner_level, bool matrix_adjust, int threads) throw(Memory::Exception)
		:	Solver<T>(A, x, b, fixed, matrix_adjust, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps),
			G(A.rows, A.columns),
			Gt(A.columns, A.rows)
		{
			try
			{
				// Initialize preconditioner
				SparseApproximateInverseSymmetric(A, G, Gt, preconditioner_threshold, preconditioner_level, this->threads);
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
				*(int*)&steps = ConjugateGradientSparseApproximateInverse(this->A, this->x, this->b, tolerance, max_steps, G, Gt, this->threads);
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
		CSRMatrix<T> G; // Sparse approximate inverse preconditioner
		CSRMatrix<T> Gt; // Sparse approximate inverse preconditioner (transpose)
};

#endif
