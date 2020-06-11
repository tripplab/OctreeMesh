// DenseConjugateGradient.h
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

#ifndef _DenseConjugateGradient_h_
#define _DenseConjugateGradient_h_

#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Container/Matrix.h>
#include <Solver/DenseSolver.h>


template <typename T>
int DenseConjugateGradient(const Matrix<T>& A, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "DenseConjugateGradient:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		int n = x.size;

		Vector<T> g(n); // Gradient
		Vector<T> p(n); // Descent direcction
		Vector<T> w(n); // w = A*p

		T gg = 0.0;
		#pragma omp parallel for default(shared) reduction(+:gg) schedule(guided) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			T* __restrict A_entry_i = A.entry[i];

			T sum = 0;
			for (register int j = 1; j <= A.columns; ++j)
			{
				sum += A_entry_i[j]*x.entry[j];
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
			// Test termination condition
			if (gg <= epsilon) // Norm(Gn) <= tolerance
			{
				break;
			}

			T pw = 0.0;
			#pragma omp parallel for default(shared) reduction(+:pw) schedule(guided) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				T* __restrict A_entry_i = A.entry[i];

				T sum = 0.0;
				for (register int j = 1; j <= A.columns; ++j)
				{
					sum += A_entry_i[j]*p.entry[j];
				}
				w.entry[i] = sum; // w = AP
				pw += p.entry[i]*w.entry[i]; // pw = p'*w
			}

			T alpha = gg/pw; // alpha = (g'*g)/(p'*w)

			T gngn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:gngn) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p.entry[i]; // Xn = x + alpha*p
				g.entry[i] += alpha*w.entry[i]; // Gn = g + alpha*w
				gngn += g.entry[i]*g.entry[i]; // gngn = Gn'*Gn
			}

			T beta = gngn/gg; // beta = (Gn'*Gn)/(g'*g)

			#pragma omp parallel for default(shared) num_threads(threads)
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
			Log(1, "-[Warning] DenseConjugateGradient did not converge in %i steps", max_steps);
		}

		return step;
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
class DenseConjugateGradientSolver : public DenseSolver<T>
{
	public:

		const int steps;


		DenseConjugateGradientSolver(Matrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, int threads) throw(Memory::Exception)
		:	DenseSolver<T>(A, x, b, fixed, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps)
		{
		}


		virtual bool Calculate() throw(Memory::Exception)
		{
			try
			{
				*(int*)&steps = DenseConjugateGradient(this->A, this->x, this->b, tolerance, max_steps, this->threads);
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
int DenseConjugateGradientJacobi(const Matrix<T>& A, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, const Vector<T>& M, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "DenseConjugateGradientJacobi:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		int n = x.size;

		Vector<T> g(n); // Gradient
		Vector<T> p(n); // Descent direcction
		Vector<T> w(n); // w = A*p
		Vector<T> q(n); // Preconditioner system solution

		T gg = 0.0;
		T gq = 0.0;
		#pragma omp parallel for default(shared) reduction(+:gg,gq) schedule(guided) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			T* __restrict A_entry_i = A.entry[i];

			T sum = 0.0;
			for (register int j = 1; j <= A.columns; ++j)
			{
				sum += A_entry_i[j]*x.entry[j];
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
			// Test termination condition
			if (gg <= epsilon) // Norm(Gn) <= tolerance
			{
				break;
			}

			T pw = 0.0;
			#pragma omp parallel for default(shared) reduction(+:pw) schedule(guided) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				T* __restrict A_entry_i = A.entry[i];

				T sum = 0.0;
				for (register int j = 1; j <= A.columns; ++j)
				{
					sum += A_entry_i[j]*p.entry[j];
				}
				w.entry[i] = sum; // w = AP
				pw += p.entry[i]*w.entry[i]; // pw = p'*w
			}

			T alpha = gq/pw; // alpha = (g'*q)/(p'*w)

			T gngn = 0.0;
			T gnqn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:gngn,gnqn) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p.entry[i]; // Xn = x + alpha*p
				g.entry[i] += alpha*w.entry[i]; // Gn = g + alpha*w
				q.entry[i] = g.entry[i]/M.entry[i]; // Solve for q: MQ = g
				gngn += g.entry[i]*g.entry[i]; // gngn = Gn'*Gn
				gnqn += g.entry[i]*q.entry[i]; // gnqn = Gn'*Qn
			}

			T beta = gnqn/gq; // beta = (Gn'*Qn)/(g'*Q)

			#pragma omp parallel for default(shared) num_threads(threads)
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
			Log(1, "-[Warning] DenseConjugateGradientJacobi did not converge in %i steps", max_steps);
		}

		return step;
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
class DenseConjugateGradientJacobiSolver : public DenseSolver<T>
{
	public:

		const int steps;


		DenseConjugateGradientJacobiSolver(Matrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, T tolerance, int max_steps, int threads) throw(Memory::Exception)
		:	DenseSolver<T>(A, x, b, fixed, threads),
			steps(0),
			tolerance(tolerance),
			max_steps(max_steps),
			M(A.rows)
		{
			// Initialize preconditioner
			for (int i = 1; i <= A.rows; ++i)
			{
				M.entry[i] = A.entry[i][i]; // M(i, i) = A(i, i)
			}
		}


		virtual bool Calculate() throw(Memory::Exception)
		{
			try
			{
				*(int*)&steps = DenseConjugateGradientJacobi(this->A, this->x, this->b, tolerance, max_steps, M, this->threads);
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

#endif
