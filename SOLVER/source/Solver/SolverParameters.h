// SolverParameters.h
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

#ifndef _SolverParameters_h_
#define _SolverParameters_h_

#include <File/File.h>
#include <Solver/BiconjugateGradient.h>
#include <Solver/Cholesky.h>
#include <Solver/ConjugateGradient.h>
#include <Solver/DenseSolver.h>
#include <Solver/DenseCholesky.h>
#include <Solver/DenseConjugateGradient.h>
#include <Solver/LU.h>
#include <Solver/Solver.h>

#include <string.h>
#include <omp.h>


#define SOLVER_SECTION_NAME_MAX_SIZE 256


enum SolverType
{
	solver_conjugate_gradient = 1,
	solver_cholesky_decomposition = 2,
	solver_cholesky2_decomposition = 3,
	solver_biconjugate_gradient = 4,
	solver_lu_decomposition = 5
};


enum PreconditionerType
{
	preconditioner_none = 0,
	preconditioner_jacobi = 1,
	preconditioner_incomplete_cholesky = 2,
	preconditioner_incomplete_cholesky2 = 3,
	preconditioner_incomplete_lu = 4,
	preconditioner_sparse_approximate_inverse = 5
};


template <typename T>
class SolverParameters
{
	public:

		// Solver
		SolverType type;
		int threads;
		T tolerance;
		int max_steps;
		PreconditionerType preconditioner;
		int preconditioner_level;
		T preconditioner_threshold;

		// Reorder to improve Cholesky and LU factorizations
		bool reorder;

		// Substructuring
		int substructuring_threads;
		T substructuring_tolerance;
		bool multiple_results;


		SolverParameters(SolverType type, int threads, T tolerance, int max_steps, PreconditionerType preconditioner, int preconditioner_level, T preconditioner_threshold)
		:	type(type),
			threads(threads),
			tolerance(tolerance),
			max_steps(max_steps),
			preconditioner(preconditioner),
			preconditioner_level(preconditioner_level),
			preconditioner_threshold(preconditioner_threshold),
			reorder((type == solver_cholesky_decomposition) || (type == solver_cholesky2_decomposition) || (type == solver_lu_decomposition)),
			substructuring_threads(),
			substructuring_tolerance(),
			multiple_results()
		{
		}


		SolverParameters(SolverType type, int threads, T tolerance, int max_steps, PreconditionerType preconditioner, int preconditioner_level, T preconditioner_threshold, int substructuring_threads, T substructuring_tolerance, bool multiple_results)
		:	type(type),
			threads(threads),
			tolerance(tolerance),
			max_steps(max_steps),
			preconditioner(preconditioner),
			preconditioner_level(preconditioner_level),
			preconditioner_threshold(preconditioner_threshold),
			reorder((type == solver_cholesky_decomposition) || (type == solver_cholesky2_decomposition) || (type == solver_lu_decomposition)),
			substructuring_threads(substructuring_threads),
			substructuring_tolerance(substructuring_tolerance),
			multiple_results(multiple_results)
		{
		}


		SolverParameters(const char* file_name) throw(File::ExceptionOpen, File::ExceptionEOF, File::ExceptionRead, File::ExceptionFormat)
		{
			try
			{
				char section_name[SOLVER_SECTION_NAME_MAX_SIZE];
				int tmp_integer;

				// Load Solver
				File file;
				file.Open(file_name);
				file.SkipComments(';');

				// Load solver parameters
				file.Read(section_name, SOLVER_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Solver}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(tmp_integer);
				type = (SolverType)tmp_integer;
				file.SkipComments(';');
				file.Read(threads);
				file.SkipComments(';');
				file.Read(tolerance);
				file.SkipComments(';');
				file.Read(max_steps);
				file.SkipComments(';');
				file.Read(tmp_integer);
				preconditioner = (PreconditionerType)tmp_integer;
				file.SkipComments(';');
				file.Read(preconditioner_level);
				file.SkipComments(';');
				file.Read(preconditioner_threshold);
				file.SkipComments(';');

				// Define if reordering is needed
				reorder =
					(type == solver_cholesky_decomposition) ||
					(type == solver_cholesky2_decomposition) ||
					(type == solver_lu_decomposition);

				// Load domain decomposition parameters
				file.Read(section_name, SOLVER_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{Substructuring}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(substructuring_threads);
				file.SkipComments(';');
				file.Read(substructuring_tolerance);
				file.SkipComments(';');
				file.Read(multiple_results);
				file.SkipComments(';');

				file.Close();
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void PrintInfo() throw()
		{
			static const char* solver_name[] =
			{
				"(Undefined)",
				"Conjugate gradient",
				"Cholesky decomposition",
				"Cholesky2 decomposition",
				"Biconjugate gradient",
				"LU decomposition"
			};

			static const char* preconditioner_name[] =
			{
				"None",
				"Jacobi",
				"Incomplete cholesky",
				"Incomplete cholesky2",
				"Incomplete lu",
				"Sparse approximate inverse"
			};

			Log(1, "Solver ---------------------------------------------------------------");
			Log(1, "-Type:              %s", solver_name[type]);
			Log(1, "-Threads:           %i", threads);
			Log(1, "-Reorder equations: %s", reorder ? "yes" : "no");
			if ((type == solver_conjugate_gradient) || (type == solver_biconjugate_gradient))
			{
				Log(1, "-Tolerance:         %.5g", tolerance);
				Log(1, "-Maximum steps:     %i", max_steps);
				Log(1, "-Preconditioner:    %s", preconditioner_name[preconditioner]);
				if ((preconditioner != preconditioner_none) || (preconditioner != preconditioner_jacobi))
				{
					Log(1, "-Level:             %i", preconditioner_level);
					if (preconditioner == preconditioner_sparse_approximate_inverse)
					{
						Log(1, "-Threshold:         %.5g", preconditioner_threshold);
					}
				}
			}
		}


		Solver<T>* Instanciate(CSRMatrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed, bool matrix_adjust = false) const throw(Memory::Exception)
		{
			switch (type)
			{
				case solver_conjugate_gradient:
				{
					switch (preconditioner)
					{
						case preconditioner_none:
						{
							return new ConjugateGradientSolver<T>(A, x, b, fixed, tolerance, max_steps, matrix_adjust, threads);
						}
						case preconditioner_jacobi:
						{
							return new ConjugateGradientJacobiSolver<T>(A, x, b, fixed, tolerance, max_steps, matrix_adjust, threads);
						}
						case preconditioner_incomplete_cholesky:
						{
							return new ConjugateGradientIncompleteCholeskySolver<T>(A, x, b, fixed, tolerance, max_steps, preconditioner_level, matrix_adjust, threads);
						}
						case preconditioner_incomplete_cholesky2:
						{
							return new ConjugateGradientIncompleteCholesky2Solver<T>(A, x, b, fixed, tolerance, max_steps, preconditioner_threshold, preconditioner_level, matrix_adjust, threads);
						}
						case preconditioner_incomplete_lu:
						{
							return (Solver<T>*)0;
						}
						case preconditioner_sparse_approximate_inverse:
						{
							return new ConjugateGradientSparseApproximateInverseSolver<T>(A, x, b, fixed, tolerance, max_steps, preconditioner_threshold, preconditioner_level, matrix_adjust, threads);
						}
					}
				}
				case solver_cholesky_decomposition:
				{
					return new CholeskySolver<T>(A, x, b, fixed, matrix_adjust, threads);
				}
				case solver_cholesky2_decomposition:
				{
					return new Cholesky2Solver<T>(A, x, b, fixed, matrix_adjust, threads);
				}
				case solver_biconjugate_gradient:
				{
					switch (preconditioner)
					{
						case preconditioner_none:
						{
							return new BiconjugateGradientSolver<T>(A, x, b, fixed, tolerance, max_steps, matrix_adjust, threads);
						}
						case preconditioner_jacobi:
						{
							return new BiconjugateGradientJacobiSolver<T>(A, x, b, fixed, tolerance, max_steps, matrix_adjust, threads);
						}
						case preconditioner_incomplete_cholesky:
						{
							return (Solver<T>*)0;
						}
						case preconditioner_incomplete_cholesky2:
						{
							return (Solver<T>*)0;
						}
						case preconditioner_incomplete_lu:
						{
							return new BiconjugateGradientIncompleteLUSolver<T>(A, x, b, fixed, tolerance, max_steps, preconditioner_level, matrix_adjust, threads);
						}
						case preconditioner_sparse_approximate_inverse:
						{
							return new BiconjugateGradientSparseApproximateInverseSolver<T>(A, x, b, fixed, tolerance, max_steps, preconditioner_threshold, preconditioner_level, matrix_adjust, threads);
						}
					}
				}
				case solver_lu_decomposition:
				{
					return new LUSolver<T>(A, x, b, fixed, matrix_adjust, threads);
				}
			}

			return (Solver<T>*)0;
		}


		DenseSolver<T>* Instanciate(Matrix<T>& A, Vector<T>& x, Vector<T>& b, Vector<bool>& fixed) const throw(Memory::Exception)
		{
			switch (type)
			{
				case solver_conjugate_gradient:
				{
					switch (preconditioner)
					{
						case preconditioner_none:
						{
							return new DenseConjugateGradientSolver<T>(A, x, b, fixed, tolerance, max_steps, threads);
						}
						case preconditioner_jacobi:
						{
							return new DenseConjugateGradientJacobiSolver<T>(A, x, b, fixed, tolerance, max_steps, threads);
						}
						case preconditioner_incomplete_cholesky:
						{
							return (DenseSolver<T>*)0;
						}
						case preconditioner_incomplete_cholesky2:
						{
							return (DenseSolver<T>*)0;
						}
						case preconditioner_incomplete_lu:
						{
							return (DenseSolver<T>*)0;
						}
						case preconditioner_sparse_approximate_inverse:
						{
							return (DenseSolver<T>*)0;
						}
					}
				}
				case solver_cholesky_decomposition:
				{
					return new DenseCholeskySolver<T>(A, x, b, fixed, threads);
				}
				case solver_cholesky2_decomposition:
				{
					return (DenseSolver<T>*)0;
				}
				case solver_biconjugate_gradient:
				{
					return (DenseSolver<T>*)0;
				}
				case solver_lu_decomposition:
				{
					return (DenseSolver<T>*)0;
				}
			}

			return (DenseSolver<T>*)0;
		}
};

#endif
