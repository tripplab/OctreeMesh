// main.cpp
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

#include <Basic/Debug.h>
#include <Basic/Log.h>
#include <Basic/Macros.h>
#include <File/MatFile.h>
#include <Container/CSRMatrix.h>
#include <Container/Vector.h>
#include <Solver/BiconjugateGradient.h>
#include <Solver/Cholesky.h>
#include <Solver/ConjugateGradient.h>
#include <Solver/LU.h>
#include <Solver/TriangularSystem.h>
#include <omp.h>


#define LOG_LEVEL                2
#define SOLVER_TYPE              1 // 1=Conjugate_gradient, 2=Cholesky_decomposition, 3=Cholesky2_decomposition, 4=Biconjugate_gradient, 5=LU_decomposition
#define SOLVER_THREADS           1
#define SOLVER_TOLERANCE         1e-5
#define SOLVER_MAX_STEPS         10000
#define SOLVER_AUTO_REORDER      1
#define PRECONDITIONER_TYPE      1 // 0=None, 1=Jacobi, 2=Incomplete_Cholesky, 3=Incomplete_Cholesky2, 4=Incomplete_LU, 5=Sparse_Approximate_Inverse
#define PRECONDITIONER_LEVEL     1
#define PRECONDITIONER_THRESHOLD 0.0


int main(int argc, char** argv)
{
	if (argc == 1)
	{
		puts("Linear solvers A*x=b for sparse matrices with symetric structure.\n");
		puts("Only MAT-File 4 files are supported.\n");
		puts("Input MAT file should contain a sparse matrix 'A' and a vector 'b':\n");
		puts("Output MAT file will contain a vector 'x':\n");
	}

	if (argc != 3)
	{
		fprintf(stderr, "Invalid number of arguments. Use:\n  %s <input mat file> <output mat file>\n", argv[0]);
		return 1;
	}

	try
	{
		// Load environment parameters
		log_level = GetEnvInteger(LOG_LEVEL);
		int solver_type = GetEnvInteger(SOLVER_TYPE);
		int solver_threads = GetEnvInteger(SOLVER_THREADS);
		double solver_tolerance = GetEnvFloat(SOLVER_TOLERANCE);
		int solver_max_steps = GetEnvInteger(SOLVER_MAX_STEPS);
		int solver_auto_reorder = GetEnvInteger(SOLVER_AUTO_REORDER);
		int preconditioner_type = GetEnvInteger(PRECONDITIONER_TYPE);
		int preconditioner_level = GetEnvInteger(PRECONDITIONER_LEVEL);
		double preconditioner_threshold = GetEnvFloat(PRECONDITIONER_THRESHOLD);

		Log(1, "------------------------------ MatSolver -----------------------------");
		Log(1, "Version:   %s", MacroValueToString(VERSION));

		// Load system of equations
		CSRMatrix<double> A;
		Vector<double> b;
		MatFile mat_file;
		mat_file.Open(argv[1]);
		mat_file.Retrive("A", A);
		mat_file.Retrive("b", b);
		mat_file.Close();
		Log(1, "Equations: %i", A.rows);
		Log(1, "nnz(A):    %i", A.NonZero());

		// Set number of threads for matrix operations
		omp_set_num_threads(solver_threads);

		Log(1, "Threads:   %i", solver_threads);

		// Call solver
		Vector<double> x(A.rows);
		switch (solver_type)
		{
			case 1: // Conjugate_gradient
			{
				x.Fill(0);
				switch (preconditioner_type)
				{
					case 0: // None
					{
						ConjugateGradient(A, x, b, solver_tolerance, solver_max_steps, solver_threads);
						break;
					}
					case 1: // Jacobi
					{
						Vector<double> M(A.rows);
						for (int i = 1; i <= A.rows; ++i)
						{
							M.entry[i] = A(i, i);
						}
						ConjugateGradientJacobi(A, x, b, solver_tolerance, solver_max_steps, M, solver_threads);
						break;
					}
					case 2: // Incomplete_Cholesky
					{
						CSRMatrix<double> L(A.rows, A.columns);
						CSRMatrix<double> Lt(A.columns, A.rows);
						SymbolicCholeskyDecomposition(A, L, Lt, preconditioner_level);
						FillCholeskyDecomposition(A, L, Lt, solver_threads);
						ConjugateGradientIncompleteCholesky(A, x, b, solver_tolerance, solver_max_steps, L, Lt, solver_threads);
						break;
					}
					case 3: // Incomplete_Cholesky2
					{
						CSRMatrix<double> L(A.rows, A.columns);
						Vector<double> D(A.rows);
						CSRMatrix<double> Lt(A.columns, A.rows);
						SymbolicCholeskyDecomposition(A, L, Lt, preconditioner_level);
						FillCholesky2Decomposition(A, L, D, Lt, solver_threads, preconditioner_threshold, 0.01);
						ConjugateGradientIncompleteCholesky2(A, x, b, solver_tolerance, solver_max_steps, L, D, Lt, solver_threads);
						break;
					}
					case 5: // Sparse_Approximate_Inverse
					{
						CSRMatrix<double> G(A.rows, A.columns);
						CSRMatrix<double> Gt(A.columns, A.rows);
						SparseApproximateInverseSymmetric(A, G, Gt, preconditioner_threshold, preconditioner_level, solver_threads);
						ConjugateGradientSparseApproximateInverse(A, x, b, solver_tolerance, solver_max_steps, G, Gt, solver_threads);
						break;
					}
					default:
					{
						Log(0, "Error: Invalid preconditioner %i", preconditioner_type);
						return 1;
					}
				}
				break;
			}
			case 2: // Cholesky_decomposition
			{
				if (solver_auto_reorder)
				{
					CSRMatrix<double> Ar(A.rows, A.columns);
					Vector<double> xr(x.size);
					Vector<double> br(b.size);
					Vector<int> index;
					Vector<int> inverse_index;
					FindingAnOrdering(A, index, inverse_index);
					Reorder(index, inverse_index, A, Ar, solver_threads);
					Reorder(index, b, br, solver_threads);
					CSRMatrix<double> L(Ar.rows, Ar.columns);
					CSRMatrix<double> Lt(Ar.columns, Ar.rows);
					Vector<double> cr(Ar.rows);
					SymbolicCholeskyDecomposition(Ar, L, Lt, -1);
					FillCholeskyDecomposition(Ar, L, Lt, solver_threads);
					LowerTriangularSystem(L, cr, br);
					UpperTriangularSystem(Lt, xr, cr);
					Reorder(inverse_index, xr, x, solver_threads);
				}
				else
				{
					CSRMatrix<double> L(A.rows, A.columns);
					CSRMatrix<double> Lt(A.columns, A.rows);
					Vector<double> c(A.rows);
					SymbolicCholeskyDecomposition(A, L, Lt, -1);
					FillCholeskyDecomposition(A, L, Lt, solver_threads);
					LowerTriangularSystem(L, c, b);
					UpperTriangularSystem(Lt, x, c);
				}
				break;
			}
			case 3: // Cholesky2_decomposition
			{
				if (solver_auto_reorder)
				{
					CSRMatrix<double> Ar(A.rows, A.columns);
					Vector<double> xr(x.size);
					Vector<double> br(b.size);
					Vector<int> index;
					Vector<int> inverse_index;
					FindingAnOrdering(A, index, inverse_index);
					Reorder(index, inverse_index, A, Ar, solver_threads);
					Reorder(index, b, br, solver_threads);
					CSRMatrix<double> L(Ar.rows, Ar.columns);
					CSRMatrix<double> Lt(Ar.columns, Ar.rows);
					Vector<double> Dr(Ar.rows);
					Vector<double> cr(Ar.rows);
					Vector<double> dr(Ar.rows);
					SymbolicCholeskyDecomposition(Ar, L, Lt, -1);
					FillCholesky2Decomposition(Ar, L, Dr, Lt, solver_threads);
					LowerTriangularSystem(L, cr, br);
					DiagonalSystem(Dr, dr, cr);
					UpperTriangularSystem(Lt, xr, cr);
					Reorder(inverse_index, xr, x, solver_threads);
				}
				else
				{
					CSRMatrix<double> L(A.rows, A.columns);
					CSRMatrix<double> Lt(A.columns, A.rows);
					Vector<double> D(A.rows);
					Vector<double> c(A.rows);
					Vector<double> d(A.rows);
					SymbolicCholeskyDecomposition(A, L, Lt, -1);
					FillCholesky2Decomposition(A, L, D, Lt, solver_threads);
					LowerTriangularSystem(L, c, b);
					DiagonalSystem(D, d, c);
					UpperTriangularSystem(Lt, x, d);
				}
				break;
			}
			case 4: // Biconjugate_gradient
			{
				x.Fill(0);
				switch (preconditioner_type)
				{
					case 0: // None
					{
						CSRMatrix<double> At(A.columns, A.rows);
						Transpose(A, At);
						BiconjugateGradient(A, At, x, b, solver_tolerance, solver_max_steps, solver_threads);
						break;
					}
					case 1: // Jacobi
					{
						CSRMatrix<double> At(A.columns, A.rows);
						Vector<double> M(A.rows);
						Transpose(A, At);
						for (int i = 1; i <= A.rows; ++i)
						{
							M.entry[i] = A(i, i);
						}
						BiconjugateGradientJacobi(A, At, x, b, solver_tolerance, solver_max_steps, M, solver_threads);
						break;
					}
					case 4: // Incomplete_LU
					{
						CSRMatrix<double> At(A.columns, A.rows);
						CSRMatrix<double> L(A.rows, A.columns);
						CSRMatrix<double> U(A.columns, A.rows);
						CSRMatrix<double> Lt(A.columns, A.rows);
						CSRMatrix<double> Ut(A.rows, A.columns);
						SymbolicCholeskyDecomposition(A, L, U, preconditioner_level);
						FillLUDecomposition(A, L, U, solver_threads);
						BiconjugateGradientIncompleteLU(A, At, x, b, solver_tolerance, solver_max_steps, L, U, Lt, Ut, solver_threads);
						break;
					}
					case 5: // Sparse_Approximate_Inverse
					{
						CSCMatrix<double> Ac(A.rows, A.columns);
						CSRMatrix<double> G(A.rows, A.columns);
						CSRMatrix<double> Gt(A.columns, A.rows);
						Transpose(A, Ac);
						SparseApproximateInverseUnsymmetric(A, G, Gt, preconditioner_threshold, preconditioner_level, solver_threads);
						BiconjugateGradientSparseApproximateInverse(A, Ac, x, b, solver_tolerance, solver_max_steps, G, Gt, solver_threads);
						break;
					}
					default:
					{
						Log(0, "Error: Invalid preconditioner %i", preconditioner_type);
						return 1;
					}
				}
				break;
			}
			case 5: // LU_decomposition
			{
				if (solver_auto_reorder)
				{
					CSRMatrix<double> Ar(A.rows, A.columns);
					Vector<double> xr(x.size);
					Vector<double> br(b.size);
					Vector<int> index;
					Vector<int> inverse_index;
					FindingAnOrdering(A, index, inverse_index);
					Reorder(index, inverse_index, A, Ar, solver_threads);
					Reorder(index, b, br, solver_threads);
					CSRMatrix<double> L(Ar.rows, Ar.columns);
					CSRMatrix<double> U(Ar.columns, Ar.rows);
					Vector<double> cr(Ar.rows);
					SymbolicCholeskyDecomposition(Ar, L, U, -1);
					FillLUDecomposition(Ar, L, U, solver_threads);
					LowerTriangularSystem(L, cr, br);
					UpperTriangularSystem(U, xr, cr);
					Reorder(inverse_index, xr, x, solver_threads);
				}
				else
				{
					CSRMatrix<double> L(A.rows, A.columns);
					CSRMatrix<double> U(A.columns, A.rows);
					Vector<double> c(A.rows);
					SymbolicCholeskyDecomposition(A, L, U, -1);
					FillLUDecomposition(A, L, U, solver_threads);
					LowerTriangularSystem(L, c, b);
					UpperTriangularSystem(U, x, c);
				}
				break;
			}
			default:
			{
				Log(0, "Error: Invalid solver");
				return 1;
			}
		}

		// Check for a valid solution and calculate residual
		int invalid = 0;
		double residual = 0;
		#pragma omp parallel for default(shared) reduction(+:invalid,residual)
		for (int i = 1; i <= x.size; ++i)
		{
			if ((x.entry[i] == Float<double>::infinite) || (x.entry[i] == -Float<double>::infinite) || Float<double>::IsNaN(x.entry[i]))
			{
				++invalid;
			}

			double sum = 0;
			int k_max = A.Count(i);
			for (register int k = 1; k <= k_max; ++k)
			{
				sum += A.entry[i][k]*x.entry[A.index[i][k]];
			}
			sum -= b.entry[i];
			residual += sum*sum;
		}
		Log(1, "Solution: %s", (invalid == 0) ? "valid" : "INVALID");
		Log(1, "Residual: %g", sqrt(residual));

		// Save result
		MatFile output_mat_file;
		output_mat_file.Create(argv[2]);
		output_mat_file.Store("x", x);
		output_mat_file.Close();

		if (Memory::memory_usage)
		{
			Log(1, "Peak allocated memory: %lu bytes", (unsigned long)Memory::peak_usage);
		}
	}
	catch (Exception&)
	{
		DebugPosition("Catch fatal exception");
	}

	if (Memory::current_usage != 0)
	{
		fprintf(stderr, "[Error] Memory leak: %lu bytes\n", (unsigned long)Memory::current_usage);
	}

	return 0;
}
