// libFEMSolver.cpp
// Copyright (C) 2014 Miguel Vargas (miguel.vargas@gmail.com)
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
#include <FiniteElement/Assembler.h>
#include <FiniteElement/Mesh.h>
#include <Solver/SolverParameters.h>

#include <stdlib.h>


#define NUMBER_OF_SOLVERS 10

#if defined(_WIN32)
	#define FEMSolverInit FEMSOLVERINIT
	#define FEMSolverEnd FEMSOLVEREND
	#define FEMSolverSetConnectivity FEMSOLVERSETCONNECTIVITY
	#define FEMSolverFillA FEMSOLVERFILLA
	#define FEMSolverSetAe FEMSOLVERSETAE
	#define FEMSolverSetFixed FEMSOLVERSETFIXED
	#define FEMSolverSetB FEMSOLVERSETB
	#define FEMSolverCompensateFixed FEMSOLVERCOMPENSATEFIXED
	#define FEMSolverRun FEMSOLVERRUN
#else
	#define FEMSolverInit femsolverinit_
	#define FEMSolverEnd femsolverend_
	#define FEMSolverSetConnectivity femsolversetconnectivity_
	#define FEMSolverFillA femsolverfilla_
	#define FEMSolverSetAe femsolversetae_
	#define FEMSolverSetFixed femsolversetfixed_
	#define FEMSolverSetB femsolversetb_
	#define FEMSolverCompensateFixed femsolvercompensatefixed_
	#define FEMSolverRun femsolverrun_
#endif


extern "C"
{
	struct Data
	{
		// Solver data
		SolverParameters<double>* solver_parameters;
		Solver<double>* solver;

		// Connectivity data
		int N; // Number of nodes
		int E; // Number of elements
		int V; // Nodes per element
		int D; // Degrees of freedom per node
		int M; // Total number of degrees of freedom
		int Ne; // Degrees of freedon of the elemental matrix
		Mesh mesh;
		Vector<int> element_index;
		Vector<int> node_index;
		Assembler* assembler;

		// System of equations
		CSRMatrix<double> A;
		Vector<double> x;
		Vector<double> b;
		Vector<bool> fixed;
		Matrix<double> Ae;
	};


	Data all_data[NUMBER_OF_SOLVERS];


	// solver_type
	//     1=Conjugate_gradient, 2=Cholesky_decomposition, 3=Cholesky2_decomposition, 4=Biconjugate_gradient, 5=LU_decomposition
	// preconditioner_type
	//     0=None, 1=Jacobi, 2=Incomplete_Cholesky, 3=Incomplete_Cholesky2, 4=Incomplete_LU, 5=Sparse_Approximate_Inverse
	void FEMSolverInit(int* id, int* solver_type, int* solver_threads, double* solver_tolerance, int* solver_max_steps, int* preconditioner_type, int* preconditioner_level, double* preconditioner_threshold, int* message_level) throw()
	{
		try
		{
			Assert((*id >= 1) && (*id <= NUMBER_OF_SOLVERS));
			Data& data = all_data[*id - 1];

			data.solver_parameters = new SolverParameters<double>((SolverType)*solver_type, *solver_threads, *solver_tolerance, *solver_max_steps, (PreconditionerType)*preconditioner_type, *preconditioner_level, *preconditioner_threshold);
			data.solver = (Solver<double>*)0;
			data.assembler = (Assembler*)0;
			log_level = *message_level;
		}
		catch (Memory::Exception&)
		{
			Log(0, "Error allocating memory (FEMSolverInit)");
			exit(1);
		}
		catch (...)
		{
			Log(0, "Unexpected exception (FEMSolverInit)");
			exit(1);
		}
	}


	void FEMSolverEnd(int* id) throw()
	{
		Assert((*id >= 1) && (*id <= NUMBER_OF_SOLVERS));
		Data& data = all_data[*id - 1];

		delete data.solver;
		data.solver = (Solver<double>*)0;

		delete data.assembler;
		data.assembler = (Assembler*)0;

		delete data.solver_parameters;
		data.solver_parameters = (SolverParameters<double>*)0;
	}


	// element_type
	//   2=Triangle, 3=Quadrilateral, 4=Tetrahedra, 5=Hexahedra
	void FEMSolverSetConnectivity(int* id, int* number_of_nodes, int* number_of_elements, int* element_type, int* nodes_per_element, int* degrees_of_freedom, int* connectivity) throw()
	{
		Assert((*id >= 1) && (*id <= NUMBER_OF_SOLVERS));
		try
		{
			Data& data = all_data[*id - 1];

			data.N = *number_of_nodes;
			data.E = *number_of_elements;
			data.V = *nodes_per_element;
			data.D = *degrees_of_freedom;

			data.mesh.Resize((ShapeType)*element_type, data.V, data.E, data.N);
			for (int v = 1; v <= data.V; ++v)
			{
				for (int e = 1; e <= data.E; ++e, ++connectivity)
				{
					data.mesh.connectivity.entry[e][v] = *connectivity;
				}
			}

			data.element_index.Resize(data.E);
			data.element_index.FillSeries(1, 1); // We will use all elements (no partitioning)

			data.node_index.Resize(data.N);
			data.mesh.GetElementsNodes(data.element_index, data.solver_parameters->reorder, data.node_index); // Get nodes for elements (reorder if necesary)

			data.assembler = new Assembler(data.mesh, data.element_index, data.node_index, data.D); // Assembler class
			if (!data.assembler)
			{
				Throw(Memory::exception);
			}

			data.M = data.N*data.D;
			data.A.Resize(data.M, data.M);
			data.assembler->AllocateMatrix(data.A); // Allocate memory for stiffness matrix (reordered if necesary)
			data.x.Resize(data.M);
			data.b.Resize(data.M);
			data.fixed.Resize(data.M);
			data.fixed.Fill(false);

			data.Ne = data.V*data.D;
			data.Ae.Resize(data.Ne, data.Ne);

			data.mesh.PrintInfo();
		}
		catch (Memory::Exception&)
		{
			Log(0, "Error allocating memory (FEMSolverSetConnectivity)");
			exit(1);
		}
		catch (...)
		{
			Log(0, "Unexpected exception (FEMSolverSetConnectivity)");
			exit(1);
		}
	}

	
	void FEMSolverFillA(int* id, double* value) throw()
	{
		Assert((*id >= 1) && (*id <= NUMBER_OF_SOLVERS));
		Data& data = all_data[*id - 1];

		data.A.Fill(*value);
	}


	void FEMSolverSetAe(int* id, int* element_id, double* Ke) throw()
	{
		Assert((*id >= 1) && (*id <= NUMBER_OF_SOLVERS));
		Data& data = all_data[*id - 1];

		Assert(data.assembler);
		Assert((*element_id >= 1) && (*element_id <= data.mesh.elements_count));

		for (int j = 1; j <= data.Ae.columns; ++j)
		{
			for (int i = 1; i <= data.Ae.rows; ++i, ++Ke)
			{
				data.Ae.entry[i][j] = *Ke;
			}
		}
		data.assembler->AssembleAe(*element_id, data.Ae, data.A);
	}


	void FEMSolverSetFixed(int* id, bool* global_fixed) throw()
	{
		Assert((*id >= 1) && (*id <= NUMBER_OF_SOLVERS));

		try
		{
			Data& data = all_data[*id - 1];

			Assert(data.assembler);

			Vector<bool> g_fixed(global_fixed, data.fixed.size);
			data.assembler->AssembleV(g_fixed, data.fixed);

			data.solver = data.solver_parameters->Instanciate(data.A, data.x, data.b, data.fixed, true);
			data.solver_parameters->PrintInfo();
		}
		catch (Memory::Exception&)
		{
			Log(0, "Error allocating memory (FEMSolverSetFixed)");
			exit(1);
		}
		catch (...)
		{
			Log(0, "Unexpected exception (FEMSolverSetFixed)");
			exit(1);
		}
	}


	void FEMSolverSetB(int* id, double* global_b) throw()
	{
		Assert((*id >= 1) && (*id <= NUMBER_OF_SOLVERS));

		Data& data = all_data[*id - 1];

		Assert(data.assembler);
		Assert(data.solver);

		Vector<double> g_b(global_b, data.b.size);
		data.assembler->AssembleV(g_b, data.b);
	}


	void FEMSolverCompensateFixed(int* id, double* global_fixed_values, double* global_b) throw()
	{
		Assert((*id >= 1) && (*id <= NUMBER_OF_SOLVERS));

		Data& data = all_data[*id - 1];

		Assert(data.assembler);
		Assert(data.solver);

		Vector<double> g_x(global_fixed_values, data.x.size);
		data.assembler->AssembleV(g_x, data.x);

		Vector<double> g_b(global_b, data.b.size);
		data.assembler->AssembleV(g_b, data.b);

		data.solver->CompensateFixed();
	}


	void FEMSolverRun(int* id, double* global_x, bool* valid_solution) throw()
	{
		Assert((*id >= 1) && (*id <= NUMBER_OF_SOLVERS));

		try
		{
			Data& data = all_data[*id - 1];

			*valid_solution = data.solver->Calculate();
			Log(1, "Solution: %s", valid_solution ? "valid" : "INVALID");

			Vector<double> g_x(global_x, data.x.size);

			// Store solution
			for (int n = 1; n <= data.N; ++n)
			{
				int i = (n - 1)*data.D;
				int gi = (data.node_index.entry[n] - 1)*data.D;
				for (int d = 1; d <= data.D; ++d)
				{
					g_x.entry[gi + d] = data.x.entry[i + d]; // Fill solution (inverse order if necesary)
				}
			}
		}
		catch (Memory::Exception&)
		{
			Log(0, "Error allocating memory (FEMSolverRun)");
			exit(1);
		}
		catch (...)
		{
			Log(0, "Unexpected exception (FEMSolverRun)");
			exit(1);
		}
	}
}
