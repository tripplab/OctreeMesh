// SchurNS.h
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

#ifndef _SchurNS_h_
#define _SchurNS_h_

#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Communication/MPI.h>
#include <Container/CSRMatrix.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>
#include <FiniteElement/Partition.h>
#include <Solver/LU.h>
#include <Solver/SolverParameters.h>
#include <Solver/TriangularSystem.h>


enum
{
	tag_schur_Mi,
	tag_schur_Bi,
	tag_schur_K_count,
	tag_schur_K_index,
	tag_schur_K_entry,
	tag_schur_u,
	tag_schur_f,
	tag_schur_fixed,
	tag_schur_fBi_adjust,
	tag_schur_barfBi,
	tag_schur_MBi,
	tag_schur_vBi,
	tag_schur_wBi,
	tag_schur_uBi,
	tag_schur_uIi,
	tag_schur_end,
	tag_schur_memory_peak_usage = 17
};


template <typename T>
class SchurNSMaster
{
	public:

		const MPI& mpi;
		const int partitions_count;

		const Vector<Substructure>& substructures;
		const Boundary& boundary;
		const int degrees_of_freedom;

		const int B;
		Vector<T> fB_adjust;
		Vector<T> MB;


		SchurNSMaster(const MPI& mpi, const Vector<Substructure>& substructures, const Boundary& boundary, const int degrees_of_freedom) throw(Memory::Exception, MPI::Exception)
		:	mpi(mpi),
			partitions_count(mpi.size - 1),
			substructures(substructures),
			boundary(boundary),
			degrees_of_freedom(degrees_of_freedom),
			B(boundary.node_index.size*degrees_of_freedom),
			fB_adjust(B),
			MB(B)
		{
		}


		void SendSystemOfEquations(const int partition_id, const CSRMatrix<T>& K, const Vector<T>& u, const Vector<T>& f, const Vector<bool>& fixed) const throw(MPI::Exception)
		{
			Assert(partition_id >= 1);
			Assert(partition_id <= partitions_count);
			Assert(u.size == K.rows);
			Assert(f.size == K.rows);
			Assert(fixed.size == K.rows);

			try
			{
				int Mi = K.rows;
				int Bi = substructures.entry[partition_id].boundary_count*degrees_of_freedom;

				// Send system size
				#pragma omp critical
				{
					mpi.Send(partition_id, tag_schur_Mi, Mi);
					mpi.Send(partition_id, tag_schur_Bi, Bi);
				}

				// Send vectors
				#pragma omp critical
				{
					mpi.Send(partition_id, tag_schur_u, u);
				}
				#pragma omp critical
				{
					mpi.Send(partition_id, tag_schur_f, f);
				}
				#pragma omp critical
				{
					mpi.Send(partition_id, tag_schur_fixed, fixed);
				}

				// Send matrix
				for (int i = 1; i <= Mi; ++i)
				{
					#pragma omp critical
					{
						mpi.Send(partition_id, tag_schur_K_count, K.Count(i));
						mpi.Send(partition_id, tag_schur_K_index, K.Count(i), &K.index[i][1]);
						mpi.Send(partition_id, tag_schur_K_entry, K.Count(i), &K.entry[i][1]);
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void InitializeSolver(CSRMatrix<T>& KBB, Vector<T>& uB, Vector<T>& fB, Vector<bool>& fixedB) throw(Memory::Exception, MPI::Exception)
		{
			Assert(KBB.rows == B);
			Assert(KBB.rows == KBB.columns);

			try
			{
				// Apply fixed conditions and compensate in fB_adjust
				Vector<int> fixed(B);
				fB_adjust.Resize(B);
				fB_adjust.Fill(0);
				for (int i = 1; i <= B; ++i)
				{
					if (fixedB.entry[i])
					{
						KBB.AllocateRow(i, 1);
						KBB.index[i][1] = i;
						KBB.entry[i][1] = 1.0;
					}
					else
					{
						int fixed_count = 0;
						int k_max = KBB.Count(i);
						for (register int k = 1; k <= k_max; ++k)
						{
							register int j = KBB.index[i][k];
							if (fixedB.entry[j])
							{
								fB_adjust.entry[i] += KBB.entry[i][k]*uB.entry[j];
								fixed.entry[++fixed_count] = j;
							}
						}
						for (register int f = 1; f <= fixed_count; ++f)
						{
							KBB.RemoveEntry(i, fixed.entry[f]);
						}
					}
				}

				// Assemble preconditioner for KBB
				for (int i = 1; i <= B; ++i)
				{
					MB.entry[i] = KBB(i, i);
				}

				// Complete fBi_adjust, barfBi and MBi
				for (int p = 1; p <= partitions_count; ++p)
				{
					int node_count = substructures.entry[p].node_index.size;
					int boundary_count = substructures.entry[p].boundary_count;
					int internal_count = node_count - boundary_count;

					int Bi = boundary_count*degrees_of_freedom;

					// Receive fBi_adjust, barfBi and MBi
					Vector<T> fBi_adjust(Bi);
					Vector<T> barfBi(Bi);
					Vector<T> MBi(Bi);
					mpi.Receive(p, tag_schur_fBi_adjust, fBi_adjust);
					mpi.Receive(p, tag_schur_barfBi, barfBi);
					mpi.Receive(p, tag_schur_MBi, MBi);
					for (int i = 1; i <= boundary_count; ++i)
					{
						int n = substructures.entry[p].node_index.entry[internal_count + i];
						int bB = (boundary.inverse_node_index.entry[n] - 1)*degrees_of_freedom;
						int bBi = (i - 1)*degrees_of_freedom;

						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							fB_adjust.entry[bB + d] += fBi_adjust.entry[bBi + d];
							fB.entry[bB + d] -= barfBi.entry[bBi + d];
							MB.entry[bB + d] -= MBi.entry[bBi + d];
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void ComplesateFixed(const Vector<T>& uB, Vector<T>& fB, const Vector<bool>& fixedB, int threads) throw()
		{
			#pragma omp parallel for default(shared) num_threads(threads)
			for (int i = 1; i <= fB.size; ++i)
			{
				if (!fixedB.entry[i])
				{
					fB.entry[i] -= fB_adjust.entry[i];
				}
				else
				{
					fB.entry[i] = uB.entry[i];
				}
			}
		}


		void MultiplyKBB(const CSRMatrix<T>& KBB, const Vector<T>& vB, Vector<T>& wB, int threads) throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				// Contribution of KBB
				#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
				for (int i = 1; i <= B; ++i)
				{
					int* __restrict KBB_index_i = KBB.index[i];
					T* __restrict KBB_entry_i = KBB.entry[i];

					register T sum = 0.0;
					const int k_max = KBB.Count(i);
					for (register int k = 1; k <= k_max; ++k)
					{
						sum += KBB_entry_i[k]*vB.entry[KBB_index_i[k]];
					}
					wB.entry[i] = sum;
				}

				// Contribution of barKBBi
				Vector<MPI::DataRequest> request_vBi(partitions_count);
				Vector<Vector<T> > partition_vBi(partitions_count);
				bool end = false;
				for (int p = 1; p <= partitions_count; ++p)
				{
					int boundary_count = substructures.entry[p].boundary_count;
					int internal_count = substructures.entry[p].node_index.size - boundary_count;

					int Bi = boundary_count*degrees_of_freedom;

					// Fill fixedBi
					Vector<T>& vBi = partition_vBi.entry[p];
					vBi.Resize(Bi);
					for (int b = 1; b <= boundary_count; ++b)
					{
						int n = substructures.entry[p].node_index.entry[internal_count + b];
						int bB = (boundary.inverse_node_index.entry[n] - 1)*degrees_of_freedom;
						int bBi = (b - 1)*degrees_of_freedom;

						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							vBi.entry[bBi + d] = vB.entry[bB + d];
						}
					}

					// Send non stop condition
					mpi.Send(p, tag_schur_end, end);

					// Send vBi
					mpi.Send(p, tag_schur_vBi, vBi);

					// Request to receive vBi with result
					mpi.RequestReceive(p, tag_schur_wBi, vBi, request_vBi.entry[p]);
				}
				for (int q = 1; q <= partitions_count; ++q)
				{
					int p = mpi.WaitAnyRequest(request_vBi);

					int boundary_count = substructures.entry[p].boundary_count;
					int internal_count = substructures.entry[p].node_index.size - boundary_count;

					// Assemble preconditioner part from barKBBi
					Vector<T>& vBi = partition_vBi.entry[p];
					for (int b = 1; b <= boundary_count; ++b)
					{
						int n = substructures.entry[p].node_index.entry[internal_count + b];
						int bB = (boundary.inverse_node_index.entry[n] - 1)*degrees_of_freedom;
						int bBi = (b - 1)*degrees_of_freedom;

						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							wB.entry[bB + d] -= vBi.entry[bBi + d];
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void MultiplyKBBTanspose(const CSRMatrix<T>& KBBt, const Vector<T>& vB, Vector<T>& wB, int threads) throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				// Contribution of KBB
				#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
				for (int i = 1; i <= B; ++i)
				{
					int* __restrict KBB_index_i = KBBt.index[i];
					T* __restrict KBBt_entry_i = KBBt.entry[i];

					register T sum = 0.0;
					const int k_max = KBBt.Count(i);
					for (register int k = 1; k <= k_max; ++k)
					{
						sum += KBBt_entry_i[k]*vB.entry[KBBt_index_i[k]];
					}
					wB.entry[i] = sum;
				}

				// Contribution of barKBBti
				Vector<MPI::DataRequest> request_vBi(partitions_count);
				Vector<Vector<T> > partition_vBi(partitions_count);
				bool end = false;
				for (int p = 1; p <= partitions_count; ++p)
				{
					int boundary_count = substructures.entry[p].boundary_count;
					int internal_count = substructures.entry[p].node_index.size - boundary_count;

					int Bi = boundary_count*degrees_of_freedom;

					// Fill fixedBi
					Vector<T>& vBi = partition_vBi.entry[p];
					vBi.Resize(Bi);
					for (int b = 1; b <= boundary_count; ++b)
					{
						int n = substructures.entry[p].node_index.entry[internal_count + b];
						int bB = (boundary.inverse_node_index.entry[n] - 1)*degrees_of_freedom;
						int bBi = (b - 1)*degrees_of_freedom;

						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							vBi.entry[bBi + d] = vB.entry[bB + d];
						}
					}

					// Send vBi
					mpi.Send(p, tag_schur_vBi, vBi);

					// Request to receive vBi with result
					mpi.RequestReceive(p, tag_schur_wBi, vBi, request_vBi.entry[p]);
				}
				for (int q = 1; q <= partitions_count; ++q)
				{
					int p = mpi.WaitAnyRequest(request_vBi);

					int boundary_count = substructures.entry[p].boundary_count;
					int internal_count = substructures.entry[p].node_index.size - boundary_count;

					// Assemble preconditioner part from barKBBti
					Vector<T>& vBi = partition_vBi.entry[p];
					for (int b = 1; b <= boundary_count; ++b)
					{
						int n = substructures.entry[p].node_index.entry[internal_count + b];
						int bB = (boundary.inverse_node_index.entry[n] - 1)*degrees_of_freedom;
						int bBi = (b - 1)*degrees_of_freedom;

						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							wB.entry[bB + d] -= vBi.entry[bBi + d];
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		// BiConjugate gradient width Jacobi preconditioner
		void Solver(CSRMatrix<T>& A, CSRMatrix<T>& At, Vector<T>& x, Vector<T>& b, T tolerance, int max_steps, int threads) throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				int N = A.rows;

				Vector<T> r1(N); // Residual 1
				Vector<T> r2(N); // Residual 2
				Vector<T> p1(N); // Descent direcction 1
				Vector<T> p2(N); // Descent direcction 2
				Vector<T> z1(N); // Preconditioner system solution 1
				Vector<T> z2(N); // Preconditioner system solution 2
				Vector<T> w1(N); // w1 = A*p
				Vector<T> w2(N); // w2 = A'*p2

				MultiplyKBB(A, x, r1, threads);
				T r1r1 = 0;
				T r2z1 = 0;
				#pragma omp parallel for default(shared) reduction(+:r2z1,r1r1) schedule(guided) num_threads(threads)
				for (int i = 1; i <= N; ++i)
				{
					r1.entry[i] = b.entry[i] - r1.entry[i];  // r1 = b - Ax;
					r2.entry[i] = r1.entry[i];       // r2 = r1
					z1.entry[i] = r1.entry[i]/MB.entry[i];  // Solve for q: MB*z1 = r1
					z2.entry[i] = r2.entry[i]/MB.entry[i];  // Solve for q: MB'*z2 = r2
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
					// Test termination condition
					if (r1r1 <= epsilon) // r1'*r1 <= tolerance
					{
						break;
					}

					MultiplyKBB(A, p1, w1, threads);
					MultiplyKBBTanspose(A, p1, w1, threads);
					T p2w1 = 0.0;
					#pragma omp parallel for default(shared) reduction(+:p2w1) schedule(guided) num_threads(threads)
					for (int i = 1; i <= N; ++i)
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
					for (int i = 1; i <= N; ++i)
					{
						x.entry[i] += alpha*p1.entry[i];   // Xn = x + alpha*p
						r1.entry[i] -= alpha*w1.entry[i];  // r1n = r1 - alpha*w1
						r2.entry[i] -= alpha*w2.entry[i];  // r2n = r2 - alpha*w2
						z1.entry[i] = r1.entry[i]/MB.entry[i];    // Solve for z1: MB*z1 = r1
						z2.entry[i] = r2.entry[i]/MB.entry[i];    // Solve for z2: MB'*z2 = r2
						r2nz1n += r2.entry[i]*z1.entry[i]; // r2nz1n = r2n'*z1n
						r1r1 += r1.entry[i]*r1.entry[i]; // r1r1 = r1n'*r1n
					}

					T beta = r2nz1n/r2z1; // beta = (r2n'*r1n)/(r2'*r1)

					#pragma omp parallel for default(shared) num_threads(threads)
					for (int i = 1; i <= N; ++i)
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
					return;
				}

				return step;

				Vector<T> g(B); // Gradient
				Vector<T> p(B); // Descent direcction
				Vector<T> w(B); // w = A*p
				Vector<T> q(B); // Preconditioner system solution

				MultiplyKBB(A, x, g, threads);
				T gg = 0.0;
				T gq = 0.0;
				#pragma omp parallel for default(shared) reduction(+:gg,gq) num_threads(threads)
				for (int i = 1; i <= B; ++i)
				{
					g.entry[i] -= b.entry[i]; // g = AX - b
					q.entry[i] = g.entry[i]/MB.entry[i]; // Solve for q: MB*q = g
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

					MultiplyKBB(A, p, w, threads);
					T pw = 0.0;
					#pragma omp parallel for default(shared) reduction(+:pw) num_threads(threads)
					for (int i = 1; i <= B; ++i)
					{
						pw += p.entry[i]*w.entry[i]; // pw = p'*w
					}

					T alpha = gq/pw; // alpha = (g'*q)/(p'*w)

					T gngn = 0.0;
					T gnqn = 0.0;
					#pragma omp parallel for default(shared) reduction(+:gngn,gnqn) num_threads(threads)
					for (int i = 1; i <= B; ++i)
					{
						x.entry[i] += alpha*p.entry[i]; // Xn = x + alpha*p
						g.entry[i] += alpha*w.entry[i]; // Gn = g + alpha*w
						q.entry[i] = g.entry[i]/MB.entry[i]; // Solve for q: MQ = g
						gngn += g.entry[i]*g.entry[i]; // gngn = Gn'*Gn
						gnqn += g.entry[i]*q.entry[i]; // gnqn = Gn'*Qn
					}

					T beta = gnqn/gq; // beta = (Gn'*Qn)/(g'*Q)

					#pragma omp parallel for default(shared) num_threads(threads)
					for (int i = 1; i <= B; ++i)
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
					Log(1, "-[Warning] ConjugateGradientJacobi did not converge in %i steps", max_steps);
				}

				// Send non stop condition
				bool end = true;
				for (int p = 1; p <= partitions_count; ++p)
				{
					mpi.Send(p, tag_schur_end, end);
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void CompleteSolution(const Vector<T>& uB, Vector<T>& u) throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				for (int b = 1; b <= boundary.node_index.size; ++b)
				{
					int nd = (boundary.node_index.entry[b] - 1)*degrees_of_freedom;
					int bd = (b - 1)*degrees_of_freedom;
					for (int d = 1; d <= degrees_of_freedom; ++d)
					{
						u.entry[nd + d] = uB.entry[bd + d];
					}
				}

				Vector<MPI::DataRequest> request_uIi(partitions_count);
				Vector<Vector<T> > partition_uIi(partitions_count);
				for (int p = 1; p <= partitions_count; ++p)
				{
					int node_count = substructures.entry[p].node_index.size;
					int boundary_count = substructures.entry[p].boundary_count;
					int internal_count = node_count - boundary_count;

					int Mi = node_count*degrees_of_freedom;
					int Bi = boundary_count*degrees_of_freedom;
					int Ii = Mi - Bi;

					// Fill uBi
					Vector<T> uBi(Bi);
					for (int b = 1; b <= boundary_count; ++b)
					{
						int n = substructures.entry[p].node_index.entry[internal_count + b];
						int bB = (boundary.inverse_node_index.entry[n] - 1)*degrees_of_freedom;
						int bBi = (b - 1)*degrees_of_freedom;

						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							uBi.entry[bBi + d] = uB.entry[bB + d];
						}
					}

					// Send uBi
					mpi.Send(p, tag_schur_uBi, uBi);

					// Request to receive uIi
					Vector<T>& uIi = partition_uIi.entry[p];
					uIi.Resize(Ii);
					mpi.RequestReceive(p, tag_schur_uIi, uIi, request_uIi.entry[p]);
				}
				for (int q = 1; q <= partitions_count; ++q)
				{
					int p = mpi.WaitAnyRequest(request_uIi);

					int node_count = substructures.entry[p].node_index.size;
					int boundary_count = substructures.entry[p].boundary_count;
					int internal_count = node_count - boundary_count;

					Vector<T>& uIi = partition_uIi.entry[p];

					for (int i = 1; i <= internal_count; ++i)
					{
						int iI = (substructures.entry[p].node_index.entry[i] - 1)*degrees_of_freedom;
						int iIi = (i - 1)*degrees_of_freedom;

						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							u.entry[iI + d] = uIi.entry[iIi + d];
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void GetMemoryPeakUsage(Vector<size_t>& partition_memory_peak_usage) throw(MPI::Exception)
		{
			try
			{
				for (int p = 1; p <= partitions_count; ++p)
				{
					mpi.Receive(p, tag_schur_memory_peak_usage, partition_memory_peak_usage.entry[p]);
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


	private:

		SchurNSMaster& operator = (SchurNSMaster&)
		{
			return *this;
		}
};


template <typename T>
class SchurNSSlave
{
	public:

		const MPI& mpi;

		int Mi;

		int Ii;
		CSRMatrix<T> KIBi;
		Vector<T> uIi;
		Vector<T> fIi;
		Vector<bool> fixedIi;
		Vector<T> fIi_adjust;

		int Bi;
		CSRMatrix<T> KBIi;
		Vector<T> uBi;

		CSRMatrix<T> L;
		CSRMatrix<T> U;
		int threads;

		Matrix<T> barKBBi;
		Matrix<T> barKBBit;


		SchurNSSlave(const MPI& mpi, const int threads) throw(MPI::Exception)
		:	mpi(mpi),
			Mi(),
			Ii(),
			KIBi(),
			uIi(),
			fIi(),
			fixedIi(),
			fIi_adjust(),
			Bi(),
			KBIi(),
			uBi(),
			L(),
			U(),
			threads(threads),
			barKBBi(),
			barKBBit()
		{
			try
			{
				// Receive system size
				mpi.Receive(0, tag_schur_Mi, Mi);
				mpi.Receive(0, tag_schur_Bi, Bi);
				Ii = Mi - Bi;
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void ReceiveMatrix(const int count_tag, const int index_tag, const int entry_tag, CSRMatrix<T>& AIIi, CSRMatrix<T>& AIBi, CSRMatrix<T>& ABIi) throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				// Allocate AIIi and AIBi
				AIIi.Resize(Ii, Ii);
				AIBi.Resize(Ii, Bi);

				// Receive Ii rows and assemble AIIi and AIBi
				Vector<int> index(Mi);
				Vector<double> values(Mi);
				for (int i = 1; i <= Ii; ++i)
				{
					int k_max;
					mpi.Receive(0, count_tag, k_max);
					mpi.Receive(0, index_tag, k_max, index.data);
					mpi.Receive(0, entry_tag, k_max, values.data);

					int size_i = 0;
					int size_b = 0;
					for (int k = 1; k <= k_max; ++k)
					{
						int j = index.entry[k];
						if (j <= Ii)
						{
							++size_i;
						}
						else
						{
							++size_b;
						}
					}

					AIIi.AllocateRow(i, size_i);
					AIBi.AllocateRow(i, size_b);
					for (int k = 1, l = 0, m = 0; k <= k_max; ++k)
					{
						int j = index.entry[k];
						if (j <= Ii)
						{
							++l;
							AIIi.index[i][l] = j;
							AIIi.entry[i][l] = values.entry[k];
						}
						else
						{
							++m;
							AIBi.index[i][m] = j - Ii;
							AIBi.entry[i][m] = values.entry[k];
						}
					}
				}

				// Allocate ABIi
				ABIi.Resize(Bi, Ii);

				// Receive Bi rows and assemble ABIi
				for (int b = 1; b <= Bi; ++b)
				{
					int k_max;
					mpi.Receive(0, count_tag, k_max);
					mpi.Receive(0, index_tag, k_max, index.data);
					mpi.Receive(0, entry_tag, k_max, values.data);

					int size_i = 0;
					for (int k = 1; k <= k_max; ++k)
					{
						int j = index.entry[k];
						if (j <= Ii)
						{
							++size_i;
						}
					}

					ABIi.AllocateRow(b, size_i);
					for (int k = 1, l = 0; k <= k_max; ++k)
					{
						int j = index.entry[k];
						if (j <= Ii)
						{
							++l;
							ABIi.index[b][l] = j;
							ABIi.entry[b][l] = values.entry[k];
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		template <typename U>
		void ReceiveVector(const int tag, Vector<U>& vIi, Vector<U>& vBi) throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				Vector<U> a(Mi);
				mpi.Receive(0, tag, a);

				vIi.Resize(Ii);
				for (int i = 1; i <= Ii; ++i)
				{
					vIi.entry[i] = a.entry[i];
				}

				vBi.Resize(Bi);
				for (int b = 1; b <= Bi; ++b)
				{
					vBi.entry[b] = a.entry[Ii + b];
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		template <typename U>
		void SendVector(const int tag, Vector<U>& vIi, Vector<U>& vBi) throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				Vector<U> a(Mi);
				for (int i = 1; i <= Ii; ++i)
				{
					a.entry[i] = vIi.entry[i];
				}
				for (int b = 1; b <= Bi; ++b)
				{
					a.entry[Ii + b] = vBi.entry[b];
				}
				mpi.Send(0, tag, a);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		bool CheckSolution(const Vector<T>& x) const throw()
		{
			// Check for a valid solution
			int invalid = 0;
			#pragma omp parallel for default(shared) reduction(+:invalid) num_threads(threads)
			for (int i = 1; i <= x.size; ++i)
			{
				if ((x.entry[i] == Float<T>::infinite) || (x.entry[i] == -Float<T>::infinite) || Float<T>::IsNaN(x.entry[i]))
				{
					++invalid;
				}
			}
			return (invalid > 0) ? false : true;
		}


		void ReceiveSystemOfEquations() throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				// Receive and factorize system of equations
				{
					Vector<T> fBi(Bi);
					Vector<bool> fixedBi(Bi);
					Vector<T> fBi_adjust(Bi);

					// Receive vectors
					ReceiveVector(tag_schur_u, uIi, uBi);
					ReceiveVector(tag_schur_f, fIi, fBi);
					ReceiveVector(tag_schur_fixed, fixedIi, fixedBi);

					// Receive KIIi, KIBi and KBIi
					CSRMatrix<T> KIIi;
					ReceiveMatrix(tag_schur_K_count, tag_schur_K_index, tag_schur_K_entry, KIIi, KIBi, KBIi);

					// For removing fixed columns
					Vector<int> fixed(Ii > Bi ? Ii : Bi);

					// Apply fixed conditions compensate in fIi_adjust
					fIi_adjust.Resize(Ii);
					fIi_adjust.Fill(0);
					for (int i = 1; i <= Ii; ++i)
					{
						if (fixedIi.entry[i])
						{
							KIIi.AllocateRow(i, 1);
							KIIi.index[i][1] = i;
							KIIi.entry[i][1] = 1.0;
							KIBi.AllocateRow(i, 0);
						}
						else
						{
							int fixed_count = 0;
							int k_max = KIIi.Count(i);
							for (register int k = 1; k <= k_max; ++k)
							{
								register int j = KIIi.index[i][k];
								if (fixedIi.entry[j])
								{
									fIi_adjust.entry[i] += KIIi.entry[i][k]*uIi.entry[j];
									fixed.entry[++fixed_count] = j;
								}
							}
							for (register int f = 1; f <= fixed_count; ++f)
							{
								KIIi.RemoveEntry(i, fixed.entry[f]);
							}
							fixed_count = 0;
							k_max = KIBi.Count(i);
							for (register int k = 1; k <= k_max; ++k)
							{
								register int j = KIBi.index[i][k];
								if (fixedBi.entry[j])
								{
									fIi_adjust.entry[i] += KIBi.entry[i][k]*uBi.entry[j];
									fixed.entry[++fixed_count] = j;
								}
							}
							for (register int f = 1; f <= fixed_count; ++f)
							{
								KIBi.RemoveEntry(i, fixed.entry[f]);
							}
						}
					}

					// Apply fixed conditions compensate in fBi_adjust
					fBi_adjust.Resize(Bi);
					fBi_adjust.Fill(0);
					for (int i = 1; i <= Bi; ++i)
					{
						if (fixedBi.entry[i])
						{
							KBIi.AllocateRow(i, 0);
						}
						else
						{
							int fixed_count = 0;
							int k_max = KBIi.Count(i);
							for (register int k = 1; k <= k_max; ++k)
							{
								register int j = KBIi.index[i][k];
								if (fixedIi.entry[j])
								{
									fBi_adjust.entry[i] += KBIi.entry[i][k]*uIi.entry[j];
									fixed.entry[++fixed_count] = j;
								}
							}
							for (register int f = 1; f <= fixed_count; ++f)
							{
								KBIi.RemoveEntry(i, fixed.entry[f]);
							}
						}
					}

					// Send contribution for fB_adjust
					mpi.Send(0, tag_schur_fBi_adjust, fBi_adjust); 

					// Factorize KII
					L.Resize(Ii, Ii);
					U.Resize(Ii, Ii);
					SymbolicCholeskyDecomposition(KIIi, L, U);
					FillLUDecomposition(KIIi, L, U, threads);
				}

				// Calculate barKBBi
				barKBBi.Resize(Bi, Bi);
				Vector<T> t(Ii);
				Vector<T> c(Ii);
				Vector<T> d(Ii);

				CSCMatrix<T> KIBi_csc;
				Convert(KIBi, KIBi_csc);
				c.Fill(0);
				for (int i = 1; i <= Bi; ++i)
				{
					int k_max = KIBi_csc.Count(i);
					for (register int k = 1; k <= k_max; ++k)
					{
						int j = KIBi_csc.index[i][k];
						c.entry[j] = KIBi_csc.entry[i][k];
					}
					LowerTriangularSystem(L, d, c);
					UpperTriangularSystem(U, t, d);
					if (!CheckSolution(t))
					{
						Log(0, "[Error] Inverting KII (partition %i).", mpi.rank);
					}
					for (register int k = 1; k <= k_max; ++k)
					{
						int j = KIBi_csc.index[i][k];
						c.entry[j] = 0;
					}

					for (int b = 1; b <= Bi; ++b)
					{
						T sum = 0;
						int k_max = KBIi.count[b];
						for (register int k = 1; k <= k_max; ++k)
						{
							int j = KBIi.index[b][k];
							sum += KBIi.entry[b][k]*t.entry[j];
						}
						barKBBi.entry[b][i] = sum;
					}
				}
				
				// Send barfBi
				for (int i = 1; i <= Ii; ++i)
				{
					c.entry[i] = fIi.entry[i] - fIi_adjust.entry[i];
				}
				Vector<T> barfBi(Bi);
				LowerTriangularSystem(L, d, c);
				UpperTriangularSystem(U, t, d);
				if (!CheckSolution(t))
				{
					Log(0, "[Error] Calculating bar_fBi (partition %i).", mpi.rank);
				}
				for (int b = 1; b <= Bi; ++b)
				{
					T sum = 0;
					int k_max = KBIi.count[b];
					for (register int k = 1; k <= k_max; ++k)
					{
						int j = KBIi.index[b][k];
						sum += KBIi.entry[b][k]*t.entry[j];
					}
					barfBi.entry[b] = sum;
				}
				mpi.Send(0, tag_schur_barfBi, barfBi);

				// Send preconditioner
				Vector<T> MBi(Bi);
				for (int i = 1; i <= Bi; ++i)
				{
					MBi.entry[i] = barKBBi.entry[i][i];
				}
				mpi.Send(0, tag_schur_MBi, MBi);

				// Store barKBBit
				barKBBit.Resize(Bi, Bi);
				Transpose(barKBBi, barKBBit);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void MultiplybarKBB(Vector<T>& barvBi, Vector<T>& barwBi) throw(MPI::Exception)
		{
			try
			{
				// Receive barvBi
				mpi.Receive(0, tag_schur_vBi, barvBi);

				//barKBBi
				#pragma omp parallel for default(shared) num_threads(threads)
				for (int i = 1; i <= Bi; ++i)
				{
					const T* __restrict barKBBi_entry_i = barKBBi.entry[i];
					register T sum = 0;
					for (register int j = 1; j <= Bi; ++j)
					{
						sum += barKBBi_entry_i[j]*barvBi.entry[j];
					}
					barwBi.entry[i] = sum;
				}

				// Send barvBi
				mpi.Send(0, tag_schur_wBi, barwBi);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void MultiplybarKBBTranspose(Vector<T>& barvBi, Vector<T>& barwBi) throw(MPI::Exception)
		{
			try
			{
				// Receive barvBi
				mpi.Receive(0, tag_schur_vBi, barvBi);

				//barKBBit
				#pragma omp parallel for default(shared) num_threads(threads)
				for (int i = 1; i <= Bi; ++i)
				{
					const T* __restrict barKBBit_entry_i = barKBBit.entry[i];
					register T sum = 0;
					for (register int j = 1; j <= Bi; ++j)
					{
						sum += barKBBit_entry_i[j]*barvBi.entry[j];
					}
					barwBi.entry[i] = sum;
				}

				// Send barvBi
				mpi.Send(0, tag_schur_wBi, barwBi);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void Solver() throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				// Allocate space for working vectors
				Vector<T> barvBi(Bi);
				Vector<T> barwBi(Bi);

				// Multiplications
				for ( ; ; )
				{
					bool end;
					mpi.Receive(0, tag_schur_end, end);
					if (end)
					{
						break;
					}
					MultiplybarKBB(barvBi, barwBi);
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void CompleteSolution() throw(MPI::Exception)
		{
			try
			{
				// Receive uBi
				mpi.Receive(0, tag_schur_uBi, uBi);

				// Calculate uIi
				Vector<T> c(fIi);
				#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
				for (int i = 1; i <= Ii; ++i)
				{
					T sum = 0;
					int k_max = KIBi.Count(i);
					for (register int k = 1; k <= k_max; ++k)
					{
						register int j = KIBi.index[i][k];
						sum += KIBi.entry[i][k]*uBi.entry[j];
					}
					if (!fixedIi.entry[i])
					{
						c.entry[i] -= sum + fIi_adjust.entry[i];
					}
					else
					{
						c.entry[i] = uIi.entry[i];
					}
				}

				// Solve the system
				Vector<T> d(Ii);
				LowerTriangularSystem(L, d, c);
				UpperTriangularSystem(U, uIi, d);
				if (!CheckSolution(uIi))
				{
					Log(0, "[Error] Calculating uI (partition %i).", mpi.rank);
				}

				// Send uIi
				mpi.Send(0, tag_schur_uIi, uIi);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void ReportMemoryPeakUsage() throw(MPI::Exception)
		{
			try
			{
				// Send peak memory usage
				mpi.Send(0, tag_schur_memory_peak_usage, Memory::peak_usage);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


	private:

		SchurNSSlave& operator = (SchurNSSlave&) throw()
		{
			return *this;
		}
};

#endif
