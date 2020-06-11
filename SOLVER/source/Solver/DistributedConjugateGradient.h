// DistributedConjugateGradient.h
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

#ifndef _DistributedConjugateGradient_h_
#define _DistributedConjugateGradient_h_

#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Communication/MPI.h>
#include <Container/CSRMatrix.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>
#include <FiniteElement/Partition.h>


enum
{
	tag_dcg_m,
	tag_dcg_Bi,
	tag_dcg_K_count,
	tag_dcg_K_index,
	tag_dcg_K_entry,
	tag_dcg_u,
	tag_dcg_f,
	tag_dcg_fixed,
	tag_dcg_f_adjust,
	tag_dcg_M,
	tag_dcg_end,
	tag_dcg_vi,
	tag_dcg_wi,
	tag_dcg_memory_peak_usage = 17
};


template <typename T>
class DistributedConjugateGradientMaster
{
	public:

		const MPI& mpi;
		const int partitions_count;

		const Vector<Partition>& partitions;
		const int degrees_of_freedom;

		const int m;
		Vector<T> f_adjust;
		Vector<T> M;
		Vector<bool> fixed_on_boundary; /**/

		Vector<Vector<T> > partition_vi;
		Vector<Vector<T> > partition_wi;
		Vector<MPI::DataRequest> request_wi;


		DistributedConjugateGradientMaster(const MPI& mpi, const Vector<Partition>& partitions, const int degrees_of_freedom, const int nodes_count) throw(Memory::Exception, MPI::Exception)
		:	mpi(mpi),
			partitions_count(mpi.size - 1),
			partitions(partitions),
			degrees_of_freedom(degrees_of_freedom),
			m(nodes_count*degrees_of_freedom),
			f_adjust(m),
			M(m),
			partition_vi(partitions_count),
			partition_wi(partitions_count),
			request_wi(partitions_count)
		{
			for (int p = 1; p <= partitions_count; ++p)
			{
				partition_vi.entry[p].Resize(partitions.entry[p].node_index.size*degrees_of_freedom);
				partition_wi.entry[p].Resize(partitions.entry[p].node_index.size*degrees_of_freedom);
			}
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
				int mi = K.rows;

				#pragma omp critical
				{
					mpi.Send(partition_id, tag_dcg_m, mi);
					mpi.Send(partition_id, tag_dcg_u, u);
					mpi.Send(partition_id, tag_dcg_f, f);
					mpi.Send(partition_id, tag_dcg_fixed, fixed);
				}

				for (int i = 1; i <= mi; ++i)
				{
					#pragma omp critical
					{
						mpi.Send(partition_id, tag_dcg_K_count, K.Count(i));
						mpi.Send(partition_id, tag_dcg_K_index, K.Count(i), &K.index[i][1]);
						mpi.Send(partition_id, tag_dcg_K_entry, K.Count(i), &K.entry[i][1]);
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void InitializeSolver(const Vector<bool>& fixed) throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				fixed_on_boundary.Resize(fixed.size); /**/
				fixed_on_boundary.Fill(false); /**/
				Vector<int> count_fixed_on_boundary(fixed.size); /**/
				count_fixed_on_boundary.Fill(0); /**/

				f_adjust.Fill(0);
				M.Fill(0);

				// Complete f and M
				for (int p = 1; p <= partitions_count; ++p)
				{
					Vector<int>& node_index = partitions.entry[p].node_index;

					int mi = node_index.size*degrees_of_freedom;

					// Receive fi_adjust, Mi
					Vector<T> fi_adjust(mi);
					Vector<T> Mi(mi);
					mpi.Receive(p, tag_dcg_f_adjust, fi_adjust);
					mpi.Receive(p, tag_dcg_M, Mi);
					int max_i = node_index.size;
					for (int i = 1; i <= max_i; ++i)
					{
						int bi = (i - 1)*degrees_of_freedom;
						int n = node_index.entry[i];
						int bn = (n - 1)*degrees_of_freedom;

						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							int bnd = bn + d;
							int bid = bi + d;
							f_adjust.entry[bnd] += fi_adjust.entry[bid];
							M.entry[bnd] += Mi.entry[bid];

							if (fixed.entry[bnd])
							{
								++count_fixed_on_boundary.entry[bnd];
							}
						}
					}
				}
				for (int i = 1; i <= M.size; ++i)
				{
					if (fixed.entry[i])
					{
						M.entry[i] = 1;
						if (count_fixed_on_boundary.entry[i] > 1) /**/
						{
							fixed_on_boundary.entry[i] = true; /**/
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void ComplesateFixed(const Vector<T>& u, Vector<T>& f, const Vector<bool>& fixed, int threads) throw()
		{
			#pragma omp parallel for default(shared) num_threads(threads)
			for (int i = 1; i <= f.size; ++i)
			{
				if (!fixed.entry[i])
				{
					f.entry[i] -= f_adjust.entry[i];
				}
				else
				{
					f.entry[i] = u.entry[i];
				}
			}
		}


		void MultiplyK(const Vector<T>& v, Vector<T>& w, int threads) throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				Vector<bool> check_fixed_on_boundary(fixed_on_boundary.size); /**/
				check_fixed_on_boundary.Fill(false); /**/

				// Fill v
				bool end = false;
				for (int p = 1; p <= partitions_count; ++p)
				{
					Vector<int>& node_index = partitions.entry[p].node_index;

					Vector<T>& vi = partition_vi.entry[p];

					int max_i = node_index.size;
					#pragma omp parallel for default(shared) num_threads(threads)
					for (int i = 1; i <= max_i; ++i)
					{
						int bi = (i - 1)*degrees_of_freedom;
						int n = node_index.entry[i];
						int bn = (n - 1)*degrees_of_freedom;

						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							vi.entry[bi + d] = v.entry[bn + d];
						}
					}

					// Send non stop condition
					mpi.Send(p, tag_dcg_end, end);

					// Send vi
					mpi.Send(p, tag_dcg_vi, vi);

					// Request to receive wi with result
					Vector<T>& wi = partition_wi.entry[p];
					mpi.RequestReceive(p, tag_dcg_wi, wi, request_wi.entry[p]);
				}

				w.Fill(0);
				for (int q = 1; q <= partitions_count; ++q)
				{
					int p = mpi.WaitAnyRequest(request_wi);

					Vector<int>& node_index = partitions.entry[p].node_index;

					Vector<T>& wi = partition_wi.entry[p];

					int max_i = node_index.size;
					#pragma omp parallel for default(shared) num_threads(threads)
					for (int i = 1; i <= max_i; ++i)
					{
						int bi = (i - 1)*degrees_of_freedom;
						int n = node_index.entry[i];
						int bn = (n - 1)*degrees_of_freedom;

						for (int d = 1; d <= degrees_of_freedom; ++d)
						{
							int ii = bn + d;
							if (fixed_on_boundary.entry[ii]) /**/
							{
								if (check_fixed_on_boundary.entry[ii]) /**/
								{
									continue; /**/
								}
								check_fixed_on_boundary.entry[ii] = true; /**/
							}
							w.entry[bn + d] += wi.entry[bi + d];
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		// Conjugate gradient width Jacobi preconditioner
		template <typename U>
		void Solver(Vector<T>& u, Vector<T>& f, T tolerance, int max_steps, int threads) throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				Vector<T> g(m); // Gradient
				Vector<T> p(m); // Descent direcction
				Vector<T> w(m); // w = A*p
				Vector<T> q(m); // Preconditioner system solution

				MultiplyK(u, g, threads);
				U gg = 0.0;
				U gq = 0.0;
				#pragma omp parallel for default(shared) reduction(+:gg,gq) num_threads(threads)
				for (int i = 1; i <= m; ++i)
				{
					g.entry[i] -= f.entry[i]; // g = AX - fB
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

					MultiplyK(p, w, threads);
					U pw = 0.0;
					#pragma omp parallel for default(shared) reduction(+:pw) num_threads(threads)
					for (int i = 1; i <= m; ++i)
					{
						pw += p.entry[i]*w.entry[i]; // pw = p'*w
					}

					U alpha = gq/pw; // alpha = (g'*q)/(p'*w)

					U gngn = 0.0;
					U gnqn = 0.0;
					#pragma omp parallel for default(shared) reduction(+:gngn,gnqn) num_threads(threads)
					for (int i = 1; i <= m; ++i)
					{
						u.entry[i] += (T)(alpha*p.entry[i]); // Xn = u + alpha*p
						g.entry[i] += (T)(alpha*w.entry[i]); // Gn = g + alpha*w
						q.entry[i] = g.entry[i]/M.entry[i]; // Solve for q: MQ = g

						gngn += g.entry[i]*g.entry[i]; // gngn = gn'*gn
						gnqn += g.entry[i]*q.entry[i]; // gnqn = gn'*qn
					}

					U beta = gnqn/gq; // beta = (Gn'*Qn)/(g'*Q)

					#pragma omp parallel for default(shared) num_threads(threads)
					for (int i = 1; i <= m; ++i)
					{
						p.entry[i] = (T)(beta*p.entry[i] - q.entry[i]); // Pn = -q + beta*p
					}

					gg = gngn;
					gq = gnqn;

					Log(2, "%5i  %.5e", step, (double)gg);
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
					mpi.Send(p, tag_dcg_end, end);
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
					mpi.Receive(p, tag_dcg_memory_peak_usage, partition_memory_peak_usage.entry[p]);
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


	private:

		DistributedConjugateGradientMaster& operator = (DistributedConjugateGradientMaster&)
		{
			return *this;
		}
};


template <typename T>
class DistributedConjugateGradientSlave
{
	public:

		const MPI& mpi;

		int m;

		CSRMatrix<T> K;
		Vector<T> u;
		Vector<T> f;
		Vector<bool> fixed;
		Vector<T> f_adjust;

		int threads;


		DistributedConjugateGradientSlave(const MPI& mpi, const int threads) throw(MPI::Exception)
		:	mpi(mpi),
			m(),
			K(),
			u(),
			f(),
			fixed(),
			f_adjust(),
			threads(threads)
		{
			try
			{
				// Receive system size
				mpi.Receive(0, tag_dcg_m, m);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void ReceiveSystemOfEquations() throw(Memory::Exception, MPI::Exception)
		{
			try
			{
				// Receive vectors
				u.Resize(m);
				f.Resize(m);
				fixed.Resize(m);
				mpi.Receive(0, tag_dcg_u, u);
				mpi.Receive(0, tag_dcg_f, f);
				mpi.Receive(0, tag_dcg_fixed, fixed);

				// Receive K
				K.Resize(m, m);
				for (int i = 1; i <= m; ++i)
				{
					int k_max;
					mpi.Receive(0, tag_dcg_K_count, k_max);
					K.AllocateRow(i, k_max);
					mpi.Receive(0, tag_dcg_K_index, k_max, &K.index[i][1]);
					mpi.Receive(0, tag_dcg_K_entry, k_max, &K.entry[i][1]);
				}

				// Apply fixed conditions compensate in f_adjust
				Vector<int> fixed_index(m);
				f_adjust.Resize(m);
				f_adjust.Fill(0);
				for (int i = 1; i <= m; ++i)
				{
					if (fixed.entry[i])
					{
						K.AllocateRow(i, 1);
						K.index[i][1] = i;
						K.entry[i][1] = 1.0;
					}
					else
					{
						int fixed_count = 0;
						int k_max = K.Count(i);
						for (register int k = 1; k <= k_max; ++k)
						{
							register int j = K.index[i][k];
							if (fixed.entry[j])
							{
								f_adjust.entry[i] += K.entry[i][k]*u.entry[j];
								fixed_index.entry[++fixed_count] = j;
							}
						}
						for (register int f = 1; f <= fixed_count; ++f)
						{
							K.RemoveEntry(i, fixed_index.entry[f]);
						}
					}
				}

				// Send contribution for fB_adjust
				mpi.Send(0, tag_dcg_f_adjust, f_adjust); 

				// Send preconditioner
				Vector<T> M(m);
				for (int i = 1; i <= m; ++i)
				{
					M.entry[i] = K(i, i);
				}
				mpi.Send(0, tag_dcg_M, M);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		void MultiplyK(Vector<T>& v, Vector<T>& w) throw(MPI::Exception)
		{
			try
			{
				// Receive v
				mpi.Receive(0, tag_dcg_vi, v);

				// Multiply Kv
				#pragma omp parallel for default(shared) schedule(guided) num_threads(threads)
				for (int i = 1; i <= m; ++i)
				{
					int* __restrict K_index_i = K.index[i];
					T* __restrict K_entry_i = K.entry[i];

					T sum = 0.0;
					int k_max = K.Count(i);
					for (register int k = 1; k <= k_max; ++k)
					{
						sum += K_entry_i[k]*v.entry[K_index_i[k]];
					}
					w.entry[i] = sum;
				}

				// Send w
				mpi.Send(0, tag_dcg_wi, w);
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
				Vector<T> v(m);
				Vector<T> w(m);

				// Multiplications
				for ( ; ; )
				{
					bool end;
					mpi.Receive(0, tag_dcg_end, end);
					if (end)
					{
						break;
					}
					MultiplyK(v, w);
				}
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
				mpi.Send(0, tag_dcg_memory_peak_usage, Memory::peak_usage);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


	private:

		DistributedConjugateGradientSlave& operator = (DistributedConjugateGradientSlave&) throw()
		{
			return *this;
		}
};

#endif
