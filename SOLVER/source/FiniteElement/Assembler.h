// Assembler.h
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

#ifndef _Assembler_h_
#define _Assembler_h_

#include <Basic/Memory.h>
#include <Container/CSRMatrix.h>
#include <Container/Matrix.h>
#include <Container/Set.h>
#include <Container/Vector.h>
#include <FiniteElement/Mesh.h>


class Assembler
{
	public:

		const Mesh& mesh;

		const Vector<int>& element_index;

		const Vector<int>& node_index;

		const int degrees_of_freedom;

		Vector<int> global_to_local;


	public:

		Assembler(const Mesh& mesh, const Vector<int>& element_index, const Vector<int>& node_index, const int degrees_of_freedom) throw(Memory::Exception);


		template <typename T>
		void AllocateMatrix(CSRMatrix<T>& A) const throw(Memory::Exception)
		{
			Assert(A.rows == node_index.size*degrees_of_freedom);
			Assert(A.rows == A.columns);

			try
			{
				// Local adjacency for nodes
				Vector<Set<int> > adjacency(node_index.size);
				for (int i = 1; i <= element_index.size; ++i)
				{
					int e = element_index.entry[i];
					for (int c = 1; c < mesh.nodes_per_element; ++c)
					{
						int n = global_to_local.entry[mesh.connectivity.entry[e][c]];
						if (n != 0)
						{
							for (register int cr = c + 1; cr <= mesh.nodes_per_element; ++cr)
							{
								int nr = global_to_local.entry[mesh.connectivity.entry[e][cr]];
								if (nr != 0)
								{
									adjacency.entry[n].Append(nr);
									adjacency.entry[nr].Append(n);
								}
							}
						}
					}
				}

				// Generate sparse matrix
				for (int n = 1; n <= node_index.size; ++n)
				{
					adjacency.entry[n].Append(n);
					int adjacency_n_size = adjacency.entry[n].size;
					int row_size = adjacency_n_size*degrees_of_freedom;
					for (int d = 1; d <= degrees_of_freedom; ++d)
					{
						int i = (n - 1)*degrees_of_freedom + d;
						A.AllocateRow(i, row_size);
						int k = 1;
						for (register SetItem<int>* __restrict set_item = adjacency.entry[n].first; set_item; set_item = set_item->next)
						{
							register int rn = set_item->value;
							for (register int rd = 1; rd <= degrees_of_freedom; ++rd, ++k)
							{
								register int j = (rn - 1)*degrees_of_freedom + rd;
								A.index[i][k] = j;
							}
						}
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		template <typename T>
		void AssembleAe(const int element_id, const Matrix<T>& Ae, CSRMatrix<T>& A) const throw()
		{
			Assert(A.rows == node_index.size*degrees_of_freedom);
			Assert(A.rows == A.columns);
			Assert(Ae.rows == mesh.nodes_per_element*degrees_of_freedom);
			Assert(Ae.rows == Ae.columns);

			for (int ni = 1; ni <= mesh.nodes_per_element; ++ni)
			{
				int li = global_to_local.entry[mesh.connectivity.entry[element_id][ni]];
				if (li != 0)
				{
					for (int nj = 1; nj <= mesh.nodes_per_element; ++nj)
					{
						int lj = global_to_local.entry[mesh.connectivity.entry[element_id][nj]];
						if (lj != 0)
						{
							for (int di = 1; di <= degrees_of_freedom; ++di)
							{
								int i = (li - 1)*degrees_of_freedom + di;
								for (int dj = 1; dj <= degrees_of_freedom; ++dj)
								{
									int j = (lj - 1)*degrees_of_freedom + dj;
									A(i, j) += Ae.entry[(ni - 1)*degrees_of_freedom  + di][(nj - 1)*degrees_of_freedom  + dj];
								}
							}
						}
					}
				}
			}
		}


		template <typename T>
		void AssembleAe(const int element_id, const Matrix<T>& Ae, Matrix<T>& A) const throw()
		{
			Assert(A.rows == node_index.size*degrees_of_freedom);
			Assert(A.rows == A.columns);
			Assert(Ae.rows == mesh.nodes_per_element*degrees_of_freedom);
			Assert(Ae.rows == Ae.columns);

			for (int ni = 1; ni <= mesh.nodes_per_element; ++ni)
			{
				int li = global_to_local.entry[mesh.connectivity.entry[element_id][ni]];
				if (li != 0)
				{
					for (int nj = 1; nj <= mesh.nodes_per_element; ++nj)
					{
						int lj = global_to_local.entry[mesh.connectivity.entry[element_id][nj]];
						if (lj != 0)
						{
							for (int di = 1; di <= degrees_of_freedom; ++di)
							{
								int i = (li - 1)*degrees_of_freedom + di;
								for (int dj = 1; dj <= degrees_of_freedom; ++dj)
								{
									int j = (lj - 1)*degrees_of_freedom + dj;
									A.entry[i][j] += Ae.entry[(ni - 1)*degrees_of_freedom  + di][(nj - 1)*degrees_of_freedom  + dj];
								}
							}
						}
					}
				}
			}
		}


		template <typename T>
		void AssembleV(const Vector<T>& global_v, Vector<T>& v) const throw(Memory::Exception)
		{
			Assert(global_v.size == mesh.nodes_count*degrees_of_freedom);
			Assert(v.size == node_index.size*degrees_of_freedom);

			for (int l = 1; l <= node_index.size; ++l)
			{
				register int g = node_index.entry[l];
				for (register int d = 1; d <= degrees_of_freedom; ++d)
				{
					register int ig = (g - 1)*degrees_of_freedom + d;
					register int il = (l - 1)*degrees_of_freedom + d;
					v.entry[il] = global_v.entry[ig];
				}
			}
		}


	private:

		inline Assembler& operator = (const Assembler&) throw();

};

#endif
