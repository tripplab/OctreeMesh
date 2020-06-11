// TriangularSystem.h
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

#ifndef _TriangularSystem_h_
#define _TriangularSystem_h_

#include <Container/CSRMatrix.h>
#include <Container/TriangularMatrix.h>
#include <Container/Vector.h>


template <typename T>
void LowerTriangularSystem(const CSRMatrix<T>& L, Vector<T>& c, const Vector<T>& b) throw()
{
	int n = L.rows;
	for (int i = 1; i <= n; ++i)
	{
		int* __restrict L_index_i = L.index[i];
		T* __restrict L_entry_i = L.entry[i];

		T sum = b.entry[i];
		int k_max = L.Count(i);
		for (register int k = 1; k < k_max; ++k)
		{
			register int j = L_index_i[k];
			sum -= L_entry_i[k]*c.entry[j];
		}
		c.entry[i] = sum/L_entry_i[k_max];
	}
}


template <typename T>
void UpperTriangularSystem(const CSRMatrix<T>& U, Vector<T>& x, const Vector<T>& c) throw()
{
	int n = U.rows;
	for (int i = n; i >= 1; --i)
	{
		int* __restrict U_index_i = U.index[i];
		T* __restrict U_entry_i = U.entry[i];

		T sum = c.entry[i];
		int k_max = U.Count(i);
		for (register int k = 2; k <= k_max; ++k)
		{
			register int j = U_index_i[k];
			sum -= U_entry_i[k]*x.entry[j];
		}
		x.entry[i] = sum/U_entry_i[1];
	}
}

#endif
