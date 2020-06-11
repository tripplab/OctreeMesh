// DenseTriangularSystem.h
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

#ifndef _DenseTriangularSystem_h_
#define _DenseTriangularSystem_h_

#include <Container/TriangularMatrix.h>
#include <Container/Vector.h>


template <typename T>
void LowerTriangularSystem(const LowerTriangularMatrix<T>& L, Vector<T>& c, const Vector<T>& b) throw()
{
	int n = L.rows;
	for (int i = 1; i <= n; ++i)
	{
		T* __restrict L_entry_i = L.entry[i];

		register T sum = b.entry[i];
		for (register int j = 1; j < i; ++j)
		{
			sum -= L_entry_i[j]*c.entry[j];
		}
		c.entry[i] = sum/L_entry_i[i];
	}
}


template <typename T>
void UpperTriangularSystem(const UpperTriangularMatrix<T>& U, Vector<T>& x, const Vector<T>& c) throw()
{
	int n = U.rows;
	for (int i = n; i >= 1; --i)
	{
		T* __restrict U_entry_i = U.entry[i];

		register T sum = c.entry[i];
		for (register int j = i + 1; j <= n; ++j)
		{
			sum -= U_entry_i[j]*x.entry[j];
		}
		x.entry[i] = sum/U_entry_i[i];
	}
}

#endif
