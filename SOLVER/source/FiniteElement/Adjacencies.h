// Adjacencies.h
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

#ifndef _Adjacencies_h
#define _Adjacencies_h

#include <Basic/Assert.h>
#include <Container/Vector.h>


struct Adjacencies
{
	Vector<int> start;
	Vector<int> data;


	inline int Get(int element_id, int index) const throw()
	{
		Assert(element_id < start.size);

		return data.entry[start.entry[element_id] - 1 + index];
	}


	inline int Size(int element_id) const throw()
	{
		Assert(element_id < start.size);

		return start.entry[element_id + 1] - start.entry[element_id];
	}

};

#endif
