// Nodes.h
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

#ifndef _Nodes_h_
#define _Nodes_h_

#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Container/Matrix.h>


template <typename T>
class Nodes
{
	public:

		const int& dimension;

		const int& nodes_count;

		Matrix<T> coordinate;


		Nodes() throw()
		:	dimension(coordinate.columns),
			nodes_count(coordinate.rows),
			coordinate()
		{
		}


		Nodes(const Nodes& nodes) throw(Memory::Exception)
		:	dimension(coordinate.columns),
			nodes_count(coordinate.rows),
			coordinate(nodes.coordinate)
		{
		}


		Nodes(int dimension, int nodes_count) throw(Memory::Exception)
		:	dimension(coordinate.columns),
			nodes_count(coordinate.rows),
			coordinate(nodes_count, dimension)
		{
		}


		Nodes& operator = (const Nodes& nodes) throw()
		{
			if (&nodes != this)
			{
				coordinate = nodes.coordinate;
			}
			return *this;
		}


		void PrintInfo(void) throw()
		{
			Log(1, "Nodes ----------------------------------------------------------------");
			Log(1, "-Dimension:         %i", dimension);
			Log(1, "-Nodes:             %i", nodes_count);
		}
};

#endif
