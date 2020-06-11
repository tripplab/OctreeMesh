// Float.h
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

#ifndef _Float_h_
#define _Float_h_

#include <Basic/Format.h>
#include <Basic/Memory.h>


template <typename T>
class Float
{
	public:
	
		static const T maximum;
		static const T minimum;
		static const T epsilon;
		static const T infinite;
		static const unsigned int size;

		static bool IsNaN(T value) throw();
};


class FormatFloat : public Format
{
	public:

		enum Type
		{
			fixed,
			exponential,
			automatic
		};

		FormatFloat(bool space, bool sign, int width, int precision, Type type) throw(Memory::Exception);
};

#endif
