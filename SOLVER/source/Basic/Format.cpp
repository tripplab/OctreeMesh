// Format.cpp
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

#include <Basic/Debug.h>
#include <Basic/Format.h>
#include <string.h>


Format::Format() throw()
:	format_string((char*)0)
{
}


Format::Format(const Format& format) throw(Memory::Exception)
:	format_string((char*)0)
{
	if (format.format_string)
	{
		size_t size = strlen(format.format_string) + 1;
		*(char**)&format_string = new char[size];
		if (!format_string)
		{
			Throw(Memory::exception);
		}
		memcpy(*(char**)&format_string, format.format_string, size);
	}
}


Format::~Format() throw()
{
	delete [] format_string;
}


Format& Format::operator = (const Format& format) throw(Memory::Exception)
{
	if (&format != this)
	{
		if (format.format_string)
		{
			size_t size = strlen(format.format_string) + 1;
			*(char**)&format_string = new char[size];
			if (!format_string)
			{
				Throw(Memory::exception);
			}
			memcpy(*(char**)&format_string, format.format_string, size);
		}
	}
	return *this;
}
