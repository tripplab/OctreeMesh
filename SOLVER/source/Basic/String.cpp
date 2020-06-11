// String.cpp
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

#include <Basic/String.h>


String::String(const String& string) throw(Memory::Exception)
:	size(string.size),
	data(new char[size])
{
	register const char* __restrict src = string.data;
	register char* __restrict dst = data;
	for (register int i = size; i; --i)
	{
		*(dst++) = *(src++);
	}
}


String::String(int size) throw(Memory::Exception)
:	size(size),
	data(new char[size])
{
	if (!data)
	{
		Throw(Memory::exception);
	}
}


String::String(const char* characters) throw(Memory::Exception)
:	size(Length(characters) + 1),
	data(new char[size])
{
	if (!data)
	{
		Throw(Memory::exception);
	}

	register const char* __restrict src = characters;
	register char* __restrict dst = data;
	for (register int i = size; i; --i)
	{
		*(dst++) = *(src++);
	}
}


String::String(const char* characters, const char*, ...) throw(Memory::Exception)
:	size(0),
	data((char*)0)
{
	for (const char** c = &characters; *c; ++c)
	{
		size += Length(*c);
	}
	++size;
	data = new char[size];
	if (!data)
	{
		Throw(Memory::exception);
	}
	register char* __restrict dst = data;
	for (const char** c = &characters; *c; ++c)
	{
		register const char* __restrict src = *c;
		while (*src)
		{
			*(dst++) = *(src++);
		}
	}
	*dst = 0;
}


String::String(const char* characters, int length) throw(Memory::Exception)
:	size(length + 1),
	data(new char[size])
{
	if (!data)
	{
		Throw(Memory::exception);
	}

	register const char* __restrict src = characters;
	register char* __restrict dst = data;
	for (register int i = length; i && *src; --i)
	{
		*(dst++) = *(src++);
	}
	*dst = 0;
}


String& String::operator = (const String& string) throw(Memory::Exception)
{
	if (this != &string)
	{
		if (size != string.size)
		{
			char* new_data = new char[size];
			if (!new_data)
			{
				Throw(Memory::exception);
			}
			delete [] data;
			size = string.size;
			data = new_data;
		}
		register const char* __restrict src = string.data;
		register char* __restrict dst = data;
		for (register int i = size; i; --i)
		{
			*(dst++) = *(src++);
		}
	}
	return *this;
}


String& String::operator = (const char* characters) throw(Memory::Exception)
{
	Assert(characters);

	int new_size = Length(characters) + 1;
	if (size != new_size)
	{
		char* new_data = new char[new_size];
		if (!new_data)
		{
			Throw(Memory::exception);
		}
		delete [] data;
		size = new_size;
		data = new_data;
	}
	register const char* __restrict src = characters;
	register char* __restrict dst = data;
	for (register int i = size; i; --i)
	{
		*(dst++) = *(src++);
	}
	return *this;
}


int String::Length() const throw()
{
	if (data)
	{
		register const char* __restrict c = data;
		while (*c)
		{
			++c;
		}
		return (int)(c - data);
	}
	return 0;
}


int String::Length(const char* characters) throw()
{
	Assert(characters);

	register const char* __restrict c = characters;
	while (*c)
	{
		++c;
	}
	return (int)(c - characters);
}
