// String.h
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

#ifndef _String_h_
#define _String_h_

#include <Basic/Assert.h>
#include <Basic/Memory.h>


class String
{
	public:
	
		int size;
		char* data;
	

		inline String() throw()
		:	size(0),
			data((char*)0)
		{
		}


		String(const String& string) throw(Memory::Exception);


		String(int size) throw(Memory::Exception);
	

		String(const char* characters) throw(Memory::Exception);


		String(const char* characters, const char*, ...) throw(Memory::Exception);

	
		String(const char* characters, int length) throw(Memory::Exception);
		
	
		inline ~String() throw()
		{
			delete [] data;
		}
		
	
		String& operator = (const String& string) throw(Memory::Exception);
		
	
		String& operator = (const char* characters) throw(Memory::Exception);


		int Length() const throw();


		static int Length(const char* characters) throw();
};

#endif
