// Integer.h
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

#ifndef _Integer_h_
#define _Integer_h_

#include <Basic/Format.h>
#include <Basic/Memory.h>


#if defined _MSC_VER

	typedef signed __int8 sint8;
	typedef unsigned __int8 uint8;

	typedef signed __int16 sint16;
	typedef unsigned __int16 uint16;

	typedef signed __int32 sint32;
	typedef unsigned __int32 uint32;

	typedef signed __int64 sint64;
	typedef unsigned __int64 uint64;

#elif defined __INTEL_COMPILER

	typedef signed char sint8;
	typedef unsigned char uint8;

	typedef signed short sint16;
	typedef unsigned short uint16;

	typedef signed int sint32;
	typedef unsigned int uint32;

	#if defined __i386__
		typedef signed long long sint64;
		typedef unsigned long long uint64;
	#endif
	#if defined __x86_64__
		typedef signed long sint64;
		typedef unsigned long uint64;
	#endif
	#pragma warning(disable: 1572)

#elif defined __GNUG__

	typedef signed char sint8;
	typedef unsigned char uint8;

	#if (__SIZEOF_SHORT__ == 2)
		typedef signed short sint16;
		typedef unsigned short uint16;
	#elif (__SIZEOF_INT__ == 2)
		typedef signed int sint16;
		typedef unsigned int uint16;
	#endif

	#if (__SIZEOF_INT__ == 4)
		typedef signed int sint32;
		typedef unsigned int uint32;
	#elif (__SIZEOF_LONG__ == 4)
		typedef signed long sint32;
		typedef unsigned long uint32;
	#endif

	#if (__SIZEOF_LONG__ == 8)
		typedef signed long sint64;
		typedef unsigned long uint64;
	#elif (__SIZEOF_LONG_LONG__ == 8)
		typedef signed long long sint64;
		typedef unsigned long long uint64;
	#endif

#endif


template <typename T>
class Integer
{
	public:
	
		static const T maximum;
		static const T minimum;
		static const bool is_signed;
		static const int size;
};


class FormatInteger : public Format
{
	public:

		enum Type
		{
			decimal,
			hexadecimal,
			HEXADECIMAL,
			octal
		};

		FormatInteger(bool space, bool sign, int width, Type type) throw(Memory::Exception);
};

#endif
