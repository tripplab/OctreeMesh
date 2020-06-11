// Integer.cpp
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

#include <Basic/Assert.h>
#include <Basic/Debug.h>
#include <Basic/Integer.h>
#include <stdio.h>


#if defined _MSC_VER

	#ifdef _CHAR_UNSIGNED
		#define CHAR_SIGNED 0
		#define CHAR_MAXIMUM (char)(255)
		#define CHAR_MINIMUM (char)(0)
	#else
		#define CHAR_SIGNED 1
		#define CHAR_MAXIMUM (char)(127)
		#define CHAR_MINIMUM (char)(-127 -1)
	#endif

	#define SCHAR_MAXIMUM (signed char)(127)
	#define SCHAR_MINIMUM (signed char)(-127 - 1)

	#define UCHAR_MAXIMUM (unsigned char)(255)
	#define UCHAR_MINIMUM (unsigned char)(0)

	#define SSHORT_MAXIMUM (signed short)(32767)
	#define SSHORT_MINIMUM (signed short)(-32767 - 1)

	#define USHORT_MAXIMUM (unsigned short)(65535)
	#define USHORT_MINIMUM (unsigned short)(0)

	#define SINT_MAXIMUM (signed int)(2147483647)
	#define SINT_MINIMUM (signed int)(-2147483647 - 1)

	#define UINT_MAXIMUM (unsigned int)(4294967295)
	#define UINT_MINIMUM (unsigned int)(0)

	#define SLONG_MAXIMUM (signed long)(2147483647)
	#define SLONG_MINIMUM (signed long)(-2147483647 - 1)

	#define ULONG_MAXIMUM (unsigned long)(4294967295)
	#define ULONG_MINIMUM (unsigned long)(0)

	#define SLLONG_MAXIMUM (signed long long)(9223372036854775807)
	#define SLLONG_MINIMUM (signed long long)(-9223372036854775807 - 1)

	#define ULLONG_MAXIMUM (unsigned long long)(18446744073709551615)
	#define ULLONG_MINIMUM (unsigned long long)(0)

#elif defined __GNUG__

	#ifdef __CHAR_UNSIGNED__
		#define CHAR_SIGNED 0
		#define CHAR_MAXIMUM (char)(255)
		#define CHAR_MINIMUM (char)(0)
	#else
		#define CHAR_SIGNED 1
		#define CHAR_MAXIMUM (char)(127)
		#define CHAR_MINIMUM (char)(-128)
	#endif

	#define SCHAR_MAXIMUM (signed char)(127)
	#define SCHAR_MINIMUM (signed char)(-128)

	#define UCHAR_MAXIMUM (unsigned char)(255)
	#define UCHAR_MINIMUM (unsigned char)(0)

	#define SSHORT_MAXIMUM (signed short)(__SHRT_MAX__)
	#define SSHORT_MINIMUM (signed short)(-__SHRT_MAX__ - 1)

	#define USHORT_MAXIMUM (unsigned short)((__SHRT_MAX__ << 1) + 1)
	#define USHORT_MINIMUM (unsigned short)(0)

	#define SINT_MAXIMUM (signed int)(__INT_MAX__)
	#define SINT_MINIMUM (signed int)(-__INT_MAX__ - 1)

	#define UINT_MAXIMUM (unsigned int)((__INT_MAX__ << 1) + 1)
	#define UINT_MINIMUM (unsigned int)(0)

	#define SLONG_MAXIMUM (signed long)(__LONG_MAX__)
	#define SLONG_MINIMUM (signed long)(-__LONG_MAX__ - 1)

	#define ULONG_MAXIMUM (unsigned long)((__LONG_MAX__ << 1) + 1)
	#define ULONG_MINIMUM (unsigned long)(0)

	#define SLLONG_MAXIMUM (signed long long)(__LONG_LONG_MAX__)
	#define SLLONG_MINIMUM (signed long long)(-__LONG_LONG_MAX__ - 1)

	#define ULLONG_MAXIMUM (unsigned long long)((__LONG_LONG_MAX__ << 1) + 1)
	#define ULLONG_MINIMUM (unsigned long long)(0)

#endif


template class Integer<char>;
template <> const char Integer<char>::maximum = CHAR_MAXIMUM;
template <> const char Integer<char>::minimum = CHAR_MINIMUM;
template <> const bool Integer<char>::is_signed = CHAR_SIGNED;
template <> const int Integer<char>::size = sizeof(char);


template class Integer<signed char>;
template <> const signed char Integer<signed char>::maximum = SCHAR_MAXIMUM;
template <> const signed char Integer<signed char>::minimum = SCHAR_MINIMUM;
template <> const bool Integer<signed char>::is_signed = true;
template <> const int Integer<signed char>::size = sizeof(signed char);


template class Integer<unsigned char>;
template <> const unsigned char Integer<unsigned char>::maximum = UCHAR_MAXIMUM;
template <> const unsigned char Integer<unsigned char>::minimum = UCHAR_MINIMUM;
template <> const bool Integer<unsigned char>::is_signed = false;
template <> const int Integer<unsigned char>::size = sizeof(unsigned char);


template class Integer<signed short>;
template <> const signed short Integer<signed short>::maximum = SSHORT_MAXIMUM;
template <> const signed short Integer<signed short>::minimum = SSHORT_MINIMUM;
template <> const bool Integer<signed short>::is_signed = true;
template <> const int Integer<signed short>::size = sizeof(signed short);


template class Integer<unsigned short>;
template <> const unsigned short Integer<unsigned short>::maximum = USHORT_MAXIMUM;
template <> const unsigned short Integer<unsigned short>::minimum = USHORT_MINIMUM;
template <> const bool Integer<unsigned short>::is_signed = false;
template <> const int Integer<unsigned short>::size = sizeof(unsigned short);


template class Integer<signed int>;
template <> const signed int Integer<signed int>::maximum = SINT_MAXIMUM;
template <> const signed int Integer<signed int>::minimum = SINT_MINIMUM;
template <> const bool Integer<signed int>::is_signed = true;
template <> const int Integer<signed int>::size = sizeof(signed int);


template class Integer<unsigned int>;
template <> const unsigned int Integer<unsigned int>::maximum = UINT_MAXIMUM;
template <> const unsigned int Integer<unsigned int>::minimum = UINT_MINIMUM;
template <> const bool Integer<unsigned int>::is_signed = false;
template <> const int Integer<unsigned int>::size = sizeof(unsigned int);


template class Integer<signed long>;
template <> const signed long Integer<signed long>::maximum = SLONG_MAXIMUM;
template <> const signed long Integer<signed long>::minimum = SLONG_MINIMUM;
template <> const bool Integer<signed long>::is_signed = true;
template <> const int Integer<signed long>::size = sizeof(signed long);


template class Integer<unsigned long>;
template <> const unsigned long Integer<unsigned long>::maximum = ULONG_MAXIMUM;
template <> const unsigned long Integer<unsigned long>::minimum = ULONG_MINIMUM;
template <> const bool Integer<unsigned long>::is_signed = false;
template <> const int Integer<unsigned long>::size = sizeof(unsigned long);


template class Integer<signed long long>;
template <> const signed long long Integer<signed long long>::maximum = SLLONG_MAXIMUM;
template <> const signed long long Integer<signed long long>::minimum = SLLONG_MINIMUM;
template <> const bool Integer<signed long long>::is_signed = true;
template <> const int Integer<signed long long>::size = sizeof(signed long long);


template class Integer<unsigned long long>;
template <> const unsigned long long Integer<unsigned long long>::maximum = ULLONG_MAXIMUM;
template <> const unsigned long long Integer<unsigned long long>::minimum = ULLONG_MINIMUM;
template <> const bool Integer<unsigned long long>::is_signed = false;
template <> const int Integer<unsigned long long>::size = sizeof(unsigned long long);


FormatInteger::FormatInteger(bool space, bool sign, int width, Type type) throw(Memory::Exception)
:	Format()
{
	Assert(width > 0);

	int size = 3;
	if (space)
	{
		++size;
	}
	if (sign)
	{
		++size;
	}
	for (register int w = width; w; ++size)
	{
		w /= 10;
	}
	*(char**)&format_string = new char[size];
	if (!format_string )
	{
		Throw(Memory::exception);
	}
	sprintf(*(char**)&format_string, "%%%s%s%i%c", space ? " " : "", sign ? "+" : "", width, (type == decimal) ? 'd' : (type == hexadecimal) ? 'x' : (type == HEXADECIMAL) ? 'X' : 'o');
}
