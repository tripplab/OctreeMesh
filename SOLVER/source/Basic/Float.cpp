// Float.cpp
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
#include <Basic/Float.h>
#include <Basic/Integer.h>
#include <stdio.h>


#if defined _MSC_VER

	#define FLOAT_MAXIMUM 3.402823466e+38F
	#define FLOAT_MINIMUM -3.402823466e+38F
	#define FLOAT_EPSILON 1.192092896e-07F
	#define FLOAT_SIZE 4

	#define DOUBLE_MAXIMUM 1.7976931348623158e+308
	#define DOUBLE_MINIMUM -1.7976931348623158e+308
	#define DOUBLE_EPSILON 2.2204460492503131e-016
	#define DOUBLE_SIZE 8

	#define LDOUBLE_MAXIMUM 1.7976931348623158e+308
	#define LDOUBLE_MINIMUM -1.7976931348623158e+308
	#define LDOUBLE_EPSILON 2.2204460492503131e-016
	#define LDOUBLE_SIZE 8

#elif defined __GNUG__

	#define FLOAT_MAXIMUM __FLT_MAX__
	#define FLOAT_MINIMUM __FLT_MIN__
	#define FLOAT_EPSILON __FLT_EPSILON__
	#define FLOAT_SIZE __FLT_SIZE__

	#define DOUBLE_MAXIMUM __DBL_MAX__
	#define DOUBLE_MINIMUM __DBL_MIN__
	#define DOUBLE_EPSILON __DBL_EPSILON__
	#define DOUBLE_SIZE __DBL_SIZE__

	#define LDOUBLE_MAXIMUM __LDBL_MAX__
	#define LDOUBLE_MINIMUM __LDBL_MIN__
	#define LDOUBLE_EPSILON __LDBL_EPSILON__
	#define LDOUBLE_SIZE __LDBL_SIZE__

#endif


union FloatBits
{
	uint32 b;
	float f;
};


static const FloatBits float_bits_infinity = {0x7F800000};


template <>
class Float<float>
{
	public:
	
		static const float maximum;
		static const float minimum;
		static const float epsilon;
		static const float infinite;
		static const unsigned int size;


		static bool IsNaN(float value) throw();
};


const float Float<float>::maximum = FLOAT_MAXIMUM;
const float Float<float>::minimum = FLOAT_MINIMUM;
const float Float<float>::epsilon = FLOAT_EPSILON;
const float Float<float>::infinite = float_bits_infinity.f;
const unsigned int Float<float>::size = sizeof(float);


bool Float<float>::IsNaN(float value) throw()
{
	static const uint32 mask1 = 0x7F800000;
	static const uint32 mask2 = 0x007FFFFF;

	FloatBits m;
	m.f = value;
	return ((m.b & mask1) == mask1) && ((m.b & mask2) != 0);
}


union DoubleBits
{
	uint64 b;
	double f;
};


static const DoubleBits double_bits_infinity = {0x7FF0000000000000ULL};


template <>
class Float<double>
{
	public:
	
		static const double maximum;
		static const double minimum;
		static const double epsilon;
		static const double infinite;
		static const unsigned int size;


		static bool IsNaN(double value) throw();
};


const double Float<double>::maximum = DOUBLE_MAXIMUM;
const double Float<double>::minimum = DOUBLE_MINIMUM;
const double Float<double>::epsilon = DOUBLE_EPSILON;
const double Float<double>::infinite = double_bits_infinity.f;
const unsigned int Float<double>::size = sizeof(double);


bool Float<double>::IsNaN(double value) throw()
{
	static const uint64 mask1 = 0x7FE0000000000000ULL;
	static const uint64 mask2 = 0x000FFFFFFFFFFFFFULL;

	DoubleBits m;
	m.f = value;
	return ((m.b & mask1) == mask1) && ((m.b & mask2) != 0);
}


FormatFloat::FormatFloat(bool space, bool sign, int width, int precision, Type type) throw(Memory::Exception)
:	Format()
{
	Assert(width > 0);
	Assert(precision >= 0);

	int size = 4;
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
	if (precision == 0)
	{
		++size;
	}
	else
	{
		for (register int p = precision; p; ++size)
		{
			p /= 10;
		}
	}
	*(char**)&format_string = new char[size];
	if (!format_string )
	{
		Throw(Memory::exception);
	}
	sprintf(*(char**)&format_string, "%%%s%s%i.%i%c", space ? " " : "", sign ? "+" : "", width, precision, (type == fixed) ? 'f' : (type == exponential) ? 'e' : 'g');
}
