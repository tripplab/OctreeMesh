// Random.cpp
// Copyright (C) 2012 Miguel Vargas (miguelvargas@users.sourceforge.net)
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

#include <Basic/Random.h>


const uint32 Random<random_linearcongruential32>::maximum = Integer<uint32>::maximum;


Random<random_linearcongruential32>::Random(uint32 seed) throw()
:	last(seed)
{
}


uint32 Random<random_linearcongruential32>::Get() throw()
{
	return (last = (uint32)1103515245UL*last + (uint32)12345UL);
}


const uint32 Random<random_mersennetwister32>::maximum = Integer<uint32>::maximum;


Random<random_mersennetwister32>::Random(uint32 seed) throw()
:	index(624)
{
	x[0] = seed;
	register uint32* element = x;
	for (register int i = 1; i < 624; ++i)
	{
		uint32 e = (uint32)1812433253UL*(*element ^ (*element >> 30)) + (uint32)i;
		*(++element) = e;
	}
}


uint32 Random<random_mersennetwister32>::Get() throw()
{
	register uint32 y;
	if (index == 624)
	{
		static const uint32 a[2] = {0UL, (uint32)2567483615UL};

		for (register int k = 0; k < 227; ++k)
		{
			y = (x[k] & (uint32)0x80000000UL) | (x[k + 1] & (uint32)0x7fffffffUL);
			x[k] = x[k + 397] ^ (y >> 1) ^ a[y & 1UL];
		}
		for (register int k = 227; k < 623; ++k)
		{
			y = (x[k] & (uint32)0x80000000UL) | (x[k + 1] & (uint32)0x7fffffffUL);
			x[k] = x[k - 227] ^ (y >> 1) ^ a[y & 1UL];
		}
		y = (x[623] & (uint32)0x80000000UL) | (x[0] & (uint32)0x7fffffffUL);
		x[623] = x[396] ^ (y >> 1) ^ a[y & 1UL];
		index = 0;
	}
	y = x[index++];
	y ^= (y >> 11);
	y ^= (y << 7) & (uint32)2636928640UL;
	y ^= (y << 15) & (uint32)4022730752UL;
	y ^= (y >> 18);
	return y;
}


const uint32 Random<random_motherofall32>::maximum = Integer<uint32>::maximum;


Random<random_motherofall32>::Random(uint32 seed) throw()
:	x0(Init(seed)),
	x1(Init(x0)),
	x2(Init(x1)),
	x3(Init(x2)),
	x4(Init(x3))
{
	for (register int i = 19; i; --i)
	{
		Get();
	}
}


uint32 Random<random_motherofall32>::Get() throw()
{
	register uint64 sum =
		(uint64)2111111111UL*(uint64)x3 +
		(uint64)1492UL*(uint64)x2 +
		(uint64)1776UL*(uint64)x1 +
		(uint64)5115UL*(uint64)x0 +
		(uint64)x4;
	x3 = x2;
	x2 = x1;
	x1 = x0;
	x4 = (uint32)(sum >> 32);
	x0 = (uint32)sum;
	return x0;
}


uint32 Random<random_motherofall32>::Init(uint32 seed) throw()
{
	return seed*(uint32)29943829UL - (uint32)1UL;
}


const uint32 Random<random_lecuyer32>::maximum = Integer<uint32>::maximum;


Random<random_lecuyer32>::Random(uint32 seed) throw()
:	s1(Init(seed)),
	s2(Init(s1)),
	s3(Init(s2))
{
	for (register int i = 19; i; --i)
	{
		Get();
	}
}


uint32 Random<random_lecuyer32>::Get() throw()
{
	register uint32 b;
	b = ((s1 << 13) ^ s1) >> 19;
	s1 = ((s1 & (uint32)4294967294UL) << 12) ^ b;
	b = ((s2 << 2) ^ s2) >> 25;
	s2 = ((s2 & (uint32)4294967288UL) << 4) ^ b;
	b = ((s3 << 3) ^ s3) >> 11;
	s3 = ((s3 & (uint32)4294967280UL) << 17) ^ b;
	return s1 ^ s2 ^ s3;
}


uint32 Random<random_lecuyer32>::Init(uint32 seed) throw()
{
	return seed*(uint32)29943829UL - (uint32)1UL;
}
