// Random.h
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

#ifndef _Random_h_
#define _Random_h_

#include <Basic/Integer.h>


enum RandomType
{
	random_linearcongruential32,
	random_mersennetwister32,
	random_motherofall32,
	random_lecuyer32
};


template <RandomType type>
class Random;


// http://en.wikipedia.org/wiki/Linear_congruential_generator
template <>
class Random<random_linearcongruential32>
{
	public:

		static const uint32 maximum;

		Random(uint32 seed = 1) throw();

		uint32 Get() throw();


	protected:

		uint32 last;
};


// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura, freely usable
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
template <>
class Random<random_mersennetwister32>
{
	public:

		static const uint32 maximum;

		Random(uint32 seed = 1) throw();

		uint32 Get() throw();


	protected:

		int index;
		uint32 x[624];
};


// Copyright (C) 1997 - 2007, Agner Fog, GNU General Public License
// http://www.agner.org/random
template <>
class Random<random_motherofall32>
{
	public:

		static const uint32 maximum;

		Random(uint32 seed = 1) throw();

		uint32 Get() throw();


	protected:

		uint32 x0;
		uint32 x1;
		uint32 x2;
		uint32 x3;
		uint32 x4;

		uint32 Init(uint32 seed) throw();
};


// P. L'Ecuyer, "Maximally Equidistributed Combined Tausworthe Generators", Mathematics of Computation, 65, 213 (1996), 203–213.
template <>
class Random<random_lecuyer32>
{
	public:

		static const uint32 maximum;

		Random(uint32 seed = 1) throw();

		uint32 Get() throw();


	protected:

		uint32 s1;
		uint32 s2;
		uint32 s3;

		uint32 Init(uint32 seed) throw();
};

#endif
