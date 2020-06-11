// Distribution.h
// Copyright (C) 2007 Miguel Vargas (miguelvargas@users.sourceforge.net)
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

#ifndef _Distribution_h_
#define _Distribution_h_

#include <Basic/Random.h>
#include <math.h>


template <typename T = double, RandomType random_type = random_linearcongruential32, typename U = double>
class Uniform
{
	public:

		const T minimum;
		const T maximum;


		Uniform(T minimum = 0, T maximum = 1, uint32 seed = 1) throw()
		:	minimum(minimum),
			maximum(maximum),
			random(seed),
			factor(((U)maximum - (U)minimum)/((U)random.maximum + (U)1))
		{
		}


		Uniform& operator = (const Uniform& uniform)
		{
			if (this != &uniform)
			{
				*(T*)&minimum = uniform.maximum;
				*(T*)&minimum = uniform.minimum;
				random = uniform.random;
			}
			return *this;
		}


		T Get() throw()
		{
			return minimum + (T)(factor*random.Get());
		}


	protected:

		Random<random_type> random;
		U factor;
};


// http://en.wikipedia.org/wiki/Box_Muller_transform
template <typename T = double, RandomType random_type = random_linearcongruential32, typename U = double>
class Normal
{
	public:

		U mean;
		U variance;


		Normal(U mean = 0, U variance = 1, uint32 seed = 1) throw()
		:	mean(mean),
			variance(variance),
			random1(seed),
			random2(seed + 1),
			factor1((U)1/random1.maximum),
			factor2((U)6.283185307179586476925286766559L/random2.maximum)
		{
		}


		Normal& operator = (const Normal& normal)
		{
			if (this != &normal)
			{
				*(U*)&mean = normal.mean;
				*(U*)&variance = normal.variance;
				random1 = normal.random1;
				random2 = normal.random2;
			}
			return *this;
		}


		T Get() throw()
		{
			return (T)(sqrt(-(U)2*variance*log(factor1*random1.Get()))*cos(factor2*random2.Get()) + mean);
		}

	protected:

		Random<random_type> random1;
		Random<random_type> random2;
		U factor1;
		U factor2;
};


// http://en.wikipedia.org/wiki/Negative_exponential_distribution
template <typename T = double, RandomType random_type = random_linearcongruential32, typename U = double>
class Exponential
{
	public:

		U rate;


		Exponential(U rate = 1, uint32 seed = 1)
		:	rate(rate),
			random(seed),
			factor1((U)1/random.maximum),
			factor2((U)-1/rate)
		{
		}


		Exponential& operator = (const Exponential& exponential)
		{
			if (this != &exponential)
			{
				*(U*)&rate = exponential.rate;
				random = exponential.random;
			}
			return *this;
		}


		T Get()
		{
			return (T)(log(factor1*random.Get())*factor2);
		}


	protected:

		Random<random_type> random;
		U factor1;
		U factor2;
};

#endif
