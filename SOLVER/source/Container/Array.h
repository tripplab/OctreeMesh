// Array.h
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

#ifndef _Array_h_
#define _Array_h_

#include <Basic/Assert.h>


template <typename T, int S>
class Array
{
	public:

		static const int size = S;
		T data[S];
		T* entry;


		Array() throw()
		:	entry(data - 1)
		{
		}


		Array(const Array& a) throw()
		:	entry(data - 1)
		{
			Copy(a.data);
		}


		Array& operator = (const Array& a) throw()
		{
			if (this != &a)
			{
				Copy(a.data);
			}
			return *this;
		}

		
		template <typename U>
		Array& operator = (const Array<U, S>& a) throw()
		{
			register const U* __restrict src = a.data;
			register T* __restrict dst = data;
			for (register int k = S; k; --k, ++dst, ++src)
			{
				*dst = (T)*src;
			}
			return *this;
		}


		inline T& operator () (int k) throw()
		{
			Assert(k >= 1);
			Assert(k <= S);

			return entry[k];
		}


		void Fill(const T& value) throw()
		{
			register T* __restrict dst = data;
			for (register int k = S; k; --k, ++dst)
			{
				*dst = value;
			}
		}


		void Copy(const T* data) throw()
		{
			Assert(data);

			register const T* __restrict src = data;
			register T* __restrict dst = this->data;
			for (register int k = S; k; --k, ++dst, ++src)
			{
				*dst = *src;
			}
		}


		void Sort() throw()
		{
			// Combsort11 http://en.wikipedia.org/wiki/Comb_sort
			bool swapped;
			register int gap = S;
			do
			{
				gap = (gap*10)/13;
				switch (gap)
				{
					case 0:
					{
						gap = 1;
						break;
					}
					case 9:
					case 10:
					{
						gap = 11;
						break;
					}
					default:
					{
						break;
					}
				}
				swapped = false;
				for (register int j = gap + 1; j <= S; ++j)
				{
					int k = j - gap;
					if (entry[j] < entry[k])
					{
						register T t = entry[j];
						entry[j] = entry[k];
						entry[k] = t;
						swapped = true;
					}
				}
			} while (swapped || (gap > 1));
		}
};

#endif
