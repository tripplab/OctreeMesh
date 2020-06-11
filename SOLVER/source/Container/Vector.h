// Vector.h
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

#ifndef _Vector_h_
#define _Vector_h_

#include <Basic/Debug.h>
#include <Basic/Assert.h>
#include <Basic/Memory.h>
#include <Basic/System.h>


template <typename T>
class Vector
{
	public:

		const int size;
		T* const data;
		T* const entry;


	private:

		bool manage;


	public:

		Vector() throw()
		:	size(0),
			data((T*)0),
			entry((T*)0),
			manage(true)
		{
		}


		Vector(const Vector& a) throw(Memory::Exception)
		:	size(a.size),
			data(new T[size]),
			entry(data - 1),
			manage(true)
		{
			if (!data)
			{
				Throw(Memory::exception);
			}
			Copy(a.data);
		}


		Vector(int size) throw(Memory::Exception)
		:	size(size),
			data(new T[size]),
			entry((T*)0),
			manage(true)
		{
			if (!data)
			{
				Throw(Memory::exception);
			}
			*(T**)&entry = data - 1;
		}

		
		Vector(T* data, int size) throw(Memory::Exception)
		:	size(size),
			data(data),
			entry(data - 1),
			manage(false)
		{
			Assert(data);
			Assert(size >= 0);
		}


		~Vector() throw()
		{
			if (manage)
			{
				delete [] data;
			}
		}


		Vector& operator = (const Vector& a) throw()
		{
			if (this != &a)
			{
				Assert(size == a.size);

				Copy(a.data);
			}
			return *this;
		}

		
		template <typename U>
		Vector& operator = (const Vector<U>& a) throw()
		{
			Assert(size == a.size);

			register const U* __restrict src = a.data;
			register T* __restrict dst = data;
			for (register int k = size; k; --k, ++dst, ++src)
			{
				*dst = (T)*src;
			}
			return *this;
		}


		inline T& operator () (int k) const throw()
		{
			Assert(k >= 1);
			Assert(k <= size);

			return entry[k];
		}


		inline T& Entry(int k) const throw()
		{
			Assert(k >= 1);
			Assert(k <= size);

			return entry[k];
		}


		void Fill(const T& value) throw()
		{
			register T* __restrict dst = data;
			for (register int k = size; k; --k, ++dst)
			{
				*dst = value;
			}
		}


		void FillSeries(const T& start, const T& increment) throw()
		{
			register T* __restrict dst = data;
			register T value = start;
			for (register int k = size; k; --k, ++dst, value += increment)
			{
				*dst = value;
			}
		}


		void Copy(const T* data) throw()
		{
			Assert(data);

			register const T* __restrict src = data;
			register T* __restrict dst = this->data;
			for (register int k = size; k; --k, ++dst, ++src)
			{
				*dst = *src;
			}
		}


		void RemoveEntry(int k) throw()
		{
			Assert(k >= 1);
			Assert(k <= size);

			register T* __restrict src = data + k;
			register T* __restrict dst = data + k - 1;
			for (register int l = size - k; l; --l, ++dst, ++src)
			{
				*dst = *src;
			}
			--*(int*)&this->size;
		}


		void Resize(int size) throw(Memory::Exception)
		{
			Assert(manage);
			Assert(size >= 0);

			T* new_data = new T[size];
			if (!new_data)
			{
				Throw(Memory::exception);
			}

			delete [] data;
			*(T**)&data = new_data;
			*(T**)&entry = data - 1;
			*(int*)&this->size = size;
		}


		void Sort() throw()
		{
			// Combsort11 http://en.wikipedia.org/wiki/Comb_sort
			bool swapped;
			register int gap = size;
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
				for (register int j = gap + 1; j <= size; ++j)
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
