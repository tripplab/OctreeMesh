// CSVector.h
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

#ifndef _CSVector_h_
#define _CSVector_h_

#include <Basic/Assert.h>
#include <Basic/Memory.h>
#include <Basic/System.h>


template <typename T>
class CSVector
{
	public:

		const int size;
		const int count;
		T* const entry;
		int* const index;

		CSVector() throw()
		:	size(0),
			count(0),
			entry((T*)0),
			index((int*)0)
		{
			--*(T**)&entry;
			--*(int**)&index;
		}

		CSVector(const CSVector& a) throw(Memory::Exception)
		:	size(a.size),
			count(a.count),
			entry(new T[count]),
			index(new int[count])
		{
			if (!entry || !index)
			{
				delete [] index;
				delete [] entry;
				Throw(Memory::exception);
			}
			--*(T**)&entry;
			--*(int**)&index;
			for (register int k = 1; k <= count; ++k)
			{
				entry[k] = a.entry[k];
				index[k] = a.index[k];
			}
		}

		CSVector(int size) throw(Memory::Exception)
		:	size(size),
			count(0),
			entry((T*)0),
			index((int*)0)
		{
			--*(T**)&entry;
			--*(int**)&index;
		}

		~CSVector() throw()
		{
			++*(T**)&entry;
			++*(int**)&index;
			delete [] index;
			delete [] entry;
		}

		CSVector& operator = (const CSVector& a) throw()
		{
			if (this != &a)
			{
				Copy(a);
			}
			return *this;
		}

		template <typename U>
		CSVector& operator = (const CSVector<U>& a) throw()
		{
			Copy(a);
			return *this;
		}

		inline T& operator () (int k) const throw()
		{
			Assert(k >= 1);
			Assert(k <= size);

			static T zero = 0;
			register int i = Search(k);
			return i ? entry[i] : zero;
		}

		void Allocate(int count) throw(Memory::Exception)
		{
			Assert(count <= size);

			T* new_entry = new T[count];
			int* new_index = new int[count];
			if (!new_entry || !new_index)
			{
				delete [] new_index;
				delete [] new_entry;
				Throw(Memory::exception);
			}

			++*(T**)&entry;
			++*(int**)&index;
			delete [] entry;
			delete [] index;

			*(T**)&entry = new_entry;
			*(int**)&index = new_index;
			--*(T**)&entry;
			--*(int**)&index;
			*(int*)&this->count = count;
		}

		template <typename U>
		void Copy(const CSVector<U>& a) throw()
		{
			Assert(size == a.size);

			if (count != a.count)
			{
				Allocate(a.count);
			}
			for (register int k = 1; k <= count; ++k)
			{
				entry[k] = (T)a.entry[k];
				index[k] = a.index[k];
			}
		}

		inline T& Entry(int q) const throw()
		{
			Assert(q >= 1);
			Assert(q <= count);

			return entry[q];
		}

		void Fill(const T& value) throw()
		{
			for (register int k = 1; k <= count; ++k)
			{
				entry[k] = value;
			}
		}

		void FillSeries(const T& start, const T& increment) throw()
		{
			register T value = start;
			for (register int k = 1; k <= count; ++k, value += increment)
			{
				entry[k] = value;
			}
		}

		inline int& Index(int q) const throw()
		{
			Assert(q >= 1);
			Assert(q <= count);

			return index[q];
		}

		int NonZero() const throw()
		{
			return count;
		}

		void Resize(int size) throw(Memory::Exception)
		{
			Assert(size >= 0);

			++*(T**)&entry;
			++*(int**)&index;
			delete [] index;
			delete [] entry;

			*(T**)&entry = (T*)0;
			*(int**)&index = (int*)0;
			--*(T**)&entry;
			--*(int**)&index;
			*(int*)&count = 0;
			*(int*)&this->size = size;
		}

		int Search(int k) const throw()
		{
			Assert(k >= 1);
			Assert(k <= size);

			// Binary search http://en.wikipedia.org/wiki/Binary_search#Single_comparison_per_iteration
			register int low = 1;
			register int hight = count + 1;
			while (low < hight)
			{
				register int middle = low + ((hight - low) >> 1);
				if (index[middle] < k)
				{
					low = middle + 1;
				}
				else
				{
					hight = middle;
				}
			}
			if ((low <= count) && (index[low] == k))
			{
				return low;
			}
			return 0;
		}

		void SortIndexes() throw()
		{
			// Combsort11 http://en.wikipedia.org/wiki/Comb_sort
			bool swapped;
			register int gap = count;
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
				for (register int j = gap + 1; j <= count; ++j)
				{
					int k = j - gap;
					if (index[j] < index[k])
					{
						register int t = index[j];
						index[j] = index[k];
						index[k] = t;
						swapped = true;
					}
				}
			} while (swapped || (gap > 1));
		}

		void SortByIndexes() throw()
		{
			// Combsort11 http://en.wikipedia.org/wiki/Comb_sort
			bool swapped;
			register int gap = count;
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
				for (register int j = gap + 1; j <= count; ++j)
				{
					int k = j - gap;
					if (index[j] < index[k])
					{
						register int t = index[j];
						index[j] = index[k];
						index[k] = t;
						register T v = entry[j];
						entry[j] = entry[k];
						entry[k] = v;
						swapped = true;
					}
				}
			} while (swapped || (gap > 1));
		}
};


template <typename T>
T SparseDotMultiplication(const int* a_index, const T* a_entry, int a_count, const int* b_index, const T* b_entry, int b_count) throw()
{
	T sum = 0;

	register int qa = 1;
	register int qb = 1;
	register int ka = a_index[qa];
	register int kb = b_index[qb];
	for (bool next = true; next; )
	{
		while (ka < kb)
		{
			++qa;
			if (qa > a_count)
			{
				next = false;
				break;
			}
			ka = a_index[qa];
		}
		if (!next)
		{
			break;
		}
		while (ka > kb)
		{
			++qb;
			if (qb > b_count)
			{
				next = false;
				break;
			}
			kb = b_index[qb];
		}
		if (!next)
		{
			break;
		}
		while (ka == kb)
		{
			sum += a_entry[qa]*b_entry[qb];
			++qa;
			if (qa > a_count)
			{
				next = false;
				break;
			}
			++qb;
			if (qb > b_count)
			{
				next = false;
				break;
			}
			ka = a_index[qa];
			kb = b_index[qb];
		}
	}
	return sum;
}


template <typename T>
T SparseDotMultiplication(const CSVector<T>& X, const CSVector<T>& Y)
{
	return SparseDotMultiplication(X.index, X.entry, X.count, Y.index, Y.entry, Y.count);
}

#endif
