// CSCMatrix.h
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

#ifndef _CSCMatrix_h_
#define _CSCMatrix_h_

#include <Basic/Assert.h>
#include <Basic/Debug.h>
#include <Basic/Memory.h>
#include <Basic/System.h>
#include <Container/Vector.h>


template <typename T>
class CSCMatrix
{
	public:

		const int rows;
		const int columns;
		int* const count;
		T** const entry;
		int** const index;


		CSCMatrix() throw()
		:	rows(0),
			columns(0),
			count((int*)0),
			entry((T**)0),
			index((int**)0)
		{
			--*(int**)&count;
			--*(T***)&entry;
			--*(int***)&index;
		}


		CSCMatrix(const CSCMatrix& matrix) throw(Memory::Exception)
		:	rows(matrix.rows),
			columns(matrix.columns),
			count(new int[columns]),
			entry(new T*[columns]),
			index(new int*[columns])
		{
			if (!count || !entry || !index)
			{
				delete [] index;
				delete [] entry;
				delete [] count;
				Throw(Memory::exception);
			}

			--*(int**)&count;
			--*(T***)&entry;
			--*(int***)&index;
			for (register int j = 1; j <= columns; ++j)
			{
				count[j] = 0;
				entry[j] = (T*)0;
				index[j] = (int*)0;
				--entry[j];
				--index[j];
			}

			try
			{
				for (int j = 1; j <= columns; ++j)
				{
					AllocateColumn(j, matrix.count[j]);
					CopyColumnIndexes(j, matrix);
					CopyColumnValues(j, matrix);
				}
			}
			catch (Memory::Exception&)
			{
				for (register int j = 1; j <= columns; ++j)
				{
					++entry[j];
					++index[j];
					delete [] entry[j];
					delete [] index[j];
				}
				++*(int**)&count;
				++*(T***)&entry;
				++*(int***)&index;
				delete [] index;
				delete [] entry;
				delete [] count;
				ReThrow();
			}
		}


		CSCMatrix(int rows, int columns) throw(Memory::Exception)
		:	rows(rows),
			columns(columns),
			count(new int[columns]),
			entry(new T*[columns]),
			index(new int*[columns])
		{
			Assert(columns >= 0);
			Assert(rows >= 0);

			if (!count || !entry || !index)
			{
				delete [] index;
				delete [] entry;
				delete [] count;
				Throw(Memory::exception);
			}

			--*(int**)&count;
			--*(T***)&entry;
			--*(int***)&index;
			for (register int j = 1; j <= columns; ++j)
			{
				count[j] = 0;
				entry[j] = (T*)0;
				index[j] = (int*)0;
				--entry[j];
				--index[j];
			}
		}


		~CSCMatrix() throw()
		{
			for (register int j = columns; j; --j)
			{
				++entry[j];
				delete [] entry[j];
			}
			for (register int j = columns; j; --j)
			{
				++index[j];
				delete [] index[j];
			}
			++*(int**)&count;
			++*(T***)&entry;
			++*(int***)&index;
			delete [] index;
			delete [] entry;
			delete [] count;
		}


		CSCMatrix& operator = (const CSCMatrix& matrix) throw()
		{
			if (this != &matrix)
			{
				Assert(rows == matrix.rows);
				Assert(columns == matrix.columns);

				for (int j = 1; j <= columns; ++j)
				{
					AllocateColumn(j, matrix.count[j]);
					CopyColumnIndexes(j, matrix);
					CopyColumnValues(j, matrix);
				}
			}
			return *this;
		}


		T& operator () (int i, int j) const throw()
		{
			// WARNING! Column indexes must be ordered before using this function (if not, use SortIndexes())
			Assert(i >= 1);
			Assert(i <= rows);
			Assert(j >= 1);
			Assert(j <= columns);

			static T zero = 0;
			register int k = Search(i, j);
			return k ? entry[j][k] : zero;
		}


		void AllocateColumn(int j, int count) throw(Memory::Exception)
		{
			Assert(j >= 1);
			Assert(j <= columns);

			if (this->count[j] != count)
			{
				T* new_entry = new T[count];
				int* new_index = new int[count];
				if (!new_entry || !new_index)
				{
					delete [] new_index;
					delete [] new_entry;
					Throw(Memory::exception);
				}
				
				++entry[j];
				++index[j];
				delete [] entry[j];
				delete [] index[j];

				entry[j] = new_entry;
				index[j] = new_index;
				--entry[j];
				--index[j];
				this->count[j] = count;
			}
		}


		void CopyColumnIndexes(int j, const CSCMatrix& matrix) throw()
		{
			Assert(j >= 1);
			Assert(j <= columns);
			Assert(count[j] == matrix.count[j]);

			register int* __restrict dst = &index[j][1];
			register int* __restrict src = &matrix.index[j][1];
			for (register int k = count[j]; k; --k, ++dst, ++src)
			{
				*dst = *src;
			}
		}


		void CopyColumnValues(int j, const CSCMatrix& matrix) throw()
		{
			Assert(j >= 1);
			Assert(j <= columns);
			Assert(count[j] == matrix.count[j]);

			register T* __restrict dst = &entry[j][1];
			register T* __restrict src = &matrix.entry[j][1];
			for (register int k = count[j]; k; --k, ++dst, ++src)
			{
				*dst = *src;
			}
		}


		void CopyStructure(const CSCMatrix& matrix) throw(Memory::Exception)
		{
			try
			{
				if (this != &matrix)
				{
					Assert(rows == matrix.rows);
					Assert(columns == matrix.columns);

					for (int j = 1; j <= columns; ++j)
					{
						AllocateColumn(j, matrix.count[j]);
						CopyColumnIndexes(j, matrix);
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		inline int& Count(int j) const throw()
		{
			Assert(j >= 1);
			Assert(j <= columns);

			return count[j];
		}


		inline T& Entry(int j, int k) const throw()
		{
			Assert(j >= 1);
			Assert(j <= columns);
			Assert(k >= 1);
			Assert(k <= count[j]);

			return entry[j][k];
		}


		void Fill(const T& value) throw()
		{
			#pragma omp parallel for default(shared) schedule(guided)
			for(int j = 1; j <= columns; ++j)
			{
				register T* __restrict entry_row = entry[j];
				for (register int k = count[j]; k; --k)
				{
					entry_row[k] = value;
				}
			}
		}


		inline int& Index(int j, int k) const throw()
		{
			Assert(j >= 1);
			Assert(j <= columns);
			Assert(k >= 1);
			Assert(k <= count[j]);

			return index[j][k];
		}


		int NonZero() const throw()
		{
			int nnz = 0;
			for(register int j = columns; j; --j)
			{
				nnz += count[j];
			}
			return nnz;
		}


		void RemoveEntry(int i, int j) throw()
		{
			Assert(i >= 1);
			Assert(i <= rows);
			Assert(j >= 1);
			Assert(j <= columns);

			register int k = Search(i, j);
			if (k)
			{
				{
					register T* __restrict src = &entry[j][k + 1];
					register T* __restrict dst = &entry[j][k];
					for (register int l = count[j] - k; l; --l, ++dst, ++src)
					{
						*dst = *src;
					}
				}
				{
					register int* __restrict src = &index[j][k + 1];
					register int* __restrict dst = &index[j][k];
					for (register int l = count[j] - k; l; --l, ++dst, ++src)
					{
						*dst = *src;
					}
				}
				--*(int*)&count[j];
			}
		}


		void Resize(int rows, int columns) throw(Memory::Exception)
		{
			Assert(rows >= 0);
			Assert(columns >= 0);

			int* new_count = new int[columns];
			T** new_entry = new T*[columns];
			int** new_index = new int*[columns];
			if (!new_count || !new_entry || !new_index)
			{
				delete [] new_index;
				delete [] new_entry;
				delete [] new_count;
				Throw(Memory::exception);
			}

			for (register int j = this->columns; j; --j)
			{
				++entry[j];
				++index[j];
				delete [] index[j];
				delete [] entry[j];
			}
			++*(int**)&count;
			++*(T***)&entry;
			++*(int***)&index;
			delete [] index;
			delete [] entry;
			delete [] count;

			*(int**)&count = new_count;
			*(T***)&entry = new_entry;
			*(int***)&index = new_index;
			--*(int**)&count;
			--*(T***)&entry;
			--*(int***)&index;
			for (register int j = columns; j; --j)
			{
				count[j] = 0;
				entry[j] = (T*)0;
				index[j] = (int*)0;
				--entry[j];
				--index[j];
			}
			*(int*)&this->columns = columns;
			*(int*)&this->rows = rows;
		}


		int Search(int i, int j) const throw()
		{
			Assert(i >= 1);
			Assert(i <= rows);
			Assert(j >= 1);
			Assert(j <= columns);

			// Binary search http://en.wikipedia.org/wiki/Binary_search#Single_comparison_per_iteration
			int* __restrict index_j = index[j];
			int size_j = count[j];
			register int low = 1;
			register int hight = size_j + 1;
			while (low < hight)
			{
				register int middle = low + ((hight - low) >> 1);
				if (index_j[middle] < i)
				{
					low = middle + 1;
				}
				else
				{
					hight = middle;
				}
			}
			if ((low <= size_j) && (index_j[low] == i))
			{
				return low;
			}
			return 0;
		}


		void SortIndexes() throw()
		{
			#pragma omp parallel for default(shared) schedule(guided)
			for(int j = 1; j <= columns; ++j)
			{
				SortIndexes(j);
			}
		}


		void SortIndexes(int j) throw()
		{
			int* __restrict index_j = index[j];
			int size_j = count[j];

			// Combsort11 http://en.wikipedia.org/wiki/Comb_sort
			bool swapped;
			register int gap = size_j;
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
				for (register int i = gap + 1; i <= size_j; ++i)
				{
					int k = i - gap;
					if (index_j[i] < index_j[k])
					{
						register int t = index_j[i];
						index_j[i] = index_j[k];
						index_j[k] = t;
						swapped = true;
					}
				}
			} while (swapped || (gap > 1));
		}


		void SortByIndexes() throw()
		{
			#pragma omp parallel for default(shared) schedule(guided)
			for(int j = 1; j <= columns; ++j)
			{
				SortByIndexes(j);
			}
		}


		void SortByIndexes(int j) throw()
		{
			T* __restrict entry_j = entry[j];
			int* __restrict index_j = index[j];
			int size_j = count[j];

			// Combsort11 http://en.wikipedia.org/wiki/Comb_sort
			bool swapped;
			register int gap = size_j;
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
				for (register int i = gap + 1; i <= size_j; ++i)
				{
					int k = i - gap;
					if (index_j[i] < index_j[k])
					{
						register int t = index_j[i];
						index_j[i] = index_j[k];
						index_j[k] = t;
						register T v = entry_j[i];
						entry_j[i] = entry_j[k];
						entry_j[k] = v;
						swapped = true;
					}
				}
			} while (swapped || (gap > 1));
		}
};


template <typename T>
void Transpose(const CSCMatrix<T>& A, CSCMatrix<T>& At) throw(Memory::Exception)
{
	Assert(A.rows == At.columns);
	Assert(A.columns == At.rows);

	try
	{
		Vector<int> At_count(At.columns);
		At_count.Fill(0);
		for (int j = 1; j <= A.columns; ++j)
		{
			int* __restrict A_index_j = A.index[j];
			int A_count_j = A.count[j];
			for (register int c = 1; c <= A_count_j; ++c)
			{
				register int i = A_index_j[c];
				++At_count.entry[i];
			}
		}
		for (int j = 1; j <= At.columns; ++j)
		{
			At.AllocateColumn(j, At_count.entry[j]);
		}
		for (int i = 1; i <= A.columns; ++i)
		{
			T* __restrict A_entry_i = A.entry[i];
			int* __restrict A_index_i = A.index[i];
			int A_count_i = A.Count(i);
			for (register int c = 1; c <= A_count_i; ++c)
			{
				register int j = A_index_i[c];
				At.entry[j][At_count.entry[j]] = A_entry_i[c];
				At.index[j][At_count.entry[j]--] = i;
			}
		}
		At.SortByIndexes();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}

#endif
