// CSRMatrix.h
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

#ifndef _CSRMatrix_h_
#define _CSRMatrix_h_

#include <Basic/Assert.h>
#include <Basic/Memory.h>
#include <Basic/System.h>
#include <Container/CSCMatrix.h>
#include <Container/Vector.h>


template <typename T>
class CSRMatrix
{
	public:

		const int rows;
		const int columns;
		T** const entry;
		int** const index;


		CSRMatrix() throw()
		:	rows(0),
			columns(0),
			entry((T**)0),
			index((int**)0)
		{
			--*(T***)&entry;
			--*(int***)&index;
		}


		CSRMatrix(const CSRMatrix& matrix) throw(Memory::Exception)
		:	rows(matrix.rows),
			columns(matrix.columns),
			entry(new T*[rows]),
			index(new int*[rows])
		{
			if (!entry || !index)
			{
				delete [] index;
				delete [] entry;
				Throw(Memory::exception);
			}

			--*(T***)&entry;
			--*(int***)&index;
			for (register int i = 1; i <= rows; ++i)
			{
				entry[i] = (T*)0;
				index[i] = (int*)0;
				--entry[i];
			}

			try
			{
				for (int i = 1; i <= rows; ++i)
				{
					AllocateRow(i, matrix.Count(i));
					CopyRowIndexes(i, matrix);
					CopyRowValues(i, matrix);
				}
			}
			catch (Memory::Exception&)
			{
				for (register int i = 1; i <= rows; ++i)
				{
					++entry[i];
					delete [] entry[i];
					delete [] index[i];
				}
				++*(T***)&entry;
				++*(int***)&index;
				delete [] index;
				delete [] entry;
				ReThrow();
			}
		}


		CSRMatrix(int rows, int columns) throw(Memory::Exception)
		:	rows(rows),
			columns(columns),
			entry(new T*[rows]),
			index(new int*[rows])
		{
			Assert(rows >= 0);
			Assert(columns >= 0);

			if (!entry || !index)
			{
				delete [] index;
				delete [] entry;
				Throw(Memory::exception);
			}

			--*(T***)&entry;
			--*(int***)&index;
			for (register int i = 1; i <= rows; ++i)
			{
				entry[i] = (T*)0;
				index[i] = (int*)0;
				--entry[i];
			}
		}


		~CSRMatrix() throw()
		{
			for (register int i = rows; i; --i)
			{
				++entry[i];
				delete [] entry[i];
			}
			for (register int i = rows; i; --i)
			{
				delete [] index[i];
			}
			++*(T***)&entry;
			++*(int***)&index;
			delete [] index;
			delete [] entry;
		}


		CSRMatrix& operator = (const CSRMatrix& matrix) throw()
		{
			if (this != &matrix)
			{
				Assert(rows == matrix.rows);
				Assert(columns == matrix.columns);

				for (int i = 1; i <= rows; ++i)
				{
					AllocateRow(i, matrix.Count(i));
					CopyRowIndexes(i, matrix);
					CopyRowValues(i, matrix);
				}
			}
			return *this;
		}


		T& operator () (int i, int j) const throw()
		{
			// WARNING! Row indexes must be ordered before using this function (if not, use SortIndexes())
			Assert(i >= 1);
			Assert(i <= rows);
			Assert(j >= 1);
			Assert(j <= columns);

			static T zero = 0;
			register int k = Search(i, j);
			return k ? entry[i][k] : zero;
		}


		void AllocateRow(int i, int count) throw(Memory::Exception)
		{
			Assert(i >= 1);
			Assert(i <= rows);

			if (this->Count(i) != count)
			{
				T* new_entry = new T[count];
				int* new_index = new int[count + 1];
				if (!new_entry || !new_index)
				{
					delete [] new_index;
					delete [] new_entry;
					Throw(Memory::exception);
				}

				++entry[i];
				delete [] entry[i];
				delete [] index[i];

				entry[i] = new_entry;
				index[i] = new_index;
				--entry[i];
				index[i][0] = count;
			}
		}


		void CopyRowIndexes(int i, const CSRMatrix& matrix) throw()
		{
			Assert(i >= 1);
			Assert(i <= rows);
			Assert(Count(i) == matrix.Count(i));

			register int* __restrict dst = &index[i][1];
			register int* __restrict src = &matrix.index[i][1];
			register int count_i = Count(i);
			for (register int k = count_i; k; --k, ++dst, ++src)
			{
				*dst = *src;
			}
		}


		void CopyRowValues(int i, const CSRMatrix& matrix) throw()
		{
			Assert(i >= 1);
			Assert(i <= rows);
			Assert(Count(i) == matrix.Count(i));

			register T* __restrict dst = &entry[i][1];
			register T* __restrict src = &matrix.entry[i][1];
			register int count_i = Count(i);
			for (register int k = count_i; k; --k, ++dst, ++src)
			{
				*dst = *src;
			}
		}


		void CopyStructure(const CSRMatrix& matrix) throw(Memory::Exception)
		{
			try
			{
				if (this != &matrix)
				{
					Assert(rows == matrix.rows);
					Assert(columns == matrix.columns);

					for (int i = 1; i <= rows; ++i)
					{
						AllocateRow(i, matrix.Count(i));
						CopyRowIndexes(i, matrix);
					}
				}
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		inline int Count(int i) const throw()
		{
			Assert(i >= 1);
			Assert(i <= rows);

			return index[i] ? index[i][0] : 0;
		}


		void Fill(const T& value) throw()
		{
			#pragma omp parallel for default(shared) schedule(guided)
			for(int i = 1; i <= rows; ++i)
			{
				register T* __restrict entry_row = entry[i];
				register int count_i = Count(i);
				for (register int k = count_i; k; --k)
				{
					entry_row[k] = value;
				}
			}
		}


		int NonZero() const throw()
		{
			int nnz = 0;
			for(register int i = rows; i; --i)
			{
				nnz += Count(i);
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
				register int count_i = index[i][0];
				{
					register T* __restrict src = &entry[i][k + 1];
					register T* __restrict dst = &entry[i][k];
					for (register int l = count_i - k; l; --l, ++dst, ++src)
					{
						*dst = *src;
					}
				}
				{
					register int* __restrict src = &index[i][k + 1];
					register int* __restrict dst = &index[i][k];
					for (register int l = count_i - k; l; --l, ++dst, ++src)
					{
						*dst = *src;
					}
				}
				--index[i][0];
			}
		}


		void Resize(int rows, int columns) throw(Memory::Exception)
		{
			Assert(rows >= 0);
			Assert(columns >= 0);

			T** new_entry = new T*[rows];
			int** new_index = new int*[rows];
			if (!new_entry || !new_index)
			{
				delete [] new_index;
				delete [] new_entry;
				Throw(Memory::exception);
			}

			for (register int i = this->rows; i; --i)
			{
				++entry[i];
				delete [] index[i];
				delete [] entry[i];
			}
			++*(T***)&entry;
			++*(int***)&index;
			delete [] index;
			delete [] entry;

			*(T***)&entry = new_entry;
			*(int***)&index = new_index;
			--*(T***)&entry;
			--*(int***)&index;
			for (register int i = rows; i; --i)
			{
				entry[i] = (T*)0;
				index[i] = (int*)0;
				--entry[i];
			}
			*(int*)&this->rows = rows;
			*(int*)&this->columns = columns;
		}


		int Search(int i, int j) const throw()
		{
			Assert(i >= 1);
			Assert(i <= rows);
			Assert(j >= 1);
			Assert(j <= columns);

			// Binary search http://en.wikipedia.org/wiki/Binary_search#Single_comparison_per_iteration
			int* __restrict index_i = index[i];
			int count_i = Count(i);
			register int low = 1;
			register int hight = count_i + 1;
			while (low < hight)
			{
				register int middle = low + ((hight - low) >> 1);
				if (index_i[middle] < j)
				{
					low = middle + 1;
				}
				else
				{
					hight = middle;
				}
			}
			if ((low <= count_i) && (index_i[low] == j))
			{
				return low;
			}
			return 0;
		}


		void SortIndexes() throw()
		{
			#pragma omp parallel for default(shared) schedule(guided)
			for(int i = 1; i <= rows; ++i)
			{
				SortIndexes(i);
			}
		}


		void SortIndexes(int i) throw()
		{
			int* __restrict index_i = index[i];
			int count_i = Count(i);

			// Combsort11 http://en.wikipedia.org/wiki/Comb_sort
			bool swapped;
			register int gap = count_i;
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
				for (register int j = gap + 1; j <= count_i; ++j)
				{
					int k = j - gap;
					if (index_i[j] < index_i[k])
					{
						register int t = index_i[j];
						index_i[j] = index_i[k];
						index_i[k] = t;
						swapped = true;
					}
				}
			} while (swapped || (gap > 1));
		}


		void SortByIndexes() throw()
		{
			#pragma omp parallel for default(shared) schedule(guided)
			for(int i = 1; i <= rows; ++i)
			{
				SortByIndexes(i);
			}
		}


		void SortByIndexes(int i) throw()
		{
			T* __restrict entry_i = entry[i];
			int* __restrict index_i = index[i];
			int count_i = Count(i);

			// Combsort11 http://en.wikipedia.org/wiki/Comb_sort
			bool swapped;
			register int gap = count_i;
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
				for (register int j = gap + 1; j <= count_i; ++j)
				{
					int k = j - gap;
					if (index_i[j] < index_i[k])
					{
						register int t = index_i[j];
						index_i[j] = index_i[k];
						index_i[k] = t;
						register T v = entry_i[j];
						entry_i[j] = entry_i[k];
						entry_i[k] = v;
						swapped = true;
					}
				}
			} while (swapped || (gap > 1));
		}
};


template <typename T>
void Transpose(const CSRMatrix<T>& A, CSRMatrix<T>& At) throw(Memory::Exception)
{
	Assert(A.rows == At.columns);
	Assert(A.columns == At.rows);

	try
	{
		Vector<int> At_count(At.rows);
		At_count.Fill(0);
		for (int i = 1; i <= A.rows; ++i)
		{
			int* __restrict A_index_i = A.index[i];
			int A_count_i = A.Count(i);
			for (register int c = 1; c <= A_count_i; ++c)
			{
				register int j = A_index_i[c];
				++At_count.entry[j];
			}
		}
		for (int i = 1; i <= At.rows; ++i)
		{
			At.AllocateRow(i, At_count.entry[i]);
		}
		for (int j = 1; j <= A.rows; ++j)
		{
			T* __restrict A_entry_j = A.entry[j];
			int* __restrict A_index_j = A.index[j];
			int A_count_j = A.Count(j);
			for (register int c = 1; c <= A_count_j; ++c)
			{
				register int i = A_index_j[c];
				At.entry[i][At_count.entry[i]] = A_entry_j[c];
				At.index[i][At_count.entry[i]--] = j;
			}
		}
		At.SortByIndexes();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void TransposeStructure(const CSRMatrix<T>& A, CSRMatrix<T>& At) throw(Memory::Exception)
{
	Assert(A.rows == At.columns);
	Assert(A.columns == At.rows);

	try
	{
		Vector<int> At_count(At.rows);
		At_count.Fill(0);
		for (int i = 1; i <= A.rows; ++i)
		{
			int* __restrict A_index_i = A.index[i];
			int A_count_i = A.Count(i);
			for (register int c = 1; c <= A_count_i; ++c)
			{
				register int j = A_index_i[c];
				++At_count.entry[j];
			}
		}
		for (int i = 1; i <= At.rows; ++i)
		{
			At.AllocateRow(i, At_count.entry[i]);
		}
		for (int j = 1; j <= A.rows; ++j)
		{
			T* __restrict A_entry_j = A.entry[j];
			int* __restrict A_index_j = A.index[j];
			int A_count_j = A.Count(j);
			for (register int c = 1; c <= A_count_j; ++c)
			{
				register int i = A_index_j[c];
				At.index[i][At_count.entry[i]--] = j;
			}
		}
		At.SortIndexes();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void Transpose(const CSRMatrix<T>& A, CSCMatrix<T>& At) throw(Memory::Exception)
{
	Assert(A.rows == At.columns);
	Assert(A.columns == At.rows);

	try
	{
		#pragma omp parallel for default(shared) schedule(guided)
		for (int i = 1; i <= A.rows; ++i)
		{
			int count = A.Count(i);
			At.AllocateColumn(i, count);
			T* __restrict A_entry_i = A.entry[i];
			int* __restrict A_index_i = A.index[i];
			T* __restrict At_entry_i = At.entry[i];
			int* __restrict At_index_i = At.index[i];
			for (register int q = 1; q <= count; ++q)
			{
				At_entry_i[q] = A_entry_i[q];
				At_index_i[q] = A_index_i[q];
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void TransposeStructure(const CSRMatrix<T>& A, CSCMatrix<T>& At) throw(Memory::Exception)
{
	Assert(A.rows == At.columns);
	Assert(A.columns == At.rows);

	try
	{
		#pragma omp parallel for default(shared) schedule(guided)
		for (int i = 1; i <= A.rows; ++i)
		{
			int count = A.Count(i);
			At.AllocateColumn(i, count);
			int* __restrict A_index_i = A.index[i];
			int* __restrict At_index_i = At.index[i];
			for (register int q = 1; q <= count; ++q)
			{
				At_index_i[q] = A_index_i[q];
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void Transpose(const CSCMatrix<T>& A, CSRMatrix<T>& At) throw(Memory::Exception)
{
	Assert(A.rows == At.columns);
	Assert(A.columns == At.rows);

	try
	{
		#pragma omp parallel for default(shared) schedule(guided)
		for (int j = 1; j <= A.columns; ++j)
		{
			int count = A.Count(j);
			At.AllocateRow(j, count);
			T* __restrict A_entry_j = A.entry[j];
			int* __restrict A_index_j = A.index[j];
			T* __restrict At_entry_j = At.entry[j];
			int* __restrict At_index_j = At.index[j];
			for (register int q = 1; q <= count; ++q)
			{
				At_entry_j[q] = A_entry_j[q];
				At_index_j[q] = A_index_j[q];
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void TransposeStructure(const CSCMatrix<T>& A, CSRMatrix<T>& At) throw(Memory::Exception)
{
	Assert(A.rows == At.columns);
	Assert(A.columns == At.rows);

	try
	{
		#pragma omp parallel for default(shared) schedule(guided)
		for (int j = 1; j <= A.columns; ++j)
		{
			int count = A.Count(j);
			At.AllocateRow(j, count);
			int* __restrict A_index_j = A.index[j];
			int* __restrict At_index_j = At.index[j];
			for (register int q = 1; q <= count; ++q)
			{
				At_index_j[q] = A_index_j[q];
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void Convert(const CSCMatrix<T>& C, CSRMatrix<T>& R) throw(Memory::Exception)
{
	Assert(C.rows == R.rows);
	Assert(C.columns == R.columns);

	try
	{
		Vector<int> R_count(R.rows);
		R_count.Fill(0);
		for (int j = 1; j <= C.columns; ++j)
		{
			int* __restrict C_index_j = C.index[j];
			int C_count_j = C.Count(j);
			for (register int c = 1; c <= C_count_j; ++c)
			{
				register int i = C_index_j[c];
				++R_count.entry[i];
			}
		}
		for (int i = 1; i <= R.rows; ++i)
		{
			R.AllocateRow(i, R_count.entry[i]);
		}
		for (int j = 1; j <= C.columns; ++j)
		{
			T* __restrict C_entry_j = C.entry[j];
			int* __restrict C_index_j = C.index[j];
			int C_count_j = C.Count(j);
			for (register int c = 1; c <= C_count_j; ++c)
			{
				register int i = C_index_j[c];
				R.entry[i][R_count.entry[i]] = C_entry_j[c];
				R.index[i][R_count.entry[i]--] = j;
			}
		}
		R.SortByIndexes();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void Convert(const CSRMatrix<T>& R, CSCMatrix<T>& C) throw(Memory::Exception)
{
	Assert(R.rows == C.rows);
	Assert(R.columns == C.columns);

	try
	{
		Vector<int> C_count(C.columns);
		C_count.Fill(0);
		for (int i = 1; i <= R.rows; ++i)
		{
			int* __restrict R_index_i = R.index[i];
			int R_count_i = R.Count(i);
			for (register int c = 1; c <= R_count_i; ++c)
			{
				register int j = R_index_i[c];
				++C_count.entry[j];
			}
		}
		for (int j = 1; j <= C.columns; ++j)
		{
			C.AllocateColumn(j, C_count.entry[j]);
		}
		for (int i = 1; i <= R.rows; ++i)
		{
			T* __restrict R_entry_i = R.entry[i];
			int* __restrict R_index_i = R.index[i];
			int R_count_i = R.Count(i);
			for (register int c = 1; c <= R_count_i; ++c)
			{
				register int j = R_index_i[c];
				C.entry[j][C_count.entry[j]] = R_entry_i[c];
				C.index[j][C_count.entry[j]--] = i;
			}
		}
		C.SortByIndexes();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}

#endif
