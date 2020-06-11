// Matrix.h
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

#ifndef _Matrix_h_
#define _Matrix_h_

#include <Basic/Assert.h>
#include <Basic/Debug.h>
#include <Basic/Memory.h>
#include <Basic/System.h>


template <typename T>
class Matrix
{
	public:

		const int rows;
		const int columns;
		T* const data;
		T** const entry;


		Matrix() throw()
		:	rows(0),
			columns(0),
			data((T*)0),
			entry((T**)0)
		{
			--*(T***)&entry;
		}


		Matrix(const Matrix& a) throw(Memory::Exception)
		:	rows(a.rows),
			columns(a.columns),
			data(new T[(size_t)rows*columns]),
			entry(new T*[rows])
		{
			if (!data || !entry)
			{
				delete [] entry;
				delete [] data;
				Throw(Memory::exception);
			}

			--*(T***)&entry;
			register T* __restrict row = data - 1;
			for (register int i = 1; i <= rows; ++i, row += columns)
			{
				entry[i] = row;
			}
			Fill(a.data);
		}


		Matrix(int rows, int columns) throw(Memory::Exception)
		:	rows(rows),
			columns(columns),
			data(new T[(size_t)rows*columns]),
			entry(new T*[rows])
		{
			Assert(rows >= 0);
			Assert(columns >= 0);

			if (!data || !entry)
			{
				delete [] entry;
				delete [] data;
				Throw(Memory::exception);
			}

			--*(T***)&entry;
			register T* __restrict row = data - 1;
			for (register int i = 1; i <= rows; ++i, row += columns)
			{
				entry[i] = row;
			}
		}


		~Matrix() throw()
		{
			++*(T***)&entry;
			delete [] entry;
			delete [] data;
		}


		Matrix& operator = (const Matrix& a) throw()
		{
			if (this != &a)
			{
				Assert(rows == a.rows);
				Assert(columns == a.columns);

				Fill(a.data);
			}
			return *this;
		}


		template <typename U>
		Matrix& operator = (const Matrix<U>& a) throw()
		{
			Assert(rows == a.rows);
			Assert(columns == a.columns);

			register const U* __restrict src = a.data;
			register T* __restrict dst = data;
			for (register size_t k = (size_t)rows*columns; k; --k, ++dst, ++src)
			{
				*dst = *src;
			}
			return *this;
		}


		inline T& operator () (int i, int j) const throw()
		{
			Assert(i >= 1);
			Assert(i <= rows);
			Assert(j >= 1);
			Assert(j <= columns);

			return entry[i][j];
		}


		inline T& Entry(int i, int j) const throw()
		{
			Assert(i >= 1);
			Assert(i <= rows);
			Assert(j >= 1);
			Assert(j <= columns);

			return entry[i][j];
		}


		void Fill(const T& value) throw()
		{
			register T* __restrict dst = data;
			for (register size_t i = (size_t)rows*columns; i; --i, ++dst)
			{
				*dst = value;
			}
		}


		void Fill(const T* data) throw()
		{
			Assert(data);

			register const T* __restrict src = data;
			register T* __restrict dst = this->data;
			for (register size_t i = (size_t)rows*columns; i; --i, ++dst, ++src)
			{
				*dst = *src;
			}
		}


		void Resize(int rows, int columns) throw(Memory::Exception)
		{
			Assert(rows >= 0);
			Assert(columns >= 0);

			T* new_data = new T[rows*columns];
			T** new_entry = new T*[rows];
			if (!new_data || !new_entry)
			{
				delete [] new_entry;
				delete [] new_data;
				Throw(Memory::exception);
			}

			++*(T***)&entry;
			delete [] entry;
			delete [] data;

			*(T**)&data = new_data;
			*(T***)&entry = new_entry;
			--*(T***)&entry;
			register T* __restrict row = data - 1;
			for (register int i = 1; i <= rows; ++i, row += columns)
			{
				entry[i] = row;
			}
			*(int*)&this->rows = rows;
			*(int*)&this->columns = columns;
		}
};


template <typename T>
void Transpose(const Matrix<T>& A, Matrix<T>& At) throw(Memory::Exception)
{
	Assert(A.rows == At.columns);
	Assert(A.columns == At.rows);

	for (register int i = 1; i <= A.rows; ++i)
	{
		register T* __restrict A_entry_i = A.entry[i];
		for (register int j = 1; j <= A.columns; ++j)
		{
			At.entry[j][i] = A_entry_i[j];
		}
	}
}

#endif
