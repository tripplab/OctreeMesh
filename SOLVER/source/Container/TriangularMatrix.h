// TriangularMatrix.h
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

#ifndef _TriangularMatrix_h_
#define _TriangularMatrix_h_

#include <Basic/Assert.h>
#include <Basic/Memory.h>


template <typename T>
class LowerTriangularMatrix
{
	public:

		const int rows;
		const int columns;
		T* const data;
		T** const entry;


		LowerTriangularMatrix() throw()
		:	rows(0),
			columns(0),
			data((T*)0),
			entry((T**)0)
		{
			--*(T***)&entry;
		}


		LowerTriangularMatrix(const LowerTriangularMatrix& a) throw(Memory::Exception)
		:	rows(a.rows),
			columns(a.columns),
			data(new T[((size_t)rows*(rows + 1)) >> 1]),
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
			for (register int i = 1; i <= rows; row += i++)
			{
				entry[i] = row;
			}
			Fill(a.data);
		}


		LowerTriangularMatrix(int rows, int columns) throw(Memory::Exception)
		:	rows(rows),
			columns(columns),
			data(new T[((size_t)rows*(rows + 1)) >> 1]),
			entry(new T*[rows])
		{
			Assert(rows >= 0);
			Assert(rows == columns);

			if (!data || !entry)
			{
				delete [] entry;
				delete [] data;
				Throw(Memory::exception);
			}

			--*(T***)&entry;
			register T* __restrict row = data - 1;
			for (register int i = 1; i <= rows; row += i++)
			{
				entry[i] = row;
			}
		}


		~LowerTriangularMatrix() throw()
		{
			++*(T***)&entry;
			delete [] entry;
			delete [] data;
		}


		LowerTriangularMatrix& operator = (const LowerTriangularMatrix& a) throw()
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
		LowerTriangularMatrix& operator = (const LowerTriangularMatrix<U>& a) throw()
		{
			Assert(rows == a.rows);
			Assert(columns == a.columns);

			register const U* __restrict src = a.data;
			register T* __restrict dst = data;
			for (register size_t k = ((size_t)rows*(rows + 1)) >> 1; k; --k, ++dst, ++src)
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
			Assert(j <= i);

			return entry[i][j];
		}


		void Fill(const T& value) throw()
		{
			register T* __restrict dst = data;
			for (register size_t i = ((size_t)rows*(rows + 1)) >> 1; i; --i, ++dst)
			{
				*dst = value;
			}
		}


		void Fill(const T* data) throw()
		{
			Assert(data);

			register const T* __restrict src = data;
			register T* __restrict dst = this->data;
			for (register size_t i = ((size_t)rows*(rows + 1)) >> 1; i; --i, ++dst, ++src)
			{
				*dst = *src;
			}
		}


		void Resize(int rows, int columns) throw(Memory::Exception)
		{
			Assert(rows >= 0);
			Assert(rows == columns);

			T* new_data = new T[((size_t)rows*(rows + 1)) >> 1];
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
			for (register int i = 1; i <= rows; row += i++)
			{
				entry[i] = row;
			}
			*(int*)&this->rows = rows;
			*(int*)&this->columns = columns;
		}
};


template <typename T>
class UpperTriangularMatrix
{
	public:

		const int rows;
		const int columns;
		T* const data;
		T** const entry;


		UpperTriangularMatrix() throw()
		:	rows(0),
			columns(0),
			data((T*)0),
			entry((T**)0)
		{
			--*(T***)&entry;
		}


		UpperTriangularMatrix(const UpperTriangularMatrix& a) throw(Memory::Exception)
		:	rows(a.rows),
			columns(a.columns),
			data(new T[((size_t)rows*(rows + 1)) >> 1]),
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
			for (register int i = 1, j = rows - 1; i <= rows; ++i, row += j, --j)
			{
				entry[i] = row;
			}
			Fill(a.data);
		}


		UpperTriangularMatrix(int rows, int columns) throw(Memory::Exception)
		:	rows(rows),
			columns(columns),
			data(new T[((size_t)rows*(rows + 1)) >> 1]),
			entry(new T*[rows])
		{
			Assert(rows >= 0);
			Assert(rows == columns);

			if (!data || !entry)
			{
				delete [] entry;
				delete [] data;
				Throw(Memory::exception);
			}

			--*(T***)&entry;
			register T* __restrict row = data - 1;
			for (register int i = 1, j = rows - 1; i <= rows; ++i, row += j, --j)
			{
				entry[i] = row;
			}
		}


		~UpperTriangularMatrix() throw()
		{
			++*(T***)&entry;
			delete [] entry;
			delete [] data;
		}


		UpperTriangularMatrix& operator = (const UpperTriangularMatrix& a) throw()
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
		UpperTriangularMatrix& operator = (const UpperTriangularMatrix<U>& a) throw()
		{
			Assert(rows == a.rows);
			Assert(columns == a.columns);

			register const U* __restrict src = a.data;
			register T* __restrict dst = data;
			for (register size_t k = ((size_t)rows*(rows + 1)) >> 1; k; --k, ++dst, ++src)
			{
				*dst = *src;
			}
			return *this;
		}


		inline T& operator () (int i, int j) const throw()
		{
			Assert(i >= 1);
			Assert(i <= rows);
			Assert(j >= i);
			Assert(j <= columns);

			return entry[i][j];
		}


		void Fill(const T& value) throw()
		{
			register T* __restrict dst = data;
			for (register size_t i = ((size_t)rows*(rows + 1)) >> 1; i; --i, ++dst)
			{
				*dst = value;
			}
		}


		void Fill(const T* data) throw()
		{
			Assert(data);

			register const T* __restrict src = data;
			register T* __restrict dst = this->data;
			for (register size_t i = ((size_t)rows*(rows + 1)) >> 1; i; --i, ++dst, ++src)
			{
				*dst = *src;
			}
		}


		void Resize(int rows, int columns) throw(Memory::Exception)
		{
			Assert(rows >= 0);
			Assert(rows == columns);

			T* new_data = new T[((size_t)rows*(rows + 1)) >> 1];
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
			for (register int i = 1, j = rows - 1; i <= rows; ++i, row += j, --j)
			{
				entry[i] = row;
			}
			*(int*)&this->rows = rows;
			*(int*)&this->columns = columns;
		}
};

#endif
