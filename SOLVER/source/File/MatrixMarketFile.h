// MatrixMarketFile.h
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

#ifndef _MatrixMarketFile_h_
#define _MatrixMarketFile_h_

#include <Basic/Assert.h>
#include <Basic/Format.h>
#include <Basic/Integer.h>
#include <Container/CSCMatrix.h>
#include <Container/CSRMatrix.h>
#include <File/File.h>


template <typename T>
void MatrixMarketSave(const char* file_name, const CSCMatrix<T>& matrix, const Format& format) throw(File::ExceptionCreate, File::ExceptionWrite)
{
	Assert(file_name);

	try
	{
		FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

		File file;
		file.Create(file_name);
		file.Write(matrix.rows, format_integer);
		file.Write(' ');
		file.Write(matrix.columns, format_integer);
		file.Write(' ');
		file.Write(matrix.NonZero(), format_integer);
		file.Write('\n');
		for (register int j = 1; j <= matrix.columns; ++j)
		{
			int count_j = matrix.Count(j);
			for (int c = 1; c <= count_j; ++c)
			{
				int i = matrix.index[j][c];
				file.Write(i, format_integer);
				file.Write(' ');
				file.Write(j, format_integer);
				file.Write(' ');
				file.Write(matrix.entry[j][c], format);
				file.Write('\n');
			}
		}
		file.Close();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatrixMarketSave(const char* file_name, const CSRMatrix<T>& matrix, const Format& format) throw(File::ExceptionCreate, File::ExceptionWrite)
{
	Assert(file_name);

	try
	{
		FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

		File file;
		file.Create(file_name);
		file.Write(matrix.rows, format_integer);
		file.Write(' ');
		file.Write(matrix.columns, format_integer);
		file.Write(' ');
		file.Write(matrix.NonZero(), format_integer);
		file.Write('\n');
		for (register int i = 1; i <= matrix.rows; ++i)
		{
			int count_i = matrix.Count(i);
			for (int c = 1; c <= count_i; ++c)
			{
				int j = matrix.index[i][c];
				file.Write(i, format_integer);
				file.Write(' ');
				file.Write(j, format_integer);
				file.Write(' ');
				file.Write(matrix.entry[i][c], format);
				file.Write('\n');
			}
		}
		file.Close();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatrixMarketLoad(const char* file_name, const CSRMatrix<T>& matrix) throw(File::ExceptionOpen, File::ExceptionEOF, File::ExceptionRead)
{
	Assert(file_name);

	try
	{
		int rows;
		int columns;
		int nnz;

		// Read header
		File file;
		file.Open(file_name);
		file.SkipComments('%');
		file.Read(rows);
		file.Read(columns);
		file.Read(nnz);
		file.Write(matrix.rows, format_integer);
		file.Write(' ');
		file.Write(matrix.columns, format_integer);
		file.Write(' ');
		file.Write(matrix.NonZero(), format_integer);
		file.Write('\n');

		// Read entries
		Vector<int> index_i(nnz);
		Vector<int> index_j(nnz);
		Vector<T> value(nnz);
		Vector<int> count(rows);
		count.Fill(0);
		for (int n = 1; n <= nnz; ++n)
		{
			file.Read(index_i.entry[n]);
			file.Read(index_j.entry[n]);
			file.Read(value.entry[n]);
			++count[index_i.entry[n]];
		}
		file.Close();

		// Allocate space
		matrix.Resize(rows, columns);
		for (int i = 1; i <= rows; ++i)
		{
			matrix.AllocateRow(i, count.entry[i]);
			count.entry[i] = 0;
		}

		// Fill matrix
		for (int n = 1; n <= nnz; ++n)
		{
			int i = index_i.entry[n];
			int j = index_j.entry[n];
			int c = ++count.entry[i];
			matrix.index[i][c] = j;
			matrix.entry[i][c] = value.entry[n];
		}
		matrix.SortIndexes();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatrixMarketLoad(const char* file_name, const CSCMatrix<T>& matrix) throw(File::ExceptionOpen, File::ExceptionEOF, File::ExceptionRead)
{
	Assert(file_name);

	try
	{
		int rows;
		int columns;
		int nnz;

		// Read header
		File file;
		file.Open(file_name);
		file.SkipComments('%');
		file.Read(rows);
		file.Read(columns);
		file.Read(nnz);
		file.Write(matrix.rows, format_integer);
		file.Write(' ');
		file.Write(matrix.columns, format_integer);
		file.Write(' ');
		file.Write(matrix.NonZero(), format_integer);
		file.Write('\n');

		// Read entries
		Vector<int> index_i(nnz);
		Vector<int> index_j(nnz);
		Vector<T> value(nnz);
		Vector<int> count(columns);
		count.Fill(0);
		for (int n = 1; n <= nnz; ++n)
		{
			file.Read(index_i.entry[n]);
			file.Read(index_j.entry[n]);
			file.Read(value.entry[n]);
			++count[index_j.entry[n]];
		}
		file.Close();

		// Allocate space
		matrix.Resize(rows, columns);
		for (int j = 1; j <= columns; ++j)
		{
			matrix.AllocateRow(j, count.entry[j]);
			count.entry[j] = 0;
		}

		// Fill matrix
		for (int n = 1; n <= nnz; ++n)
		{
			int i = index_i.entry[n];
			int j = index_j.entry[n];
			int c = ++count.entry[j];
			matrix.index[j][c] = i;
			matrix.entry[j][c] = value.entry[n];
		}
		matrix.SortIndexes();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}

#endif
