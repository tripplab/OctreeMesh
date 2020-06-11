// CSVFile.h
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

#ifndef _CSVFile_h_
#define _CSVFile_h_

#include <Basic/Assert.h>
#include <Basic/Format.h>
#include <Container/CSCMatrix.h>
#include <Container/CSRMatrix.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>
#include <File/File.h>


template <typename T>
void CSVSave(const char* file_name, const CSCMatrix<T>& matrix, const Format& format) throw(File::ExceptionCreate, File::ExceptionWrite)
{
	Assert(file_name);

	try
	{
		File file;
		file.Create(file_name);
		for (register int i = 1; i <= matrix.rows; ++i)
		{
			for (register int j = 1; j <= matrix.columns; ++j)
			{
				file.Write(matrix(i, j), format);
				file.Write((j < matrix.columns) ? ", " : "\n");
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
void CSVSave(const char* file_name, const CSRMatrix<T>& matrix, const Format& format) throw(File::ExceptionCreate, File::ExceptionWrite)
{
	Assert(file_name);

	try
	{
		File file;
		file.Create(file_name);
		for (register int i = 1; i <= matrix.rows; ++i)
		{
			for (register int j = 1; j <= matrix.columns; ++j)
			{
				file.Write(matrix(i, j), format);
				file.Write((j < matrix.columns) ? ", " : "\n");
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
void CSVSave(const char* file_name, const Matrix<T>& matrix, const Format& format) throw(File::ExceptionCreate, File::ExceptionWrite)
{
	Assert(file_name);

	try
	{
		File file;
		file.Create(file_name);
		for (register int i = 1; i <= matrix.rows; ++i)
		{
			for (register int j = 1; j <= matrix.columns; ++j)
			{
				file.Write(matrix.entry[i][j], format);
				file.Write((j < matrix.columns) ? ", " : "\n");
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
void CSVSave(const char* file_name, const Vector<T>& vector, const Format& format, bool column_vector = true) throw(File::ExceptionWrite)
{
	Assert(file_name);

	try
	{
		File file;
		file.Create(file_name);
		if (column_vector)
		{
			for (register int i = 1; i <= vector.size; ++i)
			{
				file.Write(vector.entry[i], format);
				file.Write("\n");
			}
		}
		else
		{
			for (register int j = 1; j <= vector.size; ++j)
			{
				file.Write(vector.entry[j], format);
				file.Write((j < vector.size) ? ", " : "\n");
			}
		}
		file.Close();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}

#endif
