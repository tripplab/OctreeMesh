// PBMFile.h
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

#ifndef _PBMFile_h_
#define _PBMFile_h_

#include <Basic/Assert.h>
#include <Basic/Integer.h>
#include <Basic/Memory.h>
#include <Container/CSCMatrix.h>
#include <Container/CSRMatrix.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>
#include <File/File.h>


template <typename T>
void PBMSave(const char* file_name, const CSCMatrix<T>& matrix) throw(Memory::Exception, File::ExceptionCreate, File::ExceptionWrite)
{
	Assert(file_name);

	try
	{
		FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

		File file;
		file.Create(file_name);
		file.Write("P4\n");
		file.Write(matrix.columns, format_integer);
		file.Write(" ");
		file.Write(matrix.rows, format_integer);
		file.Write("\n");

		Matrix<unsigned char> bits(matrix.rows, (matrix.columns >> 3) + ((matrix.columns % 8) ? 1 : 0));
		bits.Fill((unsigned char)0);
		for (int j = 1; j <= matrix.columns; ++j)
		{
			int byte = ((j - 1) >> 3) + 1;
			unsigned char mask = 0x80 >> ((j - 1) % 8);
			int k_max = matrix.Count(j);
			for (int k = 1; k <= k_max; ++k)
			{
				int i = matrix.index[j][k];
				bits.entry[i][byte] |= mask;
			}
		}
		file.Write(bits.data, bits.rows*bits.columns);
		file.Close();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void PBMSave(const char* file_name, const CSRMatrix<T>& matrix) throw(Memory::Exception, File::ExceptionCreate, File::ExceptionWrite)
{
	Assert(file_name);

	try
	{
		FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

		File file;
		file.Create(file_name);
		file.Write("P4\n");
		file.Write(matrix.columns, format_integer);
		file.Write(" ");
		file.Write(matrix.rows, format_integer);
		file.Write("\n");

		Vector<unsigned char> brow((matrix.columns >> 3) + ((matrix.columns % 8) ? 1 : 0));
		for (int i = 1; i <= matrix.rows; ++i)
		{
			brow.Fill((unsigned char)0);
			int* __restrict mindex = matrix.index[i];
			int k_max = matrix.Count(i);
			for (int k = 1; k <= k_max; ++k)
			{
				register int j = mindex[k] - 1;
				brow.entry[(j >> 3) + 1] |= 0x80 >> (j % 8);
			}
			file.Write(brow.data, brow.size);
		}
		file.Close();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void PBMSave(const char* file_name, const Matrix<T>& matrix) throw(Memory::Exception, File::ExceptionCreate, File::ExceptionWrite)
{
	Assert(file_name);

	try
	{
		FormatInteger format_integer(false, false, 1, FormatInteger::decimal);

		File file;
		file.Create(file_name);
		file.Write("P4\n");
		file.Write(matrix.columns, format_integer);
		file.Write(" ");
		file.Write(matrix.rows, format_integer);
		file.Write("\n");

		Vector<unsigned char> brow((matrix.columns >> 3) + ((matrix.columns % 8) ? 1 : 0));
		for (register int i = 1; i <= matrix.rows; ++i)
		{
			brow.Fill((unsigned char)0);
			for (register int j = 1; j <= matrix.columns; ++j)
			{
				if (matrix.entry[i][j])
				{
					brow.entry[((j - 1) >> 3) + 1] |= 0x80 >> ((j - 1) % 8);
				}
			}
			file.Write(brow.data, brow.size);
		}
		file.Close();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}

#endif
