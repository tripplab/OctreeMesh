// MatFile.h
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

#ifndef _MatFile_h_
#define _MatFile_h_

#include <Basic/Assert.h>
#include <Basic/Integer.h>
#include <Basic/Memory.h>
#include <Container/CSCMatrix.h>
#include <Container/CSRMatrix.h>
#include <Container/Matrix.h>
#include <Container/TriangularMatrix.h>
#include <Container/Vector.h>
#include <File/File.h>
#include <string.h>


class MatFile : public File
{
	public:

		MatFile() throw();
		virtual ~MatFile() throw();

		template <typename T>
		void Retrive(const char* name, CSCMatrix<T>& matrix) throw(Memory::Exception, ExceptionEOF, ExceptionFormat, ExceptionSeek, ExceptionRead);

		template <typename T>
		void Retrive(const char* name, CSRMatrix<T>& matrix) throw(Memory::Exception, ExceptionEOF, ExceptionFormat, ExceptionSeek, ExceptionRead);

		template <typename T>
		void Retrive(const char* name, Matrix<T>& matrix) throw(Memory::Exception, ExceptionEOF, ExceptionFormat, ExceptionSeek, ExceptionRead);

		template <typename T>
		void Retrive(const char* name, Vector<T>& vector) throw(Memory::Exception, ExceptionEOF, ExceptionFormat, ExceptionSeek, ExceptionRead);

		template <typename T>
		void Store(const char* name, const CSCMatrix<T>& matrix) throw(Memory::Exception, ExceptionSeek, ExceptionWrite);

		template <typename T>
		void Store(const char* name, const CSRMatrix<T>& matrix) throw(Memory::Exception, ExceptionSeek, ExceptionWrite);

		template <typename T>
		void Store(const char* name, const Matrix<T>& matrix) throw(Memory::Exception, ExceptionSeek, ExceptionWrite);

		template <typename T>
		void Store(const char* name, const Vector<T>& vector, bool column_vector = true) throw(ExceptionSeek, ExceptionWrite);

		template <typename T>
		void Store(const char* name, const LowerTriangularMatrix<T>& matrix) throw(Memory::Exception, ExceptionSeek, ExceptionWrite);

		template <typename T>
		void Store(const char* name, const UpperTriangularMatrix<T>& matrix) throw(Memory::Exception, ExceptionSeek, ExceptionWrite);

	private:

		template <typename T>
		uint32 TypeFullMatrix() throw();

		template <typename T>
		uint32 TypeSparseMatrix() throw();

		size_t BytesPerElement(uint32 type) throw(ExceptionFormat);

		template <typename T, typename U>
		void ReadByType(CSCMatrix<T>& matrix, int nnz) throw(Memory::Exception, ExceptionEOF, ExceptionRead);

		template <typename T, typename U>
		void ReadByType(CSRMatrix<T>& matrix, int nnz) throw(Memory::Exception, ExceptionEOF, ExceptionRead);

		template <typename T, typename U>
		void ReadByType(Matrix<T>& matrix) throw(Memory::Exception, ExceptionEOF, ExceptionRead);

		template <typename T, typename U>
		void ReadByType(Vector<T>& vector) throw(Memory::Exception, ExceptionEOF, ExceptionRead);
};


template <typename T>
void MatFile::Retrive(const char* name, CSCMatrix<T>& matrix) throw(Memory::Exception, ExceptionEOF, ExceptionFormat, ExceptionSeek, ExceptionRead)
{
	Assert(name);

	try
	{
		GoBegin();
		uint32 header[5];
		for ( ; ; )
		{
			Read(header, 5);
			char* header_name = new char[header[4] + 1];
			if (!header_name)
			{
				Throw(Memory::exception);
			}
			Read(header_name, header[4]);
			bool found = strcmp(name, header_name) == 0;
			delete [] header_name;
			if (found)
			{
				if ((header[0] % 10) != 2) // 0 = Full matrix, 1 = Text matrix, 2 = Sparse matrix
				{
					Throw(File::exception_format);
				}
				if (header[2] != 3) // 3 = real matrix, 4 = complex matrix
				{
					Throw(File::exception_format);
				}

				int nnz = header[1] - 1;
				switch ((header[0]/10) % 10) // Data format (P)
				{
					case 0:
					{
						ReadByType<T, double>(matrix, nnz);
						break;
					}
					case 1:
					{
						ReadByType<T, float>(matrix, nnz);
						break;
					}
					case 2:
					{
						ReadByType<T, int>(matrix, nnz);
						break;
					}
					case 3:
					{
						ReadByType<T, short>(matrix, nnz);
						break;
					}
					case 4:
					{
						ReadByType<T, unsigned short>(matrix, nnz);
						break;
					}
					case 5:
					{
						ReadByType<T, unsigned char>(matrix, nnz);
						break;
					}
					default:
					{
						Throw(File::exception_format);
					}
				}
				break;
			}
			else
			{
				GoTo((long)header[1]*(long)header[2]*(long)BytesPerElement(header[0]));
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatFile::Retrive(const char* name, CSRMatrix<T>& matrix) throw(Memory::Exception, ExceptionEOF, ExceptionFormat, ExceptionSeek, ExceptionRead)
{
	Assert(name);

	try
	{
		GoBegin();
		uint32 header[5];
		for ( ; ; )
		{
			Read(header, 5);
			char* header_name = new char[header[4] + 1];
			if (!header_name)
			{
				Throw(Memory::exception);
			}
			Read(header_name, header[4]);
			bool found = strcmp(name, header_name) == 0;
			delete [] header_name;
			if (found)
			{
				if ((header[0] % 10) != 2) // 0 = Full matrix, 1 = Text matrix, 2 = Sparse matrix
				{
					Throw(File::exception_format);
				}
				if (header[2] != 3) // 3 = real matrix, 4 = complex matrix
				{
					Throw(File::exception_format);
				}

				int nnz = header[1] - 1;
				switch ((header[0]/10) % 10) // Data format (P)
				{
					case 0:
					{
						ReadByType<T, double>(matrix, nnz);
						break;
					}
					case 1:
					{
						ReadByType<T, float>(matrix, nnz);
						break;
					}
					case 2:
					{
						ReadByType<T, int>(matrix, nnz);
						break;
					}
					case 3:
					{
						ReadByType<T, short>(matrix, nnz);
						break;
					}
					case 4:
					{
						ReadByType<T, unsigned short>(matrix, nnz);
						break;
					}
					case 5:
					{
						ReadByType<T, unsigned char>(matrix, nnz);
						break;
					}
					default:
					{
						Throw(File::exception_format);
					}
				}
				break;
			}
			else
			{
				GoTo((long)header[1]*(long)header[2]*(long)BytesPerElement(header[0]));
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatFile::Retrive(const char* name, Matrix<T>& matrix) throw(Memory::Exception, ExceptionEOF, ExceptionFormat, ExceptionSeek, ExceptionRead)
{
	Assert(name);

	try
	{
		GoBegin();
		uint32 header[5];
		for ( ; ; )
		{
			Read(header, 5);
			char* header_name = new char[header[4] + 1];
			if (!header_name)
			{
				Throw(Memory::exception);
			}
			Read(header_name, header[4]);
			bool found = strcmp(name, header_name) == 0;
			delete [] header_name;
			if (found)
			{
				if ((header[0] % 10) != 0) // 0 = Full matrix, 1 = Text matrix, 2 = Sparse matrix
				{
					Throw(File::exception_format);
				}
				matrix.Resize(header[1], header[2]);
				switch ((header[0]/10) % 10) // Data format (P)
				{
					case 0:
					{
						ReadByType<T, double>(matrix);
						break;
					}
					case 1:
					{
						ReadByType<T, float>(matrix);
						break;
					}
					case 2:
					{
						ReadByType<T, int>(matrix);
						break;
					}
					case 3:
					{
						ReadByType<T, short>(matrix);
						break;
					}
					case 4:
					{
						ReadByType<T, unsigned short>(matrix);
						break;
					}
					case 5:
					{
						ReadByType<T, unsigned char>(matrix);
						break;
					}
					default:
					{
						Throw(File::exception_format);
					}
				}
				break;
			}
			else
			{
				GoTo((long)header[1]*(long)header[2]*(long)BytesPerElement(header[0]));
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatFile::Retrive(const char* name, Vector<T>& vector) throw(Memory::Exception, ExceptionEOF, ExceptionFormat, ExceptionSeek, ExceptionRead)
{
	Assert(name);

	try
	{
		GoBegin();
		uint32 header[5];
		for ( ; ; )
		{
			Read(header, 5);
			char* header_name = new char[header[4] + 1];
			if (!header_name)
			{
				Throw(Memory::exception);
			}
			Read(header_name, header[4]);
			bool found = strcmp(name, header_name) == 0;
			delete [] header_name;
			if (found)
			{
				if ((header[0] % 10) != 0) // 0 = Full matrix, 1 = Text matrix, 2 = Sparse matrix
				{
					Throw(File::exception_format);
				}
				if (header[1] == 1)
				{
					vector.Resize(header[2]);
				}
				else if (header[2] == 1)
				{
					vector.Resize(header[1]);
				}
				else
				{
					Throw(File::exception_format);
				}
				switch ((header[0]/10) % 10) // Data format (P)
				{
					case 0:
					{
						ReadByType<T, double>(vector);
						break;
					}
					case 1:
					{
						ReadByType<T, float>(vector);
						break;
					}
					case 2:
					{
						ReadByType<T, int>(vector);
						break;
					}
					case 3:
					{
						ReadByType<T, short>(vector);
						break;
					}
					case 4:
					{
						ReadByType<T, unsigned short>(vector);
						break;
					}
					case 5:
					{
						ReadByType<T, unsigned char>(vector);
						break;
					}
					default:
					{
						Throw(File::exception_format);
					}
				}
				break;
			}
			else
			{
				GoTo((long)header[1]*(long)header[2]*(long)BytesPerElement(header[0]));
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatFile::Store(const char* name, const CSCMatrix<T>& matrix) throw(Memory::Exception, ExceptionSeek, ExceptionWrite)
{
	Assert(name);

	try
	{
		int nnz = 0;
		for (register int i = matrix.columns; i; --i)
		{
			nnz += matrix.Count(i);
		}

		uint32 header[5];
		header[0] = TypeSparseMatrix<T>();
		header[1] = nnz + 1;
		header[2] = 3; // 3 = real matrix, 4 = complex matrix
		header[3] = 0;
		header[4] = (uint32)strlen(name) + 1;
		GoEnd();
		Write(header, 5);
		Write(name, header[4]);

		Vector<T> ridx(nnz);
		Vector<T> cidx(nnz);
		Vector<T> value(nnz);
		register int r = 0;
		for (int j = 1; j <= matrix.columns; ++j)
		{
			int matrix_count_j = matrix.Count(j);
			for (register int k = 1; k <= matrix_count_j; ++k)
			{
				int i = matrix.index[j][k];
				++r;
				ridx.entry[r] = i;
				cidx.entry[r] = j;
				value.entry[r] = matrix.entry[j][k];
			}
		}
		Write(ridx.data, nnz);
		T rows = matrix.rows;
		Write(&rows, 1);			
		Write(cidx.data, nnz);
		T columns = matrix.columns;
		Write(&columns, 1);
		Write(value.data, nnz);
		T tmp = 0;
		Write(&tmp, 1);
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatFile::Store(const char* name, const CSRMatrix<T>& matrix) throw(Memory::Exception, ExceptionSeek, ExceptionWrite)
{
	Assert(name);

	try
	{
		CSCMatrix<T> csr(matrix.rows, matrix.columns);
		Convert(matrix, csr);
		Store(name, csr);
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatFile::Store(const char* name, const Matrix<T>& matrix) throw(Memory::Exception, ExceptionSeek, ExceptionWrite)
{
	Assert(name);

	try
	{
		uint32 header[5];
		header[0] = TypeFullMatrix<T>();
		header[1] = matrix.rows;
		header[2] = matrix.columns;
		header[3] = 0;
		header[4] = (uint32)strlen(name) + 1;
		GoEnd();
		Write(header, 5);
		Write(name, header[4]);
		Matrix<T> real(matrix.columns, matrix.rows);
		for (int i = 1; i <= matrix.rows; ++i)
		{
			for (int j = 1; j <= matrix.columns; ++j)
			{
				real.entry[j][i] = matrix.entry[i][j];
			}
		}
		Write(real.data, real.rows*real.columns);
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatFile::Store(const char* name, const Vector<T>& vector, bool column_vector) throw(ExceptionSeek, ExceptionWrite)
{
	Assert(name);

	try
	{
		uint32 header[5];
		header[0] = TypeFullMatrix<T>();
		if (column_vector)
		{
			header[1] = vector.size;
			header[2] = 1;
		}
		else
		{
			header[1] = 1;
			header[2] = vector.size;
		}
		header[3] = 0;
		header[4] = (uint32)strlen(name) + 1;
		GoEnd();
		Write(header, 5);
		Write(name, header[4]);
		Write(vector.data, vector.size);
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatFile::Store(const char* name, const LowerTriangularMatrix<T>& matrix) throw(Memory::Exception, ExceptionSeek, ExceptionWrite)
{
	Assert(name);

	try
	{
		uint32 header[5];
		header[0] = TypeFullMatrix<T>();
		header[1] = matrix.rows;
		header[2] = matrix.columns;
		header[3] = 0;
		header[4] = (uint32)strlen(name) + 1;
		GoEnd();
		Write(header, 5);
		Write(name, header[4]);
		Matrix<T> real(matrix.columns, matrix.rows);
		for (int i = 1; i <= matrix.rows; ++i)
		{
			for (int j = 1; j <= i; ++j)
			{
				real.entry[j][i] = matrix.entry[i][j];
			}
			for (int j = i + 1; j <= matrix.columns; ++j)
			{
				real.entry[j][i] = 0;
			}
		}
		Write(real.data, real.rows*real.columns);
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T>
void MatFile::Store(const char* name, const UpperTriangularMatrix<T>& matrix) throw(Memory::Exception, ExceptionSeek, ExceptionWrite)
{
	Assert(name);

	try
	{
		uint32 header[5];
		header[0] = TypeFullMatrix<T>();
		header[1] = matrix.rows;
		header[2] = matrix.columns;
		header[3] = 0;
		header[4] = (uint32)strlen(name) + 1;
		GoEnd();
		Write(header, 5);
		Write(name, header[4]);
		Matrix<T> real(matrix.columns, matrix.rows);
		for (int i = 1; i <= matrix.rows; ++i)
		{
			for (int j = 1; j < i; ++j)
			{
				real.entry[j][i] = 0;
			}
			for (int j = i; j <= matrix.columns; ++j)
			{
				real.entry[j][i] = matrix.entry[i][j];
			}
		}
		Write(real.data, real.rows*real.columns);
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T, typename U>
void MatFile::ReadByType(CSCMatrix<T>& matrix, int nnz) throw(Memory::Exception, ExceptionEOF, ExceptionRead)
{
	try
	{
		U rows;
		U columns;
		Vector<U> ridx(nnz);
		Vector<U> cidx(nnz);
		Vector<U> value(nnz);

		Read(ridx.data, nnz);
		Read(&rows, 1);			
		Read(cidx.data, nnz);
		Read(&columns, 1);
		Read(value.data, nnz);
		GoTo(sizeof(U));

		Vector<int> column_k((int)columns);
		column_k.Fill(0);
		for (register int r = 1; r <= nnz; ++r)
		{
			register int j = (int)cidx.entry[r];
			++column_k.entry[j];
		}
		matrix.Resize((int)rows, (int)columns);
		for (register int j = 1; j <= matrix.columns; ++j)
		{
			if (column_k.entry[j] > 0)
			{
				matrix.AllocateColumn(j, column_k.entry[j]);
			}
		}

		column_k.Fill(0);
		for (int r = 1; r <= nnz; ++r)
		{
			register int i = (int)ridx.entry[r];
			register int j = (int)cidx.entry[r];
			register int k = ++column_k.entry[j];
			matrix.index[j][k] = i;
			matrix.entry[j][k] = (T)value.entry[r];
		}
		matrix.SortIndexes();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T, typename U>
void MatFile::ReadByType(CSRMatrix<T>& matrix, int nnz) throw(Memory::Exception, ExceptionEOF, ExceptionRead)
{
	try
	{
		U rows;
		U columns;
		Vector<U> ridx(nnz);
		Vector<U> cidx(nnz);
		Vector<U> value(nnz);

		Read(ridx.data, nnz);
		Read(&rows, 1);			
		Read(cidx.data, nnz);
		Read(&columns, 1);
		Read(value.data, nnz);
		GoTo(sizeof(U));

		Vector<int> row_k((int)rows);
		row_k.Fill(0);
		for (register int r = 1; r <= nnz; ++r)
		{
			register int i = (int)ridx.entry[r];
			++row_k.entry[i];
		}
		matrix.Resize((int)rows, (int)columns);
		for (register int i = 1; i <= matrix.rows; ++i)
		{
			if (row_k.entry[i] > 0)
			{
				matrix.AllocateRow(i, row_k.entry[i]);
			}
		}

		row_k.Fill(0);
		for (int r = 1; r <= nnz; ++r)
		{
			register int i = (int)ridx.entry[r];
			register int j = (int)cidx.entry[r];
			register int k = ++row_k.entry[i];
			matrix.index[i][k] = j;
			matrix.entry[i][k] = (T)value.entry[r];
		}
		matrix.SortIndexes();
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T, typename U>
void MatFile::ReadByType(Matrix<T>& matrix) throw(Memory::Exception, ExceptionEOF, ExceptionRead)
{
	try
	{
		Matrix<U> value(matrix.columns, matrix.rows);
		Read(value.data, value.rows*value.columns);
		for (register int i = 1; i <= matrix.rows; ++i)
		{
			for (register int j = 1; j <= matrix.columns; ++j)
			{
				matrix.entry[i][j] = (T)value.entry[j][i];
			}
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}


template <typename T, typename U>
void MatFile::ReadByType(Vector<T>& vector) throw(Memory::Exception, ExceptionEOF, ExceptionRead)
{
	try
	{
		Vector<U> value(vector.size);
		Read(value.data, value.size);
		for (register int i = vector.size; i; --i)
		{
			vector.entry[i] = (T)value.entry[i];
		}
	}
	catch (Exception&)
	{
		ReThrow();
	}
}

#endif
