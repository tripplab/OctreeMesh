// MatFile.cpp
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

#include <File/MatFile.h>


static const int endian_test = 1;


template <>
uint32 MatFile::TypeFullMatrix<double>() throw()
{
	return (*(char*)&endian_test == 1) ? 0 : 1000;
}


template <>
uint32 MatFile::TypeFullMatrix<float>() throw()
{
	return (*(char*)&endian_test == 1) ? 10 : 1010;
}


template <>
uint32 MatFile::TypeFullMatrix<int>() throw()
{
	return (*(char*)&endian_test == 1) ? 20 : 1020;
}


template <>
uint32 MatFile::TypeFullMatrix<short>() throw()
{
	return (*(char*)&endian_test == 1) ? 30 : 1030;
}


template <>
uint32 MatFile::TypeFullMatrix<unsigned short>() throw()
{
	return (*(char*)&endian_test == 1) ? 40 : 1040;
}


template <>
uint32 MatFile::TypeFullMatrix<unsigned char>() throw()
{
	return (*(char*)&endian_test == 1) ? 50 : 1050;
}


template <>
uint32 MatFile::TypeFullMatrix<bool>() throw()
{
	return (*(char*)&endian_test == 1) ? 50 : 1050;
}


template <>
uint32 MatFile::TypeSparseMatrix<double>() throw()
{
	return (*(char*)&endian_test == 1) ? 2 : 1002;
}


template <>
uint32 MatFile::TypeSparseMatrix<float>() throw()
{
	return (*(char*)&endian_test == 1) ? 12 : 1012;
}


template <>
uint32 MatFile::TypeSparseMatrix<int>() throw()
{
	return (*(char*)&endian_test == 1) ? 22 : 1022;
}


template <>
uint32 MatFile::TypeSparseMatrix<short>() throw()
{
	return (*(char*)&endian_test == 1) ? 32 : 1032;
}


template <>
uint32 MatFile::TypeSparseMatrix<unsigned short>() throw()
{
	return (*(char*)&endian_test == 1) ? 42 : 1042;
}


template <>
uint32 MatFile::TypeSparseMatrix<unsigned char>() throw()
{
	return (*(char*)&endian_test == 1) ? 52 : 1052;
}


size_t MatFile::BytesPerElement(uint32 type) throw(ExceptionFormat)
{
	switch (type)
	{
		case 0:
		case 2:
		case 1000:
		case 1002:
			return sizeof(double);
		case 10:
		case 12:
		case 1010:
		case 1012:
			return sizeof(float);
		case 20:
		case 22:
		case 1020:
		case 1022:
			return sizeof(int);
		case 30:
		case 32:
		case 1030:
		case 1032:
			return sizeof(short);
		case 40:
		case 42:
		case 1040:
		case 1042:
			return sizeof(unsigned short);
		case 50:
		case 52:
		case 1050:
		case 1052:
			return sizeof(unsigned char);
	}
	Throw(File::exception_format);
}


MatFile::MatFile() throw()
:	File()
{
}


MatFile::~MatFile() throw()
{
}
