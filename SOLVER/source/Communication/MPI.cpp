// MPI.cpp
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

#include <Basic/Integer.h>
#include <Communication/MPI.h>


MPI::Exception MPI::exception;


template <>
struct MPI::DataType<char>
{
	static const MPI_Datatype id;
	static const int size;
};

const MPI_Datatype MPI::DataType<char>::id = MPI_BYTE;
const int MPI::DataType<char>::size = 1;


template <>
struct MPI::DataType<bool>
{
	static const MPI_Datatype id;
	static const int size;
};

const MPI_Datatype MPI::DataType<bool>::id = MPI_BYTE;
const int MPI::DataType<bool>::size = sizeof(bool);


template <>
struct MPI::DataType<int>
{
	static const MPI_Datatype id;
	static const int size;
};

const MPI_Datatype MPI::DataType<int>::id = MPI_INT;
const int MPI::DataType<int>::size = 1;


template <>
struct MPI::DataType<unsigned int>
{
	static const MPI_Datatype id;
	static const int size;
};

const MPI_Datatype MPI::DataType<unsigned int>::id = MPI_UNSIGNED;
const int MPI::DataType<unsigned int>::size = 1;


template <>
struct MPI::DataType<long>
{
	static const MPI_Datatype id;
	static const int size;
};

const MPI_Datatype MPI::DataType<long>::id = MPI_LONG;
const int MPI::DataType<long>::size = 1;


template <>
struct MPI::DataType<unsigned long>
{
	static const MPI_Datatype id;
	static const int size;
};

const MPI_Datatype MPI::DataType<unsigned long>::id = MPI_UNSIGNED_LONG;
const int MPI::DataType<unsigned long>::size = 1;


template <>
struct MPI::DataType<long long>
{
	static const MPI_Datatype id;
	static const int size;
};

const MPI_Datatype MPI::DataType<long long>::id =  MPI_LONG_LONG;
const int MPI::DataType<long long>::size = 1;


template <>
struct MPI::DataType<unsigned long long>
{
	static const MPI_Datatype id;
	static const int size;
};

const MPI_Datatype MPI::DataType<unsigned long long>::id = MPI_UNSIGNED_LONG_LONG;
const int MPI::DataType<unsigned long long>::size = 1;


template <>
struct MPI::DataType<float>
{
	static const MPI_Datatype id;
	static const int size;
};

const MPI_Datatype MPI::DataType<float>::id = MPI_FLOAT;
const int MPI::DataType<float>::size = 1;


template <>
struct MPI::DataType<double>
{
	static const MPI_Datatype id;
	static const int size;
};

const MPI_Datatype MPI::DataType<double>::id = MPI_DOUBLE;
const int MPI::DataType<double>::size = 1;


MPI::MPI(int argc, char** argv) throw(Exception)
:	size(0),
	rank(0)
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
	{
		Throw(MPI::exception);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}


MPI::~MPI()
{
	MPI_Finalize();
}


int MPI::WaitAnyRequest(const Vector<DataRequest>& request) const throw(Exception)
{
	int p;
	if (MPI_Waitany(request.size, request.data, &p, MPI_STATUS_IGNORE) != MPI_SUCCESS)
	{
		Throw(MPI::exception);
	}
	return p + 1;
}


void MPI::WaitAllRequests(const Vector<DataRequest>& request) const throw(Exception)
{
	if (MPI_Waitall(request.size, request.data, MPI_STATUSES_IGNORE) != MPI_SUCCESS)
	{
		Throw(MPI::exception);
	}
}


bool MPI::Test(DataRequest& request) const throw(Exception)
{
	int flag;
	if (MPI_Test(&request, &flag, MPI_STATUS_IGNORE) != MPI_SUCCESS)
	{
		Throw(MPI::exception);
	}
	return (flag != 0);
}
