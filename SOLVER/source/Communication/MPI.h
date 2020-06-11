// MPI.h
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

#ifndef _MPI_h_
#define _MPI_h_

#include <Basic/Assert.h>
#include <Basic/Debug.h>
#include <Basic/Exception.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>

#define MPICH_SKIP_MPICXX // For MPICH2
#define OMPI_SKIP_MPICXX  // For OpenMPI
#define MPI_NO_CPPBIND    // For SGI MPI
#include <mpi.h>


class MPI
{
	public:

		class Exception : public ::Exception {};

		static Exception exception;

		typedef MPI_Request DataRequest;

		int size;
		int rank;

		MPI(int argc, char** argv) throw(Exception);

		~MPI();

		template <typename T>
		void Send(const int node, const int tag, const T& value) const throw(Exception)
		{
			if (MPI_Send((void*)&value, DataType<T>::size, DataType<T>::id, node, tag, MPI_COMM_WORLD) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void Send(const int node, const int tag, const int count, const T* const data) const throw(Exception)
		{
			Assert(data);
			if (MPI_Send((void*)data, DataType<T>::size*count, DataType<T>::id, node, tag, MPI_COMM_WORLD) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void Send(const int node, const int tag, const Vector<T>& vector) const throw(Exception)
		{
			Assert(vector.data);
			if (MPI_Send((void*)vector.data, DataType<T>::size*vector.size, DataType<T>::id, node, tag, MPI_COMM_WORLD) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void Send(const int node, const int tag, const Matrix<T>& matrix) const throw(Exception)
		{
			Assert(matrix.data);
			if (MPI_Send((void*)matrix.data, DataType<T>::size*matrix.rows*matrix.columns, DataType<T>::id, node, tag, MPI_COMM_WORLD) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void RequestSend(const int node, const int tag, const T& value, DataRequest& request) const throw(Exception)
		{
			if (MPI_Isend((void*)&value, DataType<T>::size, DataType<T>::id, node, tag, MPI_COMM_WORLD, &request) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void RequestSend(const int node, const int tag, const int count, const T* const data, DataRequest& request) const throw(Exception)
		{
			Assert(data);
			if (MPI_Isend((void*)data, DataType<T>::size*count, DataType<T>::id, node, tag, MPI_COMM_WORLD, &request) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void RequestSend(const int node, const int tag, const Vector<T>& vector, DataRequest& request) const throw(Exception)
		{
			if (MPI_Isend((void*)vector.data, DataType<T>::size*vector.size, DataType<T>::id, node, tag, MPI_COMM_WORLD, &request) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void RequestSend(const int node, const int tag, const Matrix<T>& matrix, DataRequest& request) const throw(Exception)
		{
			if (MPI_Isend((void*)matrix.data, DataType<T>::size*matrix.rows*matrix.columns, DataType<T>::id, node, tag, MPI_COMM_WORLD, &request) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void Receive(const int node, const int tag, T& value) const throw(Exception)
		{
			if (MPI_Recv((void*)&value, DataType<T>::size, DataType<T>::id, node, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void Receive(const int node, const int tag, const int count, T* data) const throw(Exception)
		{
			if (MPI_Recv((void*)data, DataType<T>::size*count, DataType<T>::id, node, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void Receive(const int node, const int tag, Vector<T>& vector) const throw(Exception)
		{
			if (MPI_Recv((void*)vector.data, DataType<T>::size*vector.size, DataType<T>::id, node, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void Receive(const int node, const int tag, Matrix<T>& matrix) const throw(Exception)
		{
			if (MPI_Recv((void*)matrix.data, DataType<T>::size*matrix.rows*matrix.columns, DataType<T>::id, node, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void RequestReceive(const int node, const int tag, T& value, DataRequest& request) const throw(Exception)
		{
			if (MPI_Irecv((void*)&value, DataType<T>::size, DataType<T>::id, node, tag, MPI_COMM_WORLD, &request) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void RequestReceive(const int node, const int tag, Vector<T>& vector, DataRequest& request) const throw(Exception)
		{
			if (MPI_Irecv((void*)vector.data, DataType<T>::size*vector.size, DataType<T>::id, node, tag, MPI_COMM_WORLD, &request) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		template <typename T>
		void RequestReceive(const int node, const int tag, Matrix<T>& matrix, DataRequest& request) const throw(Exception)
		{
			if (MPI_Irecv((void*)matrix.data, DataType<T>::size*matrix.rows*matrix.columns, DataType<T>::id, node, tag, MPI_COMM_WORLD, &request) != MPI_SUCCESS)
			{
				Throw(MPI::exception);
			}
		}

		int WaitAnyRequest(const Vector<DataRequest>& request) const throw(Exception);

		void WaitAllRequests(const Vector<DataRequest>& request) const throw(Exception);

		bool Test(DataRequest& request) const throw(Exception);

	protected:

		template <typename T>
		struct DataType
		{
			static const MPI_Datatype id;
			static const int size;
		};
};

#endif
