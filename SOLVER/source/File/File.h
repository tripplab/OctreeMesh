// File.h
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

#ifndef _File_h_
#define _File_h_

#include <Basic/Assert.h>
#include <Basic/Debug.h>
#include <Basic/Exception.h>
#include <Basic/Format.h>
#include <Basic/String.h>
#include <stdio.h>

#if defined WIN32
	#if !defined(_OFF64_T_DEFINED) && !defined(__off_t_defined)
		typedef __int64 off_t;
	#endif
#else
	#include <sys/types.h>
#endif


enum NamedPipeType
{
	named_pipe_read_write = 0,
	named_pipe_read_only  = 1,
	named_pipe_write_only = 2
};


class File
{
	public:

		class Exception : public ::Exception {};

		class ExceptionCreate : public ::Exception {};

		class ExceptionEOF : public ::Exception {};

		class ExceptionFlush : public ::Exception {};

		class ExceptionFormat : public ::Exception {};

		class ExceptionOpen : public ::Exception {};

		class ExceptionRead : public ::Exception {};

		class ExceptionSeek : public ::Exception {};

		class ExceptionTell : public ::Exception {};

		class ExceptionWrite : public ::Exception {};


		static Exception exception;

		static ExceptionCreate exception_create;

		static ExceptionEOF exception_eof;

		static ExceptionFlush exception_flush;

		static ExceptionFormat exception_format;

		static ExceptionOpen exception_open;

		static ExceptionRead exception_read;

		static ExceptionSeek exception_seek;

		static ExceptionTell exception_tell;

		static ExceptionWrite exception_write;


		File() throw();

		File(const File& file) throw();

		~File() throw();

		File& operator = (const File& file) throw();

		void Open(const char* file_name) throw(ExceptionOpen);

		void Create(const char* file_name) throw(ExceptionCreate);

		String CreateTemp(const char* file_name_prefix) throw(ExceptionCreate);

		String NamedPipe(const char* pipe_name_postfix, NamedPipeType type) throw(ExceptionCreate);

		void Close() throw();

		void GoBegin(off_t offset_bytes = 0) throw(ExceptionSeek);

		void GoEnd(off_t offset_bytes = 0) throw(ExceptionSeek);

		void GoTo(off_t offset_bytes) throw(ExceptionSeek);

		void Flush() throw(ExceptionEOF);

		template <typename T>
		void Read(T* data, int count) throw(ExceptionEOF, ExceptionRead);

		template <typename T>
		void Write(const T* data, int count) throw(ExceptionWrite);

		template <typename T>
		void Read(T& value) throw(ExceptionEOF, ExceptionRead);

		void Read(bool& value) throw(ExceptionEOF, ExceptionRead);

		template <typename T>
		void Write(const T& data, const Format& format) throw(ExceptionWrite);

		void Read(char& character) throw(ExceptionEOF, ExceptionRead);

		void Read(char* string, int size) throw(ExceptionEOF, ExceptionRead);

		void ReadLine(char* string, int size) throw(ExceptionEOF, ExceptionRead);

		void Write(const char& character) throw(ExceptionWrite);

		void Write(const char* string) throw(ExceptionWrite);

		void SkipComments(char delimiter) throw();

		off_t Tell() throw(ExceptionTell);

	protected:

		FILE* file_stream;

	private:

		template <typename T>
		struct ReadFormat
		{
			static const char* format_string;
		};
};


template <typename T>
void File::Read(T* data, int count) throw(ExceptionEOF, ExceptionRead)
{
	Assert(file_stream);
	Assert(data);

	if (fread(data, sizeof(T)*count, 1, file_stream) != 1)
	{
		if (feof(file_stream))
		{
			Throw(File::exception_eof);
		}
		else
		{
			Throw(File::exception_read);
		}
	}
}


template <typename T>
void File::Write(const T* data, int count) throw(ExceptionWrite)
{
	Assert(file_stream);
	Assert(data);

	if (fwrite(data, sizeof(T)*count, 1, file_stream) != 1)
	{
		Throw(File::exception_write);
	}
}


template <typename T>
void File::Read(T& value) throw(ExceptionEOF, ExceptionRead)
{
	Assert(file_stream);

	if (fscanf(file_stream, ReadFormat<T>::format_string, &value) != 1)
	{
		if (feof(file_stream))
		{
			Throw(File::exception_eof);
		}
		else
		{
			Throw(File::exception_read);
		}
	}
}


template <typename T>
void File::Write(const T& value, const Format& format) throw(ExceptionWrite)
{
	Assert(file_stream);

	if (fprintf(file_stream, format.format_string, value) < 0)
	{
		Throw(File::exception_write);
	}
}

#endif
