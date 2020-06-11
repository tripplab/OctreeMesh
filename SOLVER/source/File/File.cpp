// File.cpp
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

#include <Basic/Assert.h>
#include <Basic/Debug.h>
#include <File/File.h>

#ifdef WIN32

	#define WIN32_LEAN_AND_MEAN
	#define NOSERVICE
	#define NOMCX
	#define NOIME
	#include <windows.h>
	#include <fcntl.h>
	#include <io.h>

	#ifndef fseeko
		#define fseeko _fseeki64
	#endif
	#ifndef ftello
		#define ftello _ftelli64
	#endif

#else

	#include <sys/types.h>
	#include <sys/stat.h>
	#include <limits.h>
	#include <unistd.h>

#endif


File::Exception File::exception;

File::ExceptionCreate File::exception_create;

File::ExceptionEOF File::exception_eof;

File::ExceptionFlush File::exception_flush;

File::ExceptionFormat File::exception_format;

File::ExceptionOpen File::exception_open;

File::ExceptionRead File::exception_read;

File::ExceptionSeek File::exception_seek;

File::ExceptionTell File::exception_tell;

File::ExceptionWrite File::exception_write;


File::File() throw()
:	file_stream((FILE*)0)
{
}


File::File(const File& file) throw()
:	file_stream(file.file_stream)
{
}


File::~File() throw()
{
}


File& File::operator = (const File& file) throw()
{
	if (this != &file)
	{
		file_stream = file.file_stream;
	}
	return *this;
}


void File::Open(const char* file_name) throw(ExceptionOpen)
{
	file_stream = fopen(file_name, "rb");
	if (!file_stream)
	{
		Throw(File::exception_open);
	}
}


void File::Close() throw()
{
	Assert(file_stream);

	fclose(file_stream);
	file_stream = (FILE*)0;
}


void File::Create(const char* file_name) throw(ExceptionCreate)
{
	Assert(file_name);

	file_stream = fopen(file_name, "wb");
	if (!file_stream)
	{
		Throw(File::exception_create);
	}
}


String File::CreateTemp(const char* file_name_prefix) throw(ExceptionCreate)
{
	Assert(file_name_prefix);

	#ifdef WIN32
		char file_name_template[MAX_PATH];
		_snprintf(file_name_template, MAX_PATH - 1, "%s\\%s.tmp.XXXXXX", getenv("TEMP"), file_name_prefix);
		char* file_name = _mktemp(file_name_template);
	#else
		char file_name_template[PATH_MAX];
		snprintf(file_name_template, PATH_MAX - 1, "/tmp/%s.tmp.XXXXXX", file_name_prefix);
		char* file_name = mktemp(file_name_template);
	#endif

	file_stream = fopen(file_name, "wb");
	if (!file_stream)
	{
		Throw(File::exception_create);
	}

	return String(file_name);
}


String File::NamedPipe(const char* pipe_name_postfix, NamedPipeType type) throw(ExceptionCreate)
{
	Assert(pipe_name_postfix);

	#ifdef WIN32

		#define PIPE_BUFFER_SIZE 655360
		#define PIPE_TIMEOUT 3600000 // 3600 seconds

		char pipe_name[MAX_PATH];
		_snprintf(pipe_name, MAX_PATH - 1, "\\\\.\\pipe\\%s", pipe_name_postfix);

		DWORD open_mode = FILE_FLAG_FIRST_PIPE_INSTANCE | ((type == named_pipe_read_only) ? PIPE_ACCESS_INBOUND : (type == named_pipe_write_only) ? PIPE_ACCESS_OUTBOUND : PIPE_ACCESS_DUPLEX);
		DWORD pipe_mode =  PIPE_WAIT | ((type == named_pipe_read_only) ? PIPE_READMODE_BYTE : (type == named_pipe_write_only) ? PIPE_TYPE_BYTE : PIPE_READMODE_BYTE | PIPE_TYPE_BYTE);
		HANDLE pipe_handle = CreateNamedPipe(pipe_name, open_mode, pipe_mode, PIPE_UNLIMITED_INSTANCES, PIPE_BUFFER_SIZE, PIPE_BUFFER_SIZE, PIPE_TIMEOUT, NULL);
		if (pipe_handle == INVALID_HANDLE_VALUE)
		{
			Throw(File::exception_create);
		}
		if (ConnectNamedPipe(pipe_handle, NULL) == 0)
		{
			Throw(File::exception_create);
		}

		int flags = (type == named_pipe_read_only) ? _O_RDONLY : 0;
		int file_descriptor = _open_osfhandle((intptr_t)pipe_handle, flags);
		if (file_descriptor == -1)
		{
			Throw(File::exception_create);
		}

		const char* stream_mode = (type == named_pipe_read_only) ? "rb" : (type == named_pipe_write_only) ? "wb" : "w+b";
		file_stream = _fdopen(file_descriptor, stream_mode);
		if (!file_stream)
		{
			Throw(File::exception_open);
		}

	#else

		char pipe_name[PATH_MAX];
		snprintf(pipe_name, PATH_MAX - 1, "/tmp/%s", pipe_name_postfix);

		unlink(pipe_name);
		if (mkfifo(pipe_name, 0600) == -1)
		{
			Throw(File::exception_create);
		}
		const char* stream_mode = (type == named_pipe_read_only) ? "rb" : (type == named_pipe_write_only) ? "wb" : "w+b";
		file_stream = fopen(pipe_name, stream_mode);
		if (!file_stream)
		{
			Throw(File::exception_open);
		}

	#endif

	return String(pipe_name);
}


void File::GoBegin(off_t offset_bytes) throw(ExceptionSeek)
{
	Assert(file_stream);

	if (fseeko(file_stream, offset_bytes, SEEK_SET))
	{
		Throw(File::exception_read);
	}
}


void File::GoEnd(off_t offset_bytes) throw(ExceptionSeek)
{
	Assert(file_stream);

	if (fseeko(file_stream, offset_bytes, SEEK_END))
	{
		Throw(File::exception_seek);
	}
}


void File::GoTo(off_t offset_bytes) throw(ExceptionSeek)
{
	Assert(file_stream);

	if (fseeko(file_stream, offset_bytes, SEEK_CUR))
	{
		Throw(File::exception_seek);
	}
}


void File::Flush() throw(ExceptionEOF)
{
	Assert(file_stream);

	if (fflush(file_stream) != 0)
	{
		Throw(File::exception_flush);
	}
}


void File::Read(bool& value) throw(ExceptionEOF, ExceptionRead)
{
	Assert(file_stream);

	int integer_value;
	if (fscanf(file_stream, "%i", &integer_value) != 1)
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
	value = integer_value != 0;
}


void File::Read(char& character) throw(ExceptionEOF, ExceptionRead)
{
	Assert(file_stream);

	int integer_char = fgetc(file_stream);
	if (integer_char == EOF)
	{
		Throw(File::exception_eof);
	}
	character = (char)integer_char;
}


void File::Read(char* string, int size) throw(ExceptionEOF, ExceptionRead)
{
	Assert(file_stream);
	Assert(string);

	char format_string[32];
	sprintf(format_string, "%%%is", size);
	if (fscanf(file_stream, format_string, string) != 1)
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


void File::ReadLine(char* string, int size) throw(ExceptionEOF, ExceptionRead)
{
	Assert(file_stream);
	Assert(string);

	char format_string[128];
	sprintf(format_string, "%%%i[^\r\n]", size);
	if (fscanf(file_stream, format_string, string) != 1)
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
	fscanf(file_stream, "%*1[\r]");
	fscanf(file_stream, "%*1[\n]");
}


void File::Write(const char& character) throw(ExceptionWrite)
{
	if (fputc(character, file_stream) < 0)
	{
		Throw(File::exception_write);
	}
}


void File::Write(const char* string) throw(ExceptionWrite)
{
	Assert(string);

	if (fputs(string, file_stream) < 0)
	{
		Throw(File::exception_write);
	}
}


void File::SkipComments(char delimiter) throw()
{
	Assert(file_stream);

	char format_string[6] = "%1[ ]";
	format_string[3] = delimiter;

	char comment_start[2];
	if (fscanf(file_stream, "%*[ \t\r\n]") == EOF)
	{
		return;
	}
	while (fscanf(file_stream, format_string, comment_start) == 1)
	{
		if (fscanf(file_stream, "%*[^\r\n]%*[\r\n]") == EOF)
		{
			return;
		}
	}
}


off_t File::Tell() throw(ExceptionTell)
{
	Assert(file_stream);

	off_t offset = ftello(file_stream);
	if (offset == -1)
	{
		Throw(File::exception_tell);
	}
	return offset;
}

template <>
const char* File::ReadFormat<short int>::format_string = "%hi";


template <>
const char* File::ReadFormat<unsigned short int>::format_string = "%hu";


template <>
const char* File::ReadFormat<int>::format_string = "%i";


template <>
const char* File::ReadFormat<unsigned int>::format_string = "%u";


template <>
const char* File::ReadFormat<long int>::format_string = "%li";


template <>
const char* File::ReadFormat<unsigned long int>::format_string = "%lu";


template <>
const char* File::ReadFormat<float>::format_string = "%f";


template <>
const char* File::ReadFormat<double>::format_string = "%lf";


template <>
const char* File::ReadFormat<long double>::format_string = "%Lf";
