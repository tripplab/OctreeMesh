// Debug.cpp
// Copyright (C) 2012 Miguel Vargas (miguel.vargas@gmail.com)
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

#include <Basic/Debug.h>


#if defined WIN32

	#ifndef WINVER
		#define WINVER         0x0400
		#define _WIN32_WINNT   0x0400
		#define _WIN32_WINDOWS 0x0400
		#define _WIN32_IE      0x0400
	#endif
	#define WIN32_LEAN_AND_MEAN
	#define NOCOMM
	#define NODEFERWINDOWPOS
	#define NOHELP
	#define NOIME
	#define NOMCX
	#define NOPROFILER
	#define NOSERVICE
	#define NOWH
	#include <windows.h>


	void DebugMessage(const char* message) throw()
	{
		OutputDebugStringA(message);
		OutputDebugStringA("\n");
	}


	void DebugMessage(const char* message, const char* source, const char* line, const char* function) throw()
	{
		OutputDebugStringA(message);
		OutputDebugStringA(", source: \"");
		OutputDebugStringA(source);
		OutputDebugStringA("\", line: ");
		OutputDebugStringA(line);
		OutputDebugStringA(", function: ");
		OutputDebugStringA(function);
		OutputDebugStringA("\n");
	}


	void DebugStop() throw()
	{
		DebugBreak();
	}

#else

	#include <signal.h>
	#include <stdio.h>


	void DebugMessage(const char* message) throw()
	{
		fputs(message, stderr);
		fputs("\n", stderr);
	}


	void DebugMessage(const char* message, const char* source, const char* line, const char* function) throw()
	{
		fputs(message, stderr);
		fputs(", source: \"", stderr);
		fputs(source, stderr);
		fputs("\", line: ", stderr);
		fputs(line, stderr);
		fputs(", function: ", stderr);
		fputs(function, stderr);
		fputs("\n", stderr);
	}


	void DebugStop() throw()
	{
		raise(SIGINT);
	}

#endif
