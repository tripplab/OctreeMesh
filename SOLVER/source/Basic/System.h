// System.h
// Copyright (C) 2013 Miguel Vargas-Felix (miguel.vargas@gmail.com)
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

#ifndef _System_h_
#define _System_h_


#if defined(__APPLE__) && defined(__MACH__)

	#define OS_MacOSX

#elif defined(__CYGWIN__)

	#define OS_Cygwin

#elif defined(__FreeBSD__)

	#define OS_FreeBSD

#elif defined(__linux__)

	#define OS_Linux

#elif defined(_WIN32)

	#define OS_Windows

#else

	#error Unsupported operating system

#endif


#if defined(_MSC_VER)

	#define CC_Microsoft

#elif defined(__clang__)

	#define CC_Clang

#elif defined(__INTEL_COMPILER)

	#define CC_Intel

#elif defined(__GNUC__)

	#define CC_GNU

#else

	#error Unsupported compiler

#endif


#if defined(CC_Microsoft)

	#pragma warning(disable: 4290) // C++ exception specification ignored except to indicate a function is not __declspec(nothrow)
	#pragma warning(disable: 4324) // structure was padded due to __declspec(align())
	#pragma warning(disable: 4456) // declaration hides previous local declaration
	#pragma warning(disable: 4457) // declaration hides function parameter
	#pragma warning(disable: 4458) // declaration hides class member

#elif defined(CC_Clang)

	typedef __SIZE_TYPE__ size_t;

#elif defined(CC_Intel)

	typedef __SIZE_TYPE__ size_t;

	#pragma warning(disable:  981) // Operands are evaluated in unspecified order
	#pragma warning(disable:  383) // Value copied to temporary, reference to temporary used
	#pragma warning(disable: 1418) // External function definition with no prior declaration

#elif defined(CC_GNU)

	#if defined(OS_Windows) // MinGW
		#define __MSVCRT_VERSION__ 0x800
	#endif

#endif

#endif
