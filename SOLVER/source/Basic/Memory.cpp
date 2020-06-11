// Memory.cpp
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

#include <Basic/Memory.h>

#if defined(__FreeBSD__) || defined(__APPLE__)
	#include <stdlib.h>
	#include	<malloc/malloc.h>
#elif defined WIN32
	#include <malloc.h>
#else // glibc
	#include <malloc.h>
#endif


size_t Memory::peak_usage = 0;

size_t Memory::current_usage = 0;

#ifdef MEMORY_USAGE
	const bool Memory::memory_usage = true;
#else
	const bool Memory::memory_usage = false;
#endif

Memory::Exception Memory::exception;


void* operator new (register size_t size) throw()
{
	#ifdef MEMORY_USAGE
		void* __restrict memory = malloc(size);
		if (memory)
		{
			#if defined(__FreeBSD__) || defined(__APPLE__)
				size = malloc_size(memory);
			#elif defined WIN32
				size = _msize(memory);
			#else // glibc
				size = malloc_usable_size(memory);
			#endif

			#pragma omp atomic
			Memory::current_usage += size;

			#pragma omp critical
			if (Memory::peak_usage < Memory::current_usage)
			{
				Memory::peak_usage = Memory::current_usage;
			}
		}
		return memory;
	#else
		return malloc(size);
	#endif
}


void* operator new [] (register size_t size) throw()
{
	#ifdef MEMORY_USAGE
		void* __restrict memory = malloc(size);
		if (memory)
		{
			#if defined(__FreeBSD__) || defined(__APPLE__)
				size = malloc_size(memory);
			#elif defined WIN32
				size = _msize(memory);
			#else // glibc
				size = malloc_usable_size(memory);
			#endif

			#pragma omp atomic
			Memory::current_usage += size;

			#pragma omp critical
			if (Memory::peak_usage < Memory::current_usage)
			{
				Memory::peak_usage = Memory::current_usage;
			}
		}
		return memory;
	#else
		return malloc(size);
	#endif
}


void operator delete (register void* __restrict object) throw()
{
	#ifdef MEMORY_USAGE
		if (object)
		{
			#if defined(__FreeBSD__) || defined(__APPLE__)
				size_t size = malloc_size(object);
			#elif defined WIN32
				size_t size = _msize(object);
			#else // glibc
				size_t size = malloc_usable_size(object);
			#endif

			#pragma omp atomic
			Memory::current_usage -= size;
		}
	#endif

	free(object);
}


void operator delete [] (register void* __restrict array) throw()
{
	#ifdef MEMORY_USAGE
		if (array)
		{
			#if defined(__FreeBSD__) || defined(__APPLE__)
				size_t size = malloc_size(array);
			#elif defined WIN32
				size_t size = _msize(array);
			#else // glibc
				size_t size = malloc_usable_size(array);
			#endif

			#pragma omp atomic
			Memory::current_usage -= size;
		}
	#endif

	free(array);
}
