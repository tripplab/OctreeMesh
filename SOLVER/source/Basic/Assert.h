// Assert.h
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

#ifndef _Assert_h_
#define _Assert_h_

#ifdef NDEBUG

	#define Assert(condition) ((void)0)

#else

	#include <Basic/Debug.h>
	#include <Basic/Macros.h>

	#define Assert(condition) \
	{ \
		if (!(condition)) \
		{ \
			DebugMessage("Assert fail: " #condition, __FILE__, MacroValueToString(__LINE__), __FUNCTION__); \
			DebugStop(); \
		} \
	}

#endif

#endif
