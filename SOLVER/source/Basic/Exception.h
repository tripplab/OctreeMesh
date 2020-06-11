// Exception.h
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

#ifndef _Exception_h_
#define _Exception_h_

#include <Basic/Debug.h>
#include <Basic/Macros.h>

#if defined _MSC_VER // Visual C++ does not implement exception specification
	#pragma warning(disable: 4290)
#endif


#define Throw(type) {DebugMessage("Throw: " #type, __FILE__, MacroValueToString(__LINE__), __FUNCTION__); throw type;}

#define ReThrow() {DebugMessage("Re-throw", __FILE__, MacroValueToString(__LINE__), __FUNCTION__); throw;}


class Exception
{
};

#endif
