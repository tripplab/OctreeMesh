// Macros.h
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

#ifndef _Macros_h_
#define _Macros_h_

#include <stdlib.h>


// Convert to string
#define MacroToString(a) #a
#define MacroValueToString(a) MacroToString(a)


// Read environment variables
#define GetEnvInteger(env_name) ((getenv(#env_name) != (char*)0) ? atoi(getenv(#env_name)) : env_name)
#define GetEnvFloat(env_name) ((getenv(#env_name) != (char*)0) ? atof(getenv(#env_name)) : env_name)
// Example:
//   #define LOG_LEVEL 2
//   int log_level = GetEnvInteger(LOG_LEVEL);
// Returns the value of the environment variable LOG_LEVEL or the default value defined in LOG_LEVEL

#endif
