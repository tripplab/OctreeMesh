// Log.cpp
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
#include <Basic/Log.h>
#include <Basic/Time.h>
#include <stdarg.h>


int log_level = 0;

FILE* log_stream = stdout;


void Log(int level, const char* format, ...) throw()
{
	static Time start_time;
	static bool initialized = false;
	if (!initialized)
	{
		GetEpochTime(start_time);
		initialized = true;
	}

	Assert(log_level >= 0);

	if (level <= log_level)
	{
		Time current_time;
		GetEpochTime(current_time);

		Time time_difference;
		TimeDifference(current_time, start_time, time_difference);

		fprintf(log_stream, "[%8lu.%03lu] ", time_difference.seconds, time_difference.milliseconds);

		for (int l = 1; l < level; ++l)
		{
			fputc(' ', log_stream);
		}
		va_list arguments;
		va_start(arguments, format);
		vfprintf(log_stream, format, arguments);
		va_end(arguments);
		fputc('\n', log_stream);
		fflush(log_stream);
	}
}
