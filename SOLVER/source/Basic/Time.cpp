// Time.cpp
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

#include <Basic/Time.h>


#if defined(__FreeBSD__) || defined(__APPLE__)
	#include <sys/time.h>

	void GetEpochTime(Time& time) throw()
	{
		struct timeval time_value;
		gettimeofday(&time_value, (struct timezone*)0);
		time.seconds = (long)time_value.tv_sec;
		time.milliseconds = (long)(time_value.tv_usec/1000);
	}
#else
	#include <sys/timeb.h>

	void GetEpochTime(Time& time) throw()
	{
		timeb system_time;
		ftime(&system_time);
		time.seconds = (long)system_time.time;
		time.milliseconds = (long)system_time.millitm;
	}
#endif


void TimeDifference(const Time& a, const Time& b, Time& result) throw()
{
	unsigned long difference = ((unsigned long)(a.seconds - b.seconds)*1000 + (a.milliseconds - b.milliseconds));
	result.seconds = difference/1000;
	result.milliseconds = difference % 1000;
}


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


	void SleepSeconds(unsigned int seconds) throw()
	{
		Sleep((DWORD)seconds*1000);
	}

#else

	#include <unistd.h>


	void SleepSeconds(unsigned int seconds) throw()
	{
		sleep(seconds);
	}

#endif
