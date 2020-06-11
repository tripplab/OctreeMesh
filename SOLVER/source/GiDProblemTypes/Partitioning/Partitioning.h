// Partitioning.h
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

#ifndef _Partitioning_h_
#define _Partitioning_h_

#include <Basic/Debug.h>
#include <Basic/Log.h>
#include <File/File.h>
#include <FiniteElement/Partition.h>

#include <string.h>


#define PARTITIONING_SECTION_NAME_MAX_SIZE 256



template <typename T>
class Partitioning
{
	public:

		int partitions_count;


		Partitioning(const char* file_name) throw(Memory::Exception, File::Exception)
		:	partitions_count()
		{
			try
			{
				char section_name[PARTITIONING_SECTION_NAME_MAX_SIZE];

				// Load problem
				File file;
				file.Open(file_name);
				file.SkipComments(';');

				// Load general parameters
				file.Read(section_name, PARTITIONING_SECTION_NAME_MAX_SIZE);
				if (strcmp(section_name, "{General}") != 0)
				{
					Throw(File::exception_format);
				}
				file.SkipComments(';');
				file.Read(partitions_count);
				file.SkipComments(';');

				file.Close();

				Log(1, "Partition ------------------------------------------------------------");
				Log(1, "-Partitions count: %i", partitions_count);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


	private:

		Partitioning& operator = (const Partitioning&) throw()
		{
			return *this;
		}
};

#endif
