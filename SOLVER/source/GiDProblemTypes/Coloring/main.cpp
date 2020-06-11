// main.cpp
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

#include <Basic/Debug.h>
#include <Basic/Log.h>
#include <Basic/Macros.h>
#include <Container/Vector.h>
#include <FiniteElement/Geometry.h>

#include <stdlib.h>
#include <omp.h>

#define PREFIX_MAXIMUM_LENGTH 900
#define PATH_MAXIMUM_LENGTH 1000


int main(int argc, char** argv)
{
	if ((argc < 2) || (argc > 4))
	{
		fprintf(stderr, "Invalid number of arguments. Use:\n  %s <problem_files_prefix> [log_level] [log_file]\n", argv[0]);
		fprintf(stderr, " Example:\n  %s gid_examples/problem.gid/problem\n\n", argv[0]);
		fputs(" Required problem files are:\n", stderr);
		fputs("  <problem_file_prefix>.geometry.dat\n", stderr);
		fputs(" By default, log level is set to 2:\n", stderr);
		return 1;
	}

	char geometry_file_name[PATH_MAXIMUM_LENGTH];

	sprintf(geometry_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.geometry.dat", argv[1]);

	log_level = (argc >= 3) ? atoi(argv[2]) : 2;
	if (argc == 4)
	{
		log_stream = fopen(argv[3], "wb");
		if (!log_stream)
		{
			log_stream = stdout;
			fprintf(stderr, "[Warning] Failed to create log file: %s", argv[3]);
		}
	}

	try
	{
		Log(1, "Test ----------------------------------------------------------------");
		Log(1, "-Version:       " MacroValueToString(VERSION));
		Log(1, "-Geometry file: %s", geometry_file_name);

		// Load geometry
		Geometry<double> geometry(geometry_file_name);
		geometry.PrintInfo();

		// Coloring
		Vector<Vector<int> > color_element_index;
		geometry.mesh.ElementColoring(color_element_index);
		Log(1, "Coloring done");
		Log(2, "-Colors used: %i", color_element_index.size);

		// Save mesh partitions
		char mesh_file_name[PATH_MAXIMUM_LENGTH];
		char mesh_name[31];

		for (int c = 1; c <= color_element_index.size; ++c)
		{
			Vector<int>& element_index = color_element_index.entry[c];

			Vector<int> node_index;
			geometry.mesh.GetElementsNodes(element_index, false, node_index);

			sprintf(mesh_file_name, "%." MacroValueToString(PREFIX_MAXIMUM_LENGTH) "s.%03i.post.msh", argv[1], c);
			sprintf(mesh_name, "Color_%03i", c);
			geometry.SaveMesh(mesh_file_name, mesh_name, element_index, node_index, c);
			Log(1, "Mesh saved: %s", mesh_file_name);
		}

		if (Memory::memory_usage)
		{
			Log(1, "Peak allocated memory: %lu bytes", (unsigned long)Memory::peak_usage);
		}
	}
	catch (Exception&)
	{
		DebugPosition("Catch fatal exception");
	}

	if (Memory::current_usage != 0)
	{
		fprintf(stderr, "[Error] Memory leak: %lu bytes", (unsigned long)Memory::current_usage);
	}

	if (log_stream != stdout)
	{
		fclose(log_stream);
	}

	return 0;
}
