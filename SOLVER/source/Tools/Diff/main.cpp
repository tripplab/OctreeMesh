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
#include <Basic/Float.h>
#include <Basic/Integer.h>
#include <Container/Vector.h>
#include <File/File.h>
#include <string.h>
#include <omp.h>


int main(int argc, char** argv)
{
	if (argc != 4)
	{
		fprintf(stderr, "Invalid number of arguments. Use:\n  %s <file_1.post.res> <file_2.post.res> <file_diff.post.res>\n\n", argv[0]);
		return 1;
	}

	try
	{
		FormatInteger format_integer(false, false, 1, FormatInteger::decimal);
		FormatFloat format_float(false, true, 1, 5, FormatFloat::exponential);

		const char* left_file_name = argv[1];
		const char* right_file_name = argv[2];
		const char* diff_file_name = argv[3];

		File left_file;
		File right_file;
		File diff_file;

		left_file.Open(left_file_name);
		right_file.Open(right_file_name);
		diff_file.Create(diff_file_name);

		for ( ; ; )
		{
			char left_buffer[2048];
			char left_char;

			try
			{
				left_file.Read(left_buffer, 2048);
			}
			catch (File::ExceptionEOF&)
			{
				break;
			}
			left_file.Read(left_char);

			char right_buffer[2048];
			char right_char;

			right_file.Read(right_buffer, 2048);
			right_file.Read(right_char);

			if (strcmp(left_buffer, right_buffer) != 0)
			{
				Throw(File::exception_format);
			}
			diff_file.Write(left_buffer);
			diff_file.Write(left_char);

			if (strcmp(left_buffer, "Values") == 0)
			{
				for ( ; ; )
				{
					left_file.Read(left_buffer, 2048);
					left_file.Read(left_char);

					right_file.Read(right_buffer, 2048);
					right_file.Read(right_char);

					if (strcmp(left_buffer, right_buffer) != 0)
					{
						Throw(File::exception_format);
					}

					diff_file.Write(left_buffer);
					diff_file.Write(left_char);

					if (strcmp(left_buffer, "End") == 0)
					{
						left_file.Read(left_buffer, 2048);
						left_file.Read(left_char);

						right_file.Read(right_buffer, 2048);
						right_file.Read(right_char);

						diff_file.Write(left_buffer);
						diff_file.Write(left_char);
						break;
					}

					while (left_char != '\n')
					{
						double left_value;
						double right_value;
						double diff_value;

						left_file.Read(left_value);
						left_file.Read(left_char);

						right_file.Read(right_value);
						right_file.Read(right_char);

						diff_value = left_value - right_value;
						diff_file.Write(diff_value, format_float);
						diff_file.Write(left_char);
					}
				}
			}
		}
		diff_file.Close();
		right_file.Close();
		left_file.Close();
	}
	catch (Exception&)
	{
		DebugPosition("Catch fatal exception");
	}

	return 0;
}
