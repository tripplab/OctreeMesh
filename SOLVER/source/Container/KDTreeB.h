// KDTreeB.h
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

#ifndef _KDTreeB_h_
#define _KDTreeB_h_

#include <Basic/Assert.h>
#include <Basic/Integer.h>
#include <Basic/Memory.h>
#include <Container/Matrix.h>
#include <Container/Vector.h>


struct FloatFlip
{
	typedef uint64 type;
};


class KDTreeB
{
	public:


		KDTreeB(const Vector<int>& values, const Matrix<float>& coordinates) throw(Memory::Exception)
		:	root((Node*)0),
			nodes(new Node[coordinates.rows])
		{
			Assert(values.size == coordinates.rows);

			if (!nodes)
			{
				Throw(Memory::exception);
			}

			int size = coordinates.rows;
			Node** index = new Node*[size];
			if (!index)
			{
				delete [] nodes;
				Throw(Memory::exception);
			}

			// Convert float to an integer representation
			{
				register const int* __restrict value = values.data;
				register const float* __restrict coordinate = coordinates.data;
				register Node* __restrict node = nodes;
				for (register int i = 0; i < size; ++i, ++node)
				{
					index[i] = node;
					node->value = *(value++);
					for (register int d = 0; d < 2; ++d)
					{
						node->keys[d] = *(coordinate++);
						node->bits[d] ^= -sint32(node->bits[d] >> 31) | 0x80000000; // Michael Herf. Radix Tricks. 2001. http://stereopsis.com/radix.html
					}
				}
			}

			// Build tree
			try
			{
				root = Build(index, size, 0);
			}
			catch (Exception&)
			{
				delete [] index;
				ReThrow();
			}

			// Convert integer representation to float
			{
				register Node* __restrict node = nodes;
				for (register int i = 0; i < size; ++i, ++node)
				{
					for (register int d = 0; d < 2; ++d)
					{
						node->bits[d] ^= (((node->bits[d] >> 31) - 1) | 0x80000000); // Michael Herf. Radix Tricks. 2001. http://stereopsis.com/radix.html
					}
				}
			}

			delete [] index;
		}


		KDTreeB(const Vector<int>& values, const Matrix<double>& coordinates) throw(Memory::Exception)
		:	root((Node*)0),
			nodes(new Node[coordinates.rows])
		{
			Assert(values.size == coordinates.rows);

			if (!nodes)
			{
				Throw(Memory::exception);
			}

			int size = coordinates.rows;
			Node** index = new Node*[size];
			if (!index)
			{
				delete [] nodes;
				Throw(Memory::exception);
			}

			// Convert double to an integer representation
			{
				register const int* __restrict value = values.data;
				register const double* __restrict coordinate = coordinates.data;
				register Node* __restrict node = nodes;
				for (register int i = 0; i < size; ++i, ++node)
				{
					index[i] = node;
					node->value = *(value++);
					for (register int d = 0; d < 2; ++d)
					{
						node->keys[d] = *(coordinate++);
						node->bits[d] ^= -sint64(node->bits[d] >> 63) | 0x8000000000000000; // Michael Herf. Radix Tricks. 2001. http://stereopsis.com/radix.html
					}
				}
			}

			// Build tree
			try
			{
				root = Build(index, size, 0);
			}
			catch (Exception&)
			{
				delete [] index;
				ReThrow();
			}

			// Convert integer representation to double
			{
				register Node* __restrict node = nodes;
				for (register int i = 0; i < size; ++i, ++node)
				{
					for (register int d = 0; d < 2; ++d)
					{
						node->bits[d] ^= (((node->bits[d] >> 63) - 1) | 0x8000000000000000); // Michael Herf. Radix Tricks. 2001. http://stereopsis.com/radix.html
					}
				}
			}

			delete [] index;
		}


		~KDTreeB() throw()
		{
			delete [] nodes;
		}


		int& Nearest(const double coordinate[2], double& distance_min) throw()
		{
			Assert(root);

			distance_min = Float<double>::infinite;
			Node* nearest = (Node*)0;
			root->Nearest(coordinate, distance_min, nearest, 0);
			distance_min = sqrt(distance_min);
			return nearest->value;
		}


	private:


		struct Node
		{
			Node* left;
			Node* right;
			union
			{
				double keys[2];
				FloatFlip::type bits[2];
			};
			int value;


			void Nearest(const double point[2], double& distance_min, Node*& nearest, int depth) throw()
			{
				// A. Duch. Design and Analysis of Multidimensional Data Structures, pp. 29. Departament de Llenguatges i Sistemes Informàtics, Universitat Politècnica de Catalunya. 2004
//printf("%i\n", value);
				// Calculate distance to point
				{
					register double distance_to_point = 0;
					for (register int k = 0; k < 2; ++k)
					{
						register double value = keys[k] - point[k];
						distance_to_point += value*value;
					}
					if (distance_to_point < distance_min)
					{
						distance_min = distance_to_point;
						nearest = this;
					}
				}

				register int k = depth % 2;
				register double point_k = point[k];
				register double keys_k = keys[k];
				register double distance_k = keys_k - point_k;
				distance_k *= distance_k;
				if (point_k < keys_k)
				{
					if (left)
					{
						left->Nearest(point, distance_min, nearest, depth + 1);
					}
					if (distance_k <= distance_min)
					{
						if (right)
						{
							right->Nearest(point, distance_min, nearest, depth + 1);
						}
					}
				}
				else
				{
					if (right)
					{
						right->Nearest(point, distance_min, nearest, depth + 1);
					}
					if (distance_k <= distance_min)
					{
						if (left)
						{
							left->Nearest(point, distance_min, nearest, depth + 1);
						}
					}
				}
			}
		};


		Node* root;
		Node* nodes;


		Node* Build(Node** index, int size, int depth) throw(Memory::Exception)
		{
			int k = depth % 2;

			if (size == 1)
			{
				index[0]->left = (Node*)0;
				index[0]->right = (Node*)0;
				return index[0];
			}
			else if (size == 2)
			{
				if (index[0]->bits[k] < index[1]->bits[k])
				{
					index[0]->left = (Node*)0;
					index[0]->right = (Node*)0;
					index[1]->left = index[0];
					index[1]->right = (Node*)0;
					return index[1];
				}
				else
				{
					index[0]->left = index[1];
					index[0]->right = (Node*)0;
					index[1]->left = (Node*)0;
					index[1]->right = (Node*)0;
					return index[0];
				}
			}
			else if (size == 3)
			{
				register FloatFlip::type a = index[0]->bits[k];
				register FloatFlip::type b = index[1]->bits[k];
				register FloatFlip::type c = index[2]->bits[k];
				register uint8 i;
				register uint8 j;
				register uint8 k;
				if (a < b)
				{
					if (b < c)
					{
						i = 0;
						j = 1;
						k = 2;
					}
					else if (a < c)
					{
						i = 0;
						j = 2;
						k = 1;
					}
					else
					{
						i = 2;
						j = 0;
						k = 1;
					}
				}
				else
				{
					 if (c < b)
					 {
						i = 2;
						j = 1;
						k = 0;
					 }
					 else if (c < a)
					 {
						i = 1;
						j = 2;
						k = 0;
					 }
					 else
					 {
						i = 1;
						j = 0;
						k = 2;
					 }
				}
				index[i]->left = (Node*)0;
				index[i]->right = (Node*)0;
				index[j]->left = index[i];
				index[j]->right = index[k];
				index[k]->left = (Node*)0;
				index[k]->right = (Node*)0;
				return index[j];
			}
			else if (size <= 64)
			{
				CombsortSort(index, size, k);
				int median = size >> 1;
				int left_size = size - median - (size % 2);
				int right_size = size - 1 - left_size;
				index[median]->left = Build(index, left_size, depth + 1);
				index[median]->right = Build(index + median + 1, right_size, depth + 1);
				return index[median];
			}
			else
			{
				RadixSort(index, size, k);
				int median = size >> 1;
				int left_size = size - median - (size % 2);
				int right_size = size - 1 - left_size;
				index[median]->left = Build(index, left_size, depth + 1);
				index[median]->right = Build(index + median + 1, right_size, depth + 1);
				return index[median];
			}
		}


		void CombsortSort(Node** index, int size, int k) throw(Memory::Exception)
		{
			Assert(size <= 64);

			FloatFlip::type array[64];
			for (register int i = 0; i < size; ++i)
			{
				array[i] = index[i]->bits[k];
			}

			register bool swapped;
			register int gap = size;
			do
			{
				gap *= 93;
				gap /= 116;
				if ((gap == 9) || (gap == 10))
				{
					gap = 11;
				}
				else if (gap == 0)
				{
					gap = 1;
				}
				swapped = false;
				for (register int i = gap; i < size; ++i)
				{
					register int j = i - gap;
					if (array[i] < array[j])
					{
						register Node* tmp_index = index[i];
						index[i] = index[j];
						index[j] = tmp_index;

						register FloatFlip::type tmp_array = array[i];
						array[i] = array[j];
						array[j] = tmp_array;
						swapped = true;
					}
				}
			} while (swapped || (gap > 1));
		}


		void RadixSort(Node** index, int size, int k) throw(Memory::Exception)
		{
			uint64* array = new uint64[size];
			uint64* array_b = new uint64[size];
			Node** index_b = new Node*[size];

			if (!array || !array_b || !index_b)
			{
				delete [] index_b;
				delete [] array_b;
				delete [] array;
				Throw(Memory::exception);
			}

			int histogram0[256] = {0};
			int histogram1[256] = {0};
			int histogram2[256] = {0};
			int histogram3[256] = {0};
			int histogram4[256] = {0};
			int histogram5[256] = {0};
			int histogram6[256] = {0};
			int histogram7[256] = {0};

			for (register int i = 0; i < size; ++i)
			{
				register uint64 data = index[i]->bits[k];
				array[i] = data;
				++histogram0[data & 0xff];
				++histogram1[data >> 8 & 0xff];
				++histogram2[data >> 16 & 0xff];
				++histogram3[data >> 24 & 0xff];
				++histogram4[data >> 32 & 0xff];
				++histogram5[data >> 40 & 0xff];
				++histogram6[data >> 48 & 0xff];
				++histogram7[data >> 56];
			}

			{
				register int count0 = 0;
				register int count1 = 0;
				register int count2 = 0;
				register int count3 = 0;
				register int count4 = 0;
				register int count5 = 0;
				register int count6 = 0;
				register int count7 = 0;
				for (register sint16 h = 0; h < 256; ++h)
				{
					register int tmp;

					tmp = histogram0[h];
					histogram0[h] = count0;
					count0 += tmp;

					tmp = histogram1[h];
					histogram1[h] = count1;
					count1 += tmp;

					tmp = histogram2[h];
					histogram2[h] = count2;
					count2 += tmp;

					tmp = histogram3[h];
					histogram3[h] = count3;
					count3 += tmp;

					tmp = histogram4[h];
					histogram4[h] = count4;
					count4 += tmp;

					tmp = histogram5[h];
					histogram5[h] = count5;
					count5 += tmp;

					tmp = histogram6[h];
					histogram6[h] = count6;
					count6 += tmp;

					tmp = histogram7[h];
					histogram7[h] = count7;
					count7 += tmp;
				}
			}

			for (register int i = 0; i < size; ++i)
			{
				register int j = histogram0[array[i] & 0xff]++;
				array_b[j] = array[i];
				index_b[j] = index[i];
			}

			for (register int i = 0; i < size; ++i)
			{
				register int j = histogram1[array_b[i] >> 8 & 0xff]++;
				array[j] = array_b[i];
				index[j] = index_b[i];
			}

			for (register int i = 0; i < size; ++i)
			{
				register int j = histogram2[array[i] >> 16 & 0xff]++;
				array_b[j] = array[i];
				index_b[j] = index[i];
			}

			for (register int i = 0; i < size; ++i)
			{
				register int j = histogram3[array_b[i] >> 24 & 0xff]++;
				array[j] = array_b[i];
				index[j] = index_b[i];
			}

			for (register int i = 0; i < size; ++i)
			{
				register int j = histogram4[array[i] >> 32 & 0xff]++;
				array_b[j] = array[i];
				index_b[j] = index[i];
			}

			for (register int i = 0; i < size; ++i)
			{
				register int j = histogram5[array_b[i] >> 40 & 0xff]++;
				array[j] = array_b[i];
				index[j] = index_b[i];
			}

			for (register int i = 0; i < size; ++i)
			{
				register int j = histogram6[array[i] >> 48 & 0xff]++;
				array_b[j] = array[i];
				index_b[j] = index[i];
			}

			for (register int i = 0; i < size; ++i)
			{
				register int j = histogram7[array_b[i] >> 56 & 0xff]++;
				index[j] = index_b[i];
			}

			delete [] index_b;
			delete [] array_b;
			delete [] array;
		}
};

#endif
