// KDTree.h
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

#ifndef _KDTree_h_
#define _KDTree_h_

#include <Basic/Assert.h>
#include <Basic/Memory.h>


template <typename T, typename K, int D>
struct KDNode
{
	KDNode* left;
	KDNode* right;
	T value;
	K key[D];


	K Distance(const K key[D]) throw()
	{
		K distance = 0;
		for (int k = 0; k < D; ++k)
		{
			register K value = this->key[k] - key[k];
			distance += value*value;
		}
		return distance;
	}


	KDNode<T, K, D>* Nearest(const K key[D], K& distance_min, int depth = 0) throw()
	{
		KDNode<T, K, D>* nearest = this;
		distance_min = Distance(key);

		int k = depth % D;
		KDNode<T, K, D>* target_child;
		KDNode<T, K, D>* other_child;
		if (this->key[k] > key[k])
		{
			target_child = left;
			other_child = right;
		}
		else
		{
			target_child = right;
			other_child = left;
		}

		if (target_child)
		{
			K target_distance;
			KDNode<T, K, D>* target_nearest = target_child->Nearest(key, target_distance, depth + 1);
			if (distance_min > target_distance)
			{
				distance_min = target_distance;
				nearest = target_nearest;
			}
		}

		K key_distance = nearest->key[k] - key[k];
		key_distance *= key_distance;
		if (distance_min > key_distance)
		{
			if (other_child)
			{
				K other_distance_min;
				KDNode<T, K, D>* other_nearest = other_child->Nearest(key, other_distance_min, depth + 1);
				if (distance_min > other_distance_min)
				{
					distance_min = other_distance_min;
					nearest = other_nearest;
				}
			}
		}

		return nearest;
	}
};


template <typename T, typename K, int D, int BLOCK_SIZE = 32>
class KDTree
{
	public:

		KDNode<T, K, D>* root;
		int size;


		KDTree() throw()
		:	root((KDNode<T, K, D>*)0),
			size(0),
			last_block((KDNodesBlock*)0),
			index(BLOCK_SIZE - 1)
		{
		}

/*
		KDTree(const KDTree<T, K, D, BLOCK_SIZE>& kdtree) throw(Memory::Exception)
		:	root((KDNode<T, K, D>*)0),
			size(0),
			last_block((KDNodesBlock*)0),
			index(BLOCK_SIZE - 1)
		{
			AppendBranch(kdtree.root, (KDNode<T, K, D>*)0);
		}


		template <int OTHER_BLOCK_SIZE>
		KDTree(const KDTree<T, K, D, OTHER_BLOCK_SIZE>& kdtree) throw(Memory::Exception)
		:	root((KDNode<T, K, D>*)0),
			size(0),
			last_block((KDNodesBlock*)0),
			index(BLOCK_SIZE - 1)
		{
			AppendBranch(kdtree.root, (KDNode<T, K, D>*)0);
		}
*/

		~KDTree() throw()
		{
			while (last_block)
			{
				register KDNodesBlock* __restrict previous = last_block->previous;
				delete last_block;
				last_block = previous;
			}
		}


		KDTree<T, K, D, BLOCK_SIZE>& operator = (const KDTree<T, K, D, BLOCK_SIZE>& kdtree) throw(Memory::Exception)
		{
			if (this != &kdtree)
			{
//				Clear();
//				AppendBranch(kdtree.root, (KDNode<T, K, D>*)0);
			}
			return *this;
		}


		template <int OTHER_BLOCK_SIZE>
		KDTree<T, K, D, BLOCK_SIZE>& operator = (const KDTree<T, K, D, OTHER_BLOCK_SIZE>& kdtree) throw(Memory::Exception)
		{
			if (this != &kdtree)
			{
//				Clear();
//				AppendBranch(kdtree.root, (KDNode<T, K, D>*)0);
			}
			return *this;
		}


		KDNode<T, K, D>* Append(const T& value, const K key[D]) throw(Memory::Exception)
		{
			if (index < BLOCK_SIZE - 1)
			{
				++index;
			}
			else
			{
				register KDNodesBlock* __restrict new_block = new KDNodesBlock;
				if (!new_block)
				{
					Throw(Memory::exception);
				}
				new_block->previous = last_block;
				last_block = new_block;
				index = 0;
			}
			register KDNode<T, K, D>* __restrict new_node = &last_block->data[index];
			new_node->left = (KDNode<T, K, D>*)0;
			new_node->right = (KDNode<T, K, D>*)0;
			new_node->value = value;
			for (register int k = 0; k < D; ++k)
			{
				new_node->key[k] = key[k];
			}

			if (root)
			{
				register KDNode<T, K, D>* __restrict node = root;
				for (register int k = 0; ; k = (k + 1) % D)
				{
					register KDNode<T, K, D>*& side = (node->key[k] > key[k]) ? node->left : node->right;
					if (side)
					{
						node = side;
						continue;
					}
					side = new_node;
					break;
				}
			}
			else
			{
				root = new_node;
			}
			++size;
			return new_node;
		}


		void Clear() throw()
		{
			while (last_block)
			{
				register KDNodesBlock* __restrict previous = last_block->previous;
				delete last_block;
				last_block = previous;
			}
			index = BLOCK_SIZE - 1;
			root = (KDNode<T, K, D>*)0;
			size = 0;
		}


		KDNode<T, K, D>* Nearest(const K key[D], K& distance_min) throw()
		{
			return root->Nearest(key, distance_min);
		}


	private:

		struct KDNodesBlock
		{
			KDNode<T, K, D> data[BLOCK_SIZE];
			KDNodesBlock* previous;
		};

		KDNodesBlock* last_block;
		int index;
};

#endif
