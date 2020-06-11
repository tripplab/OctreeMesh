// Set.h
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

#ifndef _Set_h_
#define _Set_h_

#include <Basic/Assert.h>
#include <Basic/Memory.h>


template <typename T>
struct SetItem
{
	SetItem* previous;
	SetItem* next;
	T value;
};


template <typename T, int BLOCK_SIZE = 8>
class Set
{
	public:

		SetItem<T>* first;
		SetItem<T>* last;
		int size;


		Set() throw()
		:	first((SetItem<T>*)0),
			last((SetItem<T>*)0),
			size(0),
			last_block((SetBlock*)0),
			index(BLOCK_SIZE - 1)
		{
		}


		Set(const Set<T, BLOCK_SIZE>& set) throw(Memory::Exception)
		:	first((SetItem<T>*)0),
			last((SetItem<T>*)0),
			size(0),
			last_block((SetBlock*)0),
			index(BLOCK_SIZE - 1)
		{
			for (register SetItem<T>* __restrict item = set.first; item; item = item->next)
			{
				Append(item->value);
			}
		}


		template <int OTHER_BLOCK_SIZE>
		Set(const Set<T, OTHER_BLOCK_SIZE>& set) throw(Memory::Exception)
		:	first((SetItem<T>*)0),
			last((SetItem<T>*)0),
			size(0),
			last_block((SetBlock*)0),
			index(BLOCK_SIZE - 1)
		{
			for (register SetItem<T>* __restrict item = set.first; item; item = item->next)
			{
				Append(item->value);
			}
		}


		~Set() throw()
		{
			while (last_block)
			{
				register SetBlock* __restrict previous = last_block->previous;
				delete last_block;
				last_block = previous;
			}
		}


		Set<T, BLOCK_SIZE>& operator = (const Set<T, BLOCK_SIZE>& set) throw(Memory::Exception)
		{
			if (this != &set)
			{
				Clear();
				for (register SetItem<T>* __restrict item = set.first; item; item = item->next)
				{
					Append(item->value);
				}
			}
			return *this;
		}


		template <int OTHER_BLOCK_SIZE>
		Set<T, BLOCK_SIZE>& operator = (const Set<T, OTHER_BLOCK_SIZE>& set) throw(Memory::Exception)
		{
			if (this != &set)
			{
				Clear();
				for (register SetItem<T>* __restrict item = set.first; item; item = item->next)
				{
					Append(item->value);
				}
			}
			return *this;
		}


		void Append(const T& value) throw(Memory::Exception)
		{
			register SetItem<T>* __restrict search_item = first;
			while (search_item)
			{
				if (value == search_item->value)
				{
					return;
				}
				if (value < search_item->value)
				{
					break;
				}
				search_item = search_item->next;
			}
			if (index < BLOCK_SIZE - 1)
			{
				++index;
			}
			else
			{
				register SetBlock* __restrict new_block = new SetBlock;
				if (!new_block)
				{
					Throw(Memory::exception);
				}
				if (last_block)
				{
					new_block->previous = last_block;
				}
				else
				{
					new_block->previous = (SetBlock*)0;
				}
				last_block = new_block;
				index = 0;
			}
			register SetItem<T>* __restrict new_item = &last_block->data[index];
			if (search_item)
			{
				if (search_item->previous)
				{
					new_item->previous = search_item->previous;
					search_item->previous->next = new_item;
				}
				else
				{
					new_item->previous = (SetItem<T>*)0;
					first = new_item;
				}
				new_item->next = search_item;
				search_item->previous = new_item;
			}
			else
			{
				if (last)
				{
					new_item->previous = last;
					last->next = new_item;
				}
				else
				{
					new_item->previous = (SetItem<T>*)0;
					first = new_item;
				}
				new_item->next = (SetItem<T>*)0;
				last = new_item;
			}
			new_item->value = value;
			++size;
		}


		void Clear() throw()
		{
			while (last_block)
			{
				register SetBlock* __restrict previous = last_block->previous;
				delete last_block;
				last_block = previous;
			}
			index = BLOCK_SIZE - 1;
			first = (SetItem<T>*)0;
			last = (SetItem<T>*)0;
			size = 0;
		}


		void Delete(SetItem<T>* item) throw()
		{
			Assert(item);

			if (item->previous)
			{
				item->previous->next = item->next;
			}
			else
			{
				first = item->next;
			}
			if (item->next)
			{
				item->next->previous = item->previous;
			}
			else
			{
				last = item->previous;
			}
			register SetItem<T>* __restrict last_item = &last_block->data[index];
			if (item != last_item)
			{
				item->previous = last_item->previous;
				if (last_item->previous)
				{
					last_item->previous->next = item;
				}
				else
				{
					first = item;
				}
				item->next = last_item->next;
				if (last_item->next)
				{
					last_item->next->previous = item;
				}
				else
				{
					last = item;
				}
				item->value = last_item->value;
			}
			if (index == 0)
			{
				register SetBlock* __restrict delete_block = last_block;
				last_block = last_block->previous;
				delete delete_block;
				index = BLOCK_SIZE - 1;
			}
			else
			{
				--index;
			}
			--size;
		}


		SetItem<T>* Search(const T& value) throw()
		{
			register SetItem<T>* __restrict item = first;
			while (item)
			{
				if (item->value == value)
				{
					return item;
				}
				item = item->next;
			}
			return (SetItem<T>*)0;
		}


	private:

		struct SetBlock
		{
			SetItem<T> data[BLOCK_SIZE];
			SetBlock* previous;
		};

		SetBlock* last_block;
		int index;

};

#endif
