// List.h
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

#ifndef _List_h_
#define _List_h_

#include <Basic/Assert.h>
#include <Basic/Memory.h>


template <typename T>
struct ListItem
{
	ListItem* previous;
	ListItem* next;
	T value;
};


template <typename T, int BLOCK_SIZE = 8>
class List
{
	public:

		ListItem<T>* first;
		ListItem<T>* last;
		int size;


		List() throw()
		:	first((ListItem<T>*)0),
			last((ListItem<T>*)0),
			size(0),
			last_block((ListBlock*)0),
			index(BLOCK_SIZE - 1)
		{
		}


		List(const List<T, BLOCK_SIZE>& list) throw(Memory::Exception)
		:	first((ListItem<T>*)0),
			last((ListItem<T>*)0),
			size(0),
			last_block((ListBlock*)0),
			index(BLOCK_SIZE - 1)
		{
			for (register ListItem<T>* __restrict item = list.first; item; item = item->next)
			{
				AppendLast(item->value);
			}
		}


		template <int OTHER_BLOCK_SIZE>
		List(const List<T, OTHER_BLOCK_SIZE>& list) throw(Memory::Exception)
		:	first((ListItem<T>*)0),
			last((ListItem<T>*)0),
			size(0),
			last_block((ListBlock*)0),
			index(BLOCK_SIZE - 1)
		{
			for (register ListItem<T>* __restrict item = list.first; item; item = item->next)
			{
				AppendLast(item->value);
			}
		}


		~List() throw()
		{
			while (last_block)
			{
				register ListBlock* __restrict previous = last_block->previous;
				delete last_block;
				last_block = previous;
			}
		}


		List<T, BLOCK_SIZE>& operator = (const List<T, BLOCK_SIZE>& list) throw(Memory::Exception)
		{
			if (this != &list)
			{
				Clear();
				for (register ListItem<T>* __restrict item = list.first; item; item = item->next)
				{
					AppendLast(item->value);
				}
			}
			return *this;
		}


		template <int OTHER_BLOCK_SIZE>
		List<T, BLOCK_SIZE>& operator = (const List<T, OTHER_BLOCK_SIZE>& list) throw(Memory::Exception)
		{
			if (this != &list)
			{
				Clear();
				for (register ListItem<T>* __restrict item = list.first; item; item = item->next)
				{
					AppendLast(item->value);
				}
			}
			return *this;
		}


		void AppendLast(const T& value) throw(Memory::Exception)
		{
			if (index < BLOCK_SIZE - 1)
			{
				++index;
			}
			else
			{
				register ListBlock* __restrict new_block = new ListBlock;
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
					new_block->previous = (ListBlock*)0;
				}
				last_block = new_block;
				index = 0;
			}
			register ListItem<T>* __restrict new_item = &last_block->data[index];
			if (last)
			{
				new_item->previous = last;
				last->next = new_item;
			}
			else
			{
				new_item->previous = (ListItem<T>*)0;
				first = new_item;
			}
			new_item->next = (ListItem<T>*)0;
			new_item->value = value;
			last = new_item;
			++size;
		}


		void AppendFirst(const T& value) throw(Memory::Exception)
		{
			if (index < BLOCK_SIZE - 1)
			{
				++index;
			}
			else
			{
				register ListBlock* __restrict new_block = new ListBlock;
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
					new_block->previous = (ListBlock*)0;
				}
				last_block = new_block;
				index = 0;
			}
			register ListItem<T>* __restrict new_item = &last_block->data[index];
			if (first)
			{
				new_item->next = first;
				first->previous = new_item;
			}
			else
			{
				new_item->next = (ListItem<T>*)0;
				last = new_item;
			}
			new_item->previous = (ListItem<T>*)0;
			new_item->value = value;
			first = new_item;
			++size;
		}


		void Clear() throw()
		{
			while (last_block)
			{
				register ListBlock* __restrict previous = last_block->previous;
				delete last_block;
				last_block = previous;
			}
			index = BLOCK_SIZE - 1;
			first = (ListItem<T>*)0;
			last = (ListItem<T>*)0;
			size = 0;
		}


		void Delete(ListItem<T>* item) throw()
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
			register ListItem<T>* __restrict last_item = &last_block->data[index];
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
				register ListBlock* __restrict delete_block = last_block;
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


		ListItem<T>* Search(const T& value) throw()
		{
			register ListItem<T>* __restrict item = first;
			while (item)
			{
				if (item->value == value)
				{
					return item;
				}
				item = item->next;
			}
			return (ListItem<T>*)0;
		}


	private:

		struct ListBlock
		{
			ListItem<T> data[BLOCK_SIZE];
			ListBlock* previous;
		};

		ListBlock* last_block;
		int index;

};

#endif
