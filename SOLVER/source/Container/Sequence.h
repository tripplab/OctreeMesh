// Sequence.h
// Copyright (C) 2013 Miguel Vargas-Felix (miguel.vargas@gmail.com)
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

#ifndef _Sequence_h_
#define _Sequence_h_

#include <Basic/Assert.h>
#include <Basic/Memory.h>


template <typename TYPE>
struct SequenceItem
{
	SequenceItem* next;

	TYPE value;
};


template <typename TYPE, int BLOCK_SIZE = 8>
class Sequence
{
	public:

		SequenceItem<TYPE>* first;

		SequenceItem<TYPE>* last;


	private:

		struct SequenceBlock
		{
			SequenceItem<TYPE> data[BLOCK_SIZE];
			SequenceBlock* previous;
		};

		SequenceBlock* last_block;

		int index;


	public:

		int size;


		Sequence() throw()
		:	first((SequenceItem<TYPE>*)0),
			last((SequenceItem<TYPE>*)0),
			last_block((SequenceBlock*)0),
			index(BLOCK_SIZE - 1),
			size(0)
		{
		}


		Sequence(const Sequence<TYPE, BLOCK_SIZE>& sequence) throw(Memory::Exception)
		:	first((SequenceItem<TYPE>*)0),
			last((SequenceItem<TYPE>*)0),
			last_block((SequenceBlock*)0),
			index(BLOCK_SIZE - 1),
			size(0)
		{
			for (register SequenceItem<TYPE>* __restrict item = sequence.first; item; item = item->next)
			{
				AppendLast(item->value);
			}
		}


		template <int BLOCK_SIZE2>
		Sequence(const Sequence<TYPE, BLOCK_SIZE2>& sequence) throw(Memory::Exception)
		:	first((SequenceItem<TYPE>*)0),
			last((SequenceItem<TYPE>*)0),
			last_block((SequenceBlock*)0),
			index(BLOCK_SIZE - 1),
			size(0)
		{
			for (register SequenceItem<TYPE>* __restrict item = sequence.first; item; item = item->next)
			{
				AppendLast(item->value);
			}
		}


		~Sequence() throw()
		{
			while (last_block)
			{
				register SequenceBlock* __restrict previous = last_block->previous;
				delete last_block;
				last_block = previous;
			}
		}


		Sequence<TYPE, BLOCK_SIZE>& operator = (const Sequence<TYPE, BLOCK_SIZE>& sequence) throw(Memory::Exception)
		{
			if (this != &sequence)
			{
				Clear();
				for (register SequenceItem<TYPE>* __restrict item = sequence.first; item; item = item->next)
				{
					AppendLast(item->value);
				}
			}
			return *this;
		}


		template <int BLOCK_SIZE2>
		Sequence<TYPE, BLOCK_SIZE>& operator = (const Sequence<TYPE, BLOCK_SIZE2>& sequence) throw(Memory::Exception)
		{
			if (this != &sequence)
			{
				Clear();
				for (register SequenceItem<TYPE>* __restrict item = sequence.first; item; item = item->next)
				{
					AppendLast(item->value);
				}
			}
			return *this;
		}


		TYPE& Append() throw(Memory::Exception)
		{
			register SequenceItem<TYPE>* __restrict new_item;
			if (index < BLOCK_SIZE - 1)
			{
				++index;
				new_item = &last_block->data[index];
			}
			else
			{
				register SequenceBlock* __restrict new_block = new SequenceBlock;
				if (!new_block)
				{
					Throw(Memory::Exception());
				}
				if (last_block)
				{
					new_block->previous = last_block;
				}
				else
				{
					new_block->previous = (SequenceBlock*)0;
				}
				last_block = new_block;
				index = 0;
				new_item = last_block->data;
			}

			if (last)
			{
				last->next = new_item;
			}
			else
			{
				first = new_item;
			}
			new_item->next = (SequenceItem<TYPE>*)0;
			last = new_item;
			++size;
			return new_item->value;
		}


		void Clear() throw()
		{
			while (last_block)
			{
				register SequenceBlock* __restrict previous = last_block->previous;
				delete last_block;
				last_block = previous;
			}
			first = (SequenceItem<TYPE>*)0;
			last = (SequenceItem<TYPE>*)0;
			index = BLOCK_SIZE - 1;
			size = 0;
		}


		SequenceItem<TYPE>* Search(const TYPE& value) throw()
		{
			register SequenceItem<TYPE>* __restrict item = first;
			while (item)
			{
				if (item->value == value)
				{
					return item;
				}
				item = item->next;
			}
			return (SequenceItem<TYPE>*)0;
		}
};

#endif
