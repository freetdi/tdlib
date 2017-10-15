/*
 * Copyright (C) 2017 Felix Salfelder
 * Authors: Felix Salfelder
 *
 * This file is part of "freetdi", the free tree decomposition intiative
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 *
 * custom allocators
 */
/*--------------------------------------------------------------------------*/
#pragma once
/*--------------------------------------------------------------------------*/
#include <assert.h>
/*--------------------------------------------------------------------------*/
template<unsigned S>
class TRIE_SHARED_AREA{ //
public:
	TRIE_SHARED_AREA() : _seek(NULL), _size(0)
	{
		// incomplete(); hmmm.
	}
	TRIE_SHARED_AREA(const TRIE_SHARED_AREA& p)
		: _seek(p._seek), _size(p._size)
	{ untested();
	}
public: //op
	// 	TRIE_SHARED_AREA&
	// 		: _seek(p._seek), _size(p._size)
	// 	{ untested();
	// 	}
public:
	void reserve(unsigned n){
		if(_seek){ untested();
			return;
		}else{
			assert(!_seek);
			assert(!_size);
			_size = 0;
			_seek = malloc(n * S);
			if(!_seek){ untested();
				throw std::bad_alloc();
			}else{
			}
			_end = (void*) (intptr_t(_seek) + S*n);
		}
	}
	void deallocate() const{
		_seek = (void*) (intptr_t(_seek) - S*_size);
		_size = 0;
	}
	void* allocate(size_t s) const{
		assert(s==S);
		if(_seek!=_end){
			_seek = (void*) (intptr_t(_seek) + S);
			++_size;
			return (void*) (intptr_t(_seek) - S);
		}else{ untested();
			std::cerr << "memory exhausted: " << _size << "\n";
			exit(1);
		}
	}
private:
	mutable void* _seek;
	mutable void* _end;
	mutable size_t _size;
};
/*--------------------------------------------------------------------------*/
