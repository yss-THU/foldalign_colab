#ifndef STACK
#define STACK

#include "exception.cxx"

#include <string>
#include <iostream>

/******************************************************************************
*                                                                             *
*   Copyright 2004 - 2007 Jakob Hull Havgaard, hull@bioinf.ku.dk              *
*                                                                             *
*   This file is part of Foldalign                                            *
*                                                                             *
*   Foldalign is free software; you can redistribute it and/or modify         *
*   it under the terms of the GNU General Public License as published by      *
*   the Free Software Foundation; either version 2 of the License, or         *
*   (at your option) any later version.                                       *
*                                                                             *
*   Foldalign is distributed in the hope that it will be useful,              *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with Foldalign; if not, write to the Free Software                  *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA *
*                                                                             *
******************************************************************************/


// Written by Jakob Hull Havgaard, 2004, hull@bioinf.kvl.dk.
template<class item>
class stack {
// This is a stack of items

private:
	item* data; // The stack
	int max; // The size of the stack
	int next; // The next position to be used
	int prev; //Position fifo will read from
	stack(const stack&) {}; // No copy
	stack() {}; // No default

public:
	stack(int size) { // Init the stack
		max = size; // Save the size
		next = 0; // Init to first position
		prev = -1;
		data = new item[max]; // Allocate the storage
	};
	
	void push(const item& in) { // Store the in data
		if (next >= max) {
			std::string error = "Program error. Stack overflow";
			throw exception(error, false);
		}
		data[next] = in;
		next++;
	};
	
	item pop() { // Get the last data
		if (next < 1) {
			std::string error = "Program error. Stack underflow";
			throw exception(error, false);
		}
		next--;
		return(data[next]);
	};
	
	item fifo() {
		if (prev >= next) {
			std::string error = "Program error. Stack fifo overread";
			throw exception(error, false);
		}
		prev++;
		return(data[prev]);
	};
	
	int getsize() { // return the current stack size
		return(next);
	};

	int getfifosize() {
		return(next-prev-1);
	};

	~stack() {
		delete[] data;
	}
};
		
#endif /*STACK*/
