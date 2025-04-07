#ifndef MEM_STACK
#define MEM_STACK

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


template<class item, class count>
class mem_stack {

public:

	inline mem_stack(const count size) : max(size) {

		error = "Could not allocate the mem_stack array. Most likely cause: Out of memory.";
		if (max > 0) {
			try {
				data = new item*[max]; // Allocate the storage
			}
			catch ( ... ) {throw exception(error, false);}
			if (data == 0) {throw exception(error, false);}
		}
		
		error = "Mem_stack out of memory.";
			
		next = 0;
	};
	
	inline void push(item*& in) { // Store the in data

		if (next == max) {
			delete in;
		}
		else {
			data[next] = in;
			next++;
		}

	};
	
	inline item* pop() { // Get the last data

		if (next < 1) {
			item* returnVal = new item;
			if (returnVal == 0) {throw exception(error, false);}
			return returnVal;
		}

		next--;
		return(data[next]);
	};
	
	inline count getsize() const { // return the current stack size
		return(next);
	};

	inline ~mem_stack() {

		if (max > 0) {
			for(count i = 0; i < next; i++) {delete data[i];}
			delete[] data;
			data = 0;
		}
	}

private:
	item** data; // The stack
	const count max;  // The size of the stack
	count next; // The next position to be used
	mem_stack(const mem_stack&); // No copy
	mem_stack(); // No default
	mem_stack& operator=(const mem_stack& from); // No assignment
	std::string error; // The error message
};
		
#endif /*MEM_STACK*/
