#ifndef MULTISTACK
#define MULTISTACK
#include "stack.cxx"

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


// This class coordinates a number of stacks.
// The class is quite unsafe from over- and underflows
template< class first, class second>
class multistacks {
private:
	static const int n_stacks = 2;
	stack<first>*  f_stacks[n_stacks]; // The stacks
	stack<second>* s_stacks[n_stacks]; // The stacks
	int next; // The current stack size
	int prev; // The fifo read position
	
public:
	multistacks(int size) {
		next = 0; // The stacks are empty
		prev = -1;
		for(int n=0; n<n_stacks; n++) {
			f_stacks[n] = new stack<first>(size);
			s_stacks[n] = new stack<second>(size);
		}
	};
	
	void push(first i, first k, second Wi, second Wk) { // Store the values
		f_stacks[0]->push(i);
		f_stacks[1]->push(k);
		s_stacks[0]->push(Wi);
		s_stacks[1]->push(Wk);
		next++;
	};
	
	void pop(first& i, first& k, second& Wi, second& Wk) { // Get the values
		next--;
		i = f_stacks[0]->pop();
		k = f_stacks[1]->pop();
		Wi = s_stacks[0]->pop();
		Wk = s_stacks[1]->pop();
	};
	
	void fifo(first& i, first& k, second& Wi, second& Wk) { // Get the values
		prev++;
		i = f_stacks[0]->fifo();
		k = f_stacks[1]->fifo();
		Wi = s_stacks[0]->fifo();
		Wk = s_stacks[1]->fifo();
	};

	int getsize() { // Return the current stack size
		return(next);
	};

	int getfifosize() { // Return the current stack size
		return(next-prev-1);
	};

	~multistacks() {
		for(int n=0; n<n_stacks; n++) {
			delete f_stacks[n];
			delete s_stacks[n];
		}
	}

};
#endif /*MULTISTACK*/
