#ifndef TWOLINK
#define TWOLINK

#include "longCell.cxx"
#include "foldalign.hxx"
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

#define TEMP_DEF_TWOLINK template< class storeType >
#define TEMP_SPES_TWOLINK storeType

template< class storeType >
class two_link {
public:
	
	inline two_link();
	
	inline storeType* putNext(const lengthType Wi, const lengthType Wk);
	
	inline storeType* getNext(lengthType& Wi, lengthType& Wk);
	
	inline storeType* getPos(lengthType& Wi, lengthType& Wk);

	inline void resetCurrent() {Wi_curr = Wi_start; Wk_curr = Wi_curr->start;}

	inline void lastWk();

	inline void removeRest();
	
	inline static long getSize() {
		return char_size*(Wi_size + Wk_size)/mem_scale;
	}

	inline ~two_link();

private:

	struct Wknode {
		Wknode* next;
		lengthType Wk;
		storeType cell;
	};

	struct Winode {
		Winode* next;
		Wknode* start;
		lengthType Wi;
	};

	Winode* Wi_curr;
	Winode* Wi_start;

	Wknode* Wk_curr;

	inline bool nextWi();

	// Deletes all Wk nodes for a Wi nodes
	inline void deleteWknodes(Winode* Wi_del);

	inline storeType* findWk(lengthType& Wk);

	static long Wi_size;
	static long Wk_size;
	
	const long Wi_cell_size;
	const long Wk_cell_size;

	static struct Wknode Wk_empty;

	void throwOutOfMemoryError();
};

TEMP_DEF_TWOLINK long two_link< TEMP_SPES_TWOLINK >::Wi_size = 0;
TEMP_DEF_TWOLINK long two_link< TEMP_SPES_TWOLINK >::Wk_size = 0;

TEMP_DEF_TWOLINK struct two_link< TEMP_SPES_TWOLINK >::Wknode 
two_link< TEMP_SPES_TWOLINK >::Wk_empty = {0, -1};

TEMP_DEF_TWOLINK
inline two_link< TEMP_SPES_TWOLINK >::two_link()
 : Wi_cell_size(sizeof(Winode)),
   Wk_cell_size(sizeof(Wknode)) {
	
	// Sounds like you've been snacking on some of the local mushrooms. 
	
	// Set the current node to a new node and make it the start node
	Wi_curr = new Winode;
	if ( Wi_curr == 0) {throwOutOfMemoryError();}
	Wi_start = Wi_curr;

	// This start node has no length
	Wi_curr->Wi = 0;
	
	// The start node is currently not followed by any node
	Wi_curr->next = 0;
	
	// The Wk node of the start node is set to a new node
	Wi_curr->start = new Wknode;
	if ( Wi_curr->start == 0) {throwOutOfMemoryError();}
	
	// The Wk_current is set to the Wk of the Wi start node
	Wk_curr = Wi_curr->start;
	
	// It has no length and the there is no next node, its cell has invalid values
	Wk_curr->Wk = -1;
	Wk_curr->next = 0;
	Wk_curr->cell.set(big_neg, big_neg, noState);
	
	Wi_size+=Wi_cell_size;
	Wk_size+=Wk_cell_size;
}

TEMP_DEF_TWOLINK 
inline two_link< TEMP_SPES_TWOLINK >::~two_link() {

	// Reset the current node to the start node and delete all elements in the
	// structure.

	resetCurrent();

	while (Wi_curr != 0) {

		deleteWknodes(Wi_curr);

		Winode* tmp = Wi_curr->next;
		delete Wi_curr;
		Wi_curr = tmp;
		Wi_size -= Wi_cell_size;
	}
}

TEMP_DEF_TWOLINK
inline storeType* two_link< TEMP_SPES_TWOLINK >::putNext(const lengthType Wi, const lengthType Wk) {

	if (Wi < Wi_curr->Wi) {
std::cerr << "Program Error. Warning: This part of the two_link code is not working correctly" << std::endl;
// Luckely this case should never happen.
// The problem is that the since Wi < Wi_curr->Wi the code below can never
// insert the data point at the right place. Currently the data point is added
// at the end of the datastructe which is guarenteed not to be the right place.
// Note that below it assumed that the Wks comes in the right order and so it
// is may be also reasonable to assume so for this for Wi
		while (Wi < Wi_curr->Wi) {
			if (Wi_curr->next == 0) {
				break;
			}
			Wi_curr = Wi_curr->next;
		}
		Wk_curr = Wi_curr->start;
		while (Wk < Wk_curr->Wk) {
			if (Wk_curr->next == 0) {
				break;
			}
			Wk_curr = Wk_curr->next;
		}
	}

	if (Wi == Wi_curr->Wi) {
		// The window size of the new element is the same as that of the current
		// Wi node. Therefor the Wi_node is the same and a new Wk_node is added
		Wk_curr->next = new Wknode;
		if (Wk_curr->next == 0) {
			throwOutOfMemoryError();
		}
		Wk_curr = Wk_curr->next;
	}
	else {
		// Add both a new Wi node and a new Wk node
		Wi_curr->next = new Winode;
		if (Wi_curr->next == 0) {
			throwOutOfMemoryError();
		}
		Wi_curr = Wi_curr->next;
		Wi_curr->next = 0;
		Wi_curr->Wi = Wi;
		Wi_curr->start = new Wknode;
		if (Wi_curr->start == 0) {
			throwOutOfMemoryError();
		}
		Wk_curr = Wi_curr->start;
		Wi_size += Wi_cell_size;
	}
	Wk_curr->next = 0;
	Wk_curr->Wk = Wk;
	Wk_size += Wk_cell_size;
	return &(Wk_curr->cell);
}

TEMP_DEF_TWOLINK
inline storeType* two_link< TEMP_SPES_TWOLINK >::getNext(lengthType& Wi, lengthType& Wk) {

	if (Wk_curr->next == 0) {
		
		// There are no more Wk nodes for the current Wi.
		// If possible go to the next Wi. If not possible
		// return invalid values
		
		if (!nextWi()) {

			// There are no more data in the structure
			// Return invalid values
			
			Wi = 0;
			Wk = 0;
			return 0;
		}

		// There is another Wi, The return element has been set by nextWi().

	}
	else {

		// There is still at least one Wk node left
		// Get the next Wk node

		Wk_curr = Wk_curr->next;
	
	}
		
	// Return the current element
	Wi = Wi_curr->Wi;
	Wk = Wk_curr->Wk;

	return &(Wk_curr->cell);
}

TEMP_DEF_TWOLINK
inline storeType* two_link< TEMP_SPES_TWOLINK >::getPos(lengthType& Wi, lengthType& Wk) {

	// Start at top of the data structure
	resetCurrent();

	if (Wi == Wi_curr->Wi) {return findWk(Wk);}

	// Get the next Wi.
	while (nextWi()) {

		// If the current Wi is to small -> get the next one.
		if (Wi_curr->Wi < Wi) {continue;}
		
		// Is the the correct Wi?
		if (Wi == Wi_curr->Wi) {return findWk(Wk);}

		// The current Wi is to large -> the Wi search for does not exist return 0

		return 0;
	}

	// The end of the data structure has been reached with no result. Return 0
	return 0;

}


TEMP_DEF_TWOLINK
inline storeType* two_link< TEMP_SPES_TWOLINK >::findWk(lengthType& Wk) {

	// Get the next Wk node as long as the current Wk is smaller than the
	// one needed

	while (Wk_curr->Wk < Wk) {

		// There are no more Wk_nodes return 0
		if (Wk_curr->next == 0) {return 0;}
				
		// Next node.
		Wk_curr = Wk_curr->next;

	}

	// If it is the right size return the cell
	if (Wk == Wk_curr->Wk) {return &(Wk_curr->cell);}

	// else return empty
	return 0;
}

TEMP_DEF_TWOLINK
inline bool two_link< TEMP_SPES_TWOLINK >::nextWi() {
	
	if (Wi_curr->next == 0) {return false;}
	
	Wi_curr = Wi_curr->next;
	Wk_curr = Wi_curr->start;

	return true;
}

TEMP_DEF_TWOLINK
inline void two_link< TEMP_SPES_TWOLINK >::lastWk() {
	
	// Fast forward to the end of the current Wk structure.
	// This will make the next call to getNext switch to the next Wi node
	Wk_curr = &Wk_empty;
}

TEMP_DEF_TWOLINK
inline void two_link< TEMP_SPES_TWOLINK >::removeRest() {

	// Remove the rest of the data structure starting from the next Wi node

	if (Wi_curr->next != 0) {

		// Save the location of the next Wi node
		Winode* Wi_tmp = Wi_curr->next;

		// Cut the link between the current node and the next node
		Wi_curr->next = 0;
		
		// Make the 'next' node the current node
		// This node and those down the chain are to be deleted
		Wi_curr = Wi_tmp;

		while(Wi_curr != 0) {

			deleteWknodes(Wi_curr);

			// Save the next Wi node
			Winode* tmp = Wi_curr->next;

			// Delete the current Wi node
			delete Wi_curr;

			// Make the saved next Wi node the current node
			Wi_curr = tmp;
			
			Wi_size -= Wi_cell_size;

		}
	}
	
	resetCurrent();
}

TEMP_DEF_TWOLINK
inline void two_link< TEMP_SPES_TWOLINK >::deleteWknodes(Winode* Wi_del) {

	// Initialize to the start node
	Wk_curr = Wi_del->start;

	// Keep going until the end of the data structure is found
	while (Wk_curr != 0) {

		// Save the location of the next node
		Wknode* tmp = Wk_curr->next;

		// Delete the current Wk node
		delete Wk_curr;
				
		// The new current node is the saved next node
		Wk_curr = tmp;
		
		Wk_size -= Wk_cell_size;
	}
	Wi_del->start = 0;

}

TEMP_DEF_TWOLINK
inline void two_link< TEMP_SPES_TWOLINK >::throwOutOfMemoryError() {
	const std::string error = "Memory error. Most likely out of memory";
	throw exception(error, false);
}
#endif /* TWOLINK */
