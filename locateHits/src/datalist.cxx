#ifndef DATALIST
#define DATALIST

#include <iostream>

#include "nohit.hxx"

#define TEMP_DEF_DATALIST template< class cell >
#define TEMP_SPES_DATALIST cell

/*! \brief A double linked list which is index after the alignments score. For
	each score a linked list of the alignments with that score is kept.
*/
TEMP_DEF_DATALIST
class datalist {

public:
	// Note checkScore is uninitialised before the add function is called for
	// the first time
	datalist() : startScores(0), topScore(0), checkScore(0), checkElement(0) {};
	
	/*! \brief Add a new alignment to the list
	*/
	void add(scoreType score, cell alignment);
	
	/*! \brief Get the next alignment for the list. The function starts by
		returning an alignment with the best score (if there are more than one
		it is not defined which one will be returned). If the there no more
		alignments with the current best score then the best score (as indicated
		by the topScore variable) is updated to the second best score, and an
		alignment with this score is returned.
	*/
	void getNext(scoreType& returnScore, cell& returnCell);

	/*! \brief Function unknown.
	*/
	void checkNext(scoreType& returnScore, cell& returnCell);
	
	/*! \brief Initialize the checkScore for use by the checkNext
	*/
	void initCheckNext() {
		checkScore = topScore; 
		if (topScore == 0) {checkElement = 0;}
		else {checkElement = topScore->start;}
	}
	
	~datalist();

private:
	
	//! \brief Stores a cell (an alignment) and a pointer to the next element
	struct element {
		// Pointer to the next element
		element* next;
		// Stores the data
		cell data;
	};

	/*! \brief Stores an element in the double link list. The list is sorted
		after alignment score
	*/
	struct store {
		//! The next node (score)
		store* next;
		//! The previous node (score)
		store* prev;
		//! The linked list of cells (alignments) with this score
		element* start;
		//! Pointer to the last cell (alignment) with this score
		element* last;
		//! The alignment score
		scoreType score;
	};
	
	//! The root of the double linked list (lowest score)
	store* startScores;
	
	//! The top of the double linked list (highest score)
	store* topScore;
	
	//! Used by the checkNext function
	store* checkScore;
	//! Used by the checkNext function
	element* checkElement;

	//! Stores an alignment
	void putStore(store*& curr, store*& next, store*& prev, element*& start,
                 element*& last, scoreType score);


};

//TEMP_DEF_DATALIST
//inline datalist< TEMP_SPES_DATALIST >::datalist() {
//
//	startScore = 0;
//}

TEMP_DEF_DATALIST
inline void datalist< TEMP_SPES_DATALIST >
::add(scoreType score, cell alignment) {

	store* curr_store = startScores;
	store* prev_store = 0;
		
	element* curr_elem = new element();
	curr_elem->next = 0;
	curr_elem->data = alignment;

	while(1) {
		
		if (curr_store == 0) {
			// We have reached the top.
			store* store_null = 0;
			if (prev_store != 0) {
				// There exists a previous node
				curr_store = new store();

				putStore(curr_store, store_null, prev_store, curr_elem, curr_elem, score);

				prev_store->next = curr_store;
			}
			else {
				// There isnt a previous node. This must be the start node
				curr_store = new store();

				putStore(curr_store, store_null, store_null, curr_elem, curr_elem, score);
				
				startScores = curr_store;
			}
			topScore = curr_store;
			checkScore = topScore;
			break;
		}
		else if (score > curr_store->score) {
			// The score is larger than the current score go to the next level
			prev_store = curr_store;
			curr_store = curr_store->next;
			// Do not break
		}
		else if (score == curr_store->score) {
			// This is the correct score. Store the alignment
			curr_store->last->next = curr_elem;
			curr_store->last = curr_elem;
			break;
		}
		else {
			// The score is smaller than the current score. 
			if (curr_store->prev != 0) {
				//The alignment is therefore stored between the previous and the 
				//current scores.
				prev_store = new store();

				putStore(prev_store, curr_store, curr_store->prev, curr_elem,
				         curr_elem, score);

				curr_store->prev = prev_store;
				prev_store->prev->next = prev_store;
				break;
			}
			else {
				prev_store = new store();

				putStore(prev_store, curr_store, curr_store->prev, curr_elem,
				         curr_elem, score);
				
				curr_store->prev = prev_store;
				startScores = prev_store;
				
				break;
			}
		}
	}
}

TEMP_DEF_DATALIST
inline void datalist< TEMP_SPES_DATALIST >
::putStore(store*& curr, store*& next, store*& prev, element*& start,
           element*& last, scoreType score) {

	curr->next = next;
	curr->prev = prev;
	curr->start = start;
	curr->last = last;
	curr->score = score;
	
}

TEMP_DEF_DATALIST
inline void datalist< TEMP_SPES_DATALIST >
::getNext(scoreType& returnScore, cell& returnCell) {

	if (topScore == 0) {returnScore = bad; return;}

	if (topScore->start == 0) {
		
		// Delete the old topScore
		store* tmpScore = topScore;
		topScore = topScore->prev;
		delete tmpScore;

		if (topScore == 0) {
		
			// There is nothing left. Return the bad score
			returnScore = bad;
			return;
			
		}
		else {getNext(returnScore, returnCell); return;}
	}
	else {
		returnScore = topScore->score;
		returnCell = topScore->start->data;
			
		// Delete the old start
		element* tmp = topScore->start;
		topScore->start = topScore->start->next;
		delete tmp;
		
		return;
	}
}	

TEMP_DEF_DATALIST
inline void datalist< TEMP_SPES_DATALIST >
::checkNext(scoreType& returnScore, cell& returnCell) {

	if (checkElement == 0) {
	
		if (checkScore == 0 || checkScore->prev == 0) {
			returnScore = bad;
			return;
		}
		checkScore = checkScore->prev;
		checkElement = checkScore->start;
		
		checkNext(returnScore, returnCell);
		return;
	}
	
	returnScore = checkScore->score;
	returnCell = checkElement->data;
	
	checkElement = checkElement->next;
}


TEMP_DEF_DATALIST
inline datalist< TEMP_SPES_DATALIST >::~datalist() {

	while (topScore != 0) {
	
		while(topScore->start != 0) {
		
			element* tmp = topScore->start;
			topScore->start = topScore->start->next;
			delete tmp;
			
		}
		
		store* tmpScore = topScore;
		topScore = topScore->prev;
		delete tmpScore;
	}
		
}

#endif /* DATALIST */
