#ifndef CELL
#define CELL

#include "mbllist.cxx"
#include "foldalign.hxx"

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

// The memory cell class

class cell {
public:

//	inline cell(): score(big_neg), state(noState),
//								 length1(0), length2(0), length3(0), length4(0) {};

	inline cell(scoreType simSc = big_neg, scoreType enSc = big_neg,
				stateType st = noState, 
	            lengthType len1 = 0, lengthType len2 = 0,
				lengthType len3 = 0, lengthType len4 = 0)
			: similarityScore(simSc), energyScore(enSc), state(st),
			  length1(len1), length2(len2), 
		    length3(len3), length4(len4) {};

	inline void set(scoreType simSc, scoreType enSc, stateType st,
						 lengthType len1, lengthType len2,
						 lengthType len3, lengthType len4, mbllist* mblptr=0) {
		similarityScore = simSc;
		energyScore = enSc;
		state = st;
		length1 = len1;
		length2 = len2;
		length3 = len3;
		length4 = len4;
		// The mblptr is ignored on purpose
	};

	inline void setState(stateType st) {
		state = st;
	};

	// This function does nothing
	inline void setPointer(mbllist* mblpointer) {};

	inline scoreType getSimilarityScore() const {return similarityScore;};
	inline scoreType getEnergyScore() const {return energyScore;};
	inline stateType getState() const {return state;};
	inline void      getLengths(lengthType& len1, lengthType& len2,
	                            lengthType& len3, lengthType& len4) const {
		len1 = length1;
		len2 = length2;
		len3 = length3;
		len4 = length4;
	}

	inline lengthType getLength1() const {return length1;}
	inline lengthType getLength2() const {return length2;}
	inline lengthType getLength3() const {return length3;}
	inline lengthType getLength4() const {return length4;}
	
	inline mbllist* getPointer() const {return 0;}

private:
	scoreType similarityScore;    // The alignment similarity score
	scoreType energyScore;        // The structure energy score
	stateType state;    // The state
	lengthType length1; // Loop length sequence I left side
	lengthType length2; // Loop length sequence I right side
	lengthType length3; // Loop length sequence K left side
	lengthType length4; // Loop length sequence K right side
};

#endif /*CELL*/
