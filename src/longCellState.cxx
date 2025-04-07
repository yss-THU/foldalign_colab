#ifndef LONGCELLSTATE
#define LONGCELLSTATE

#include "foldalign.hxx"
#include "mbllist.cxx"

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

class longCellState {
public:

	inline longCellState(): similarityScore(big_neg), energyScore(big_neg),
		 state(noState) {};
	inline longCellState(scoreType simSc, scoreType enSc, stateType st)
		:  similarityScore(simSc), energyScore(enSc), state(st) {};

	inline void set(scoreType simSc, scoreType enSc, stateType st, mbllist* mblptr = 0) {
		similarityScore = simSc;
		energyScore = enSc;
		state = st;
		//mblptr ignored
	};

	inline void setPointer(mbllist* mblptr = 0) {};

	inline scoreType getSimilarityScore() const {return similarityScore;}
	inline scoreType getEnergyScore() const {return energyScore;}
	inline stateType getState() const {return state;}
	inline mbllist* getPointer() const {return 0;}
	
private:
	scoreType similarityScore;    // The alignment score
	scoreType energyScore;
	stateType state;    // The state
};

#endif /*LONGCELLSTATE*/
