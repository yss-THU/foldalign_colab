#ifndef MBLCELL
#define MBLCELL

#include "cell.cxx"
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

class mblcell : public cell {
public:

	inline mblcell() : cell(), mblptr(0) {};
	inline mblcell(scoreType simSc, scoreType enSc, stateType st,
	              lengthType len1, lengthType len2,
					  lengthType len3, lengthType len4)
			: cell(simSc, enSc, st, len1, len2, len3, len4), mblptr(0) {};
	inline mblcell(scoreType simSc, scoreType enSc, stateType st,
	              lengthType len1, lengthType len2,
					  lengthType len3, lengthType len4,
					  mbllist* mblpointer) 
			: cell(simSc, enSc, st, len1, len2, len3, len4), mblptr(mblpointer) {};

	inline void set(scoreType simSc, scoreType enSc, stateType st,
						 lengthType len1, lengthType len2,
						 lengthType len3, lengthType len4, mbllist* mblpointer) {
		cell::set(simSc, enSc, st, len1, len2, len3, len4);
		mblptr = mblpointer;
	};
	
	inline void setPointer(mbllist* mblpointer) {mblptr = mblpointer;};

	inline mbllist* getPointer() const {return mblptr;};
	
private:

	// Pointer to the last multibranched loop this alignment has passed.
	mbllist* mblptr; 
	
};

#endif /* MBLCELL */
