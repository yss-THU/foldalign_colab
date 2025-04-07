#ifndef MBLLIST
#define MBLLIST

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

class mbllist {
public:

	inline mbllist()
		: left_i(0), left_k(0), left_Wi(0), left_Wk(0),
		  right_i(0), right_k(0), right_Wi(0), right_Wk(0), 
		  simscore(big_neg), enscore(big_neg) {};
		  
	inline mbllist(positionType ri, positionType rk,
						lengthType rWi, lengthType rWk,
						positionType li, positionType lk,
						lengthType lWi, lengthType lWk,
						scoreType similarityScore, scoreType energyScore) :
			left_i(li), left_k(lk), left_Wi(lWi), left_Wk(lWk), 
			right_i(ri), right_k(rk), right_Wi(rWi), right_Wk(rWk), 
			simscore(similarityScore), enscore(energyScore) {};
	
	inline void get(positionType& ri, positionType& rk,
				  lengthType& rWi, lengthType& rWk, 
				  positionType& li, positionType& lk,
				  lengthType& lWi, lengthType& lWk, 
				  scoreType& similarityScore, scoreType& energyScore) {
		li = left_i; lk = left_k; lWi = left_Wi; lWk = left_Wk;
		ri = right_i; rk = right_k; rWi = right_Wi; rWk = right_Wk;
		similarityScore = simscore; energyScore = enscore;
	}
		
private:

	positionType left_i;
	positionType left_k;
	lengthType left_Wi;
	lengthType left_Wk;

	positionType right_i;
	positionType right_k;
	lengthType right_Wi;
	lengthType right_Wk;

	scoreType simscore;
	scoreType enscore;

};

#endif /* MBLLIST */
