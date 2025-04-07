#ifndef RESULTS
#define RESULTS

#include "foldalign.hxx"
#include "mbllist.cxx"

#include <iostream>
#include <math.h>

/******************************************************************************
*                                                                             *
*   Copyright 2004 - 2007 Jakob Hull Havgaard, hull@bioinf.u.dk             *
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

//************************************************************
// This class strores the results from a Foldalign alignment *
// Currently only the score and position of the best hit is  *
// stored                                                    *
// Started by Jakob Hull Havgaard 2004                       *
// hull@bioinf.kvl.dk                                        *
//************************************************************

// Oh, Coop, uh, about the uniform... replacing the quiet elegance of the dark
// suit and tie with the casual indifference of these muted earth tones is a
// form of fashion suicide, but, uh, call me crazy - on you it works. 

class results {
public:

	// Set up the result object (allocate memory for 4 position parameters)
//	results() : score(big_neg), logscore(big_neg), state(noState),
//		pos_i(0), pos_k(0), pos_Wi(0), pos_Wk(0), mblptr(0) {};

	results(scoreType new_score = big_neg, double new_logscore = big_neg,
				scoreType new_similarityScore = big_neg,
				scoreType new_energyScore = big_neg,
				stateType newstate = noState,
				positionType i = 0, positionType k = 0, 
				lengthType Wi = 0, lengthType Wk = 0, mbllist* mblpointer = 0)
		: score(new_score), logscore(new_logscore),
		  similarityScore(new_similarityScore), energyScore(new_energyScore),
		  state(newstate),
		  pos_i(i), pos_k(k), pos_Wi(Wi), pos_Wk(Wk), mblptr(mblpointer) {};

	// Return the score
	scoreType getScore() const {return score;};

	// Return the score
	double getLogScore() const {return logscore;};

	scoreType getSimilarityScore() const {return similarityScore;}

	scoreType getEnergyScore() const {return energyScore;}
	
	// Return the score
	scoreType getState() const {return state;};

	// return  mblpointer
	mbllist* getPointer() const {return mblptr;};

	// Return the positions through a reference
	void getPos(positionType& i, positionType& k, 
	            lengthType& Wi, lengthType& Wk ){
		i = pos_i; k = pos_k; Wi = pos_Wi; Wk = pos_Wk;
	};

	// Save a score and its position
	void store(scoreType new_score, double new_log_score,
			   scoreType new_similarityScore, scoreType new_energyScore,
			   stateType new_state,
	           positionType i, positionType k, 
				  lengthType Wi, lengthType Wk, mbllist* mblpointer = 0) {
		score = new_score;
		logscore = new_log_score;
		similarityScore = new_similarityScore;
		energyScore = new_energyScore;
		state = new_state;
		pos_i = i;
		pos_k = k;
		pos_Wi = Wi;
		pos_Wk = Wk;
		mblptr = mblpointer;
	};

private:
	scoreType score;     // Holds the best score
	double logscore;		// Holds the log of the score.
	scoreType similarityScore;
	scoreType energyScore;
	stateType state;     // Holds the best state
	positionType pos_i;  // Position i
	positionType pos_k;  // Position k
	lengthType pos_Wi;   // Position Wi
	lengthType pos_Wk;   // Position Wk
	mbllist* mblptr;		// The pointer to the last multibranch point.
};

#endif /*RESULTS*/
