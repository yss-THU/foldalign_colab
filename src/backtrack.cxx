#ifndef BACKTRACK
#define BACKTRACK
#include <iostream>
#include <math.h>
#include <new>
#include <queue>

#include "helper.cxx"
#include "output.cxx"
#include "fold.cxx"
#include "multistacks.cxx"
#include "arguments.cxx"
#include "results.cxx"
#include "cell.cxx"
#include "longCell.cxx"
#include "longCellState.cxx"
#include "longCellPtr.cxx"
#include "stack_ssl.cxx"
#include "mbllist.cxx"
#include "mblcell.cxx"
#include "constraints.cxx"
#include "tables.cxx"
#include "combineSimilarityEnergy.cxx"

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


class backtrack {
public:
	//*************************************
	// Realigns two subsequences and backtracks to get the structure
	//
	// The constructor initiates the first realignment and starts the following
	// steps. See note below

	inline backtrack(scoreType bestSimilarityScore,
					 scoreType bestEnergyScore,
	                 positionType begin_1,
						  positionType end_1,
						  positionType begin_2,
						  positionType end_2,
						  const sequence* const one,
						  const sequence* const two,
						  arguments& argu,
						  scorematrix& score,
						  output*& outp,
						  constraints* constr,
						  constraints* btConstr = 0);

	//*************************************
	// There is nothing to do for the destructor

//	inline ~backtrack() {};

private:

	// The backtrack/global alignment has 3 steps:
	//
	// 1: An almost normal branch alignment of the region is done. The only
	//    difference is that for coordinates a pointer to the last mbl is stored.
	//
	// 2: Using the extra mbl pointers all branch points in the final structure
	//    is located. This is used to split the backtrack into stems.
	//
	// 3: Each stem segment is realigned storing the all the information needed
	//    to backtrack the segment. The segment is then backtracked.


	struct alignPos {
		positionType i;
		positionType k;
		lengthType Wi;
		lengthType Wk;
		
		inline void set(positionType ni, positionType nk,
							 lengthType nWi, lengthType nWk) {
			i = ni; k = nk; Wi = nWi; Wk = nWk;
		}
	};
	
	struct mblPoint {
		alignPos begin;
		alignPos end;
		scoreType similarityScore;
		scoreType energyScore;
		bool branch;
		bool right;

		inline void set(alignPos nBegin, alignPos nEnd, 
						scoreType nSimilarityScore,
						scoreType nEnergyScore,
							 bool nBranch, bool nRight) {
			begin = nBegin;
			end = nEnd;
			similarityScore = nSimilarityScore;
			energyScore = nEnergyScore;
			branch = nBranch;
			right = nRight;
		}
	};
	
	struct alignmentScore {
		scoreType similarityScore;
		scoreType energyScore;
	};

	inline alignmentScore doIterativeAlignment(
					positionType begin_1, positionType end_1,
					positionType begin_2, positionType end_2, alignmentScore best_score);
	
	// Calls the fold object/function for realignment.
	template< class ltmClass, class startClass, class stmClass, 
	          bool global_align, bool realigning, bool mblrealign >
	inline alignmentScore foldfunc(results& reres,
						 		     longTermMemory< ltmClass >*& localLtm,
									  const positionType& i_dimension,
									  const positionType& k_dimension,
									  const positionType& sublen_I,
									  const positionType& sublen_K,
									  const bool& mem_global,
									  const positionType& i,
									  const positionType& j,
									  const positionType& k,
									  const positionType& l,
									  longTermMemory< startClass >*& startCoord,
									  stack_ssl<mbllist, int>*& mbl_stack);

	// Do the initial almost normal realignment
	template< class stmClass, class ltmClass, class startClass >
	inline alignmentScore mblRealign(results& reres,
						  longTermMemory< ltmClass >*& localLtm,
						  const positionType& i,
						  const positionType& j,
						  const positionType& k,
						  const positionType& l,
						  longTermMemory< startClass >*& startCoord,
						  stack_ssl<mbllist, int>*& mbl_stack);

	// After the initail realignment do a backtrack
	// which locates the stem segments
	// these can be realigned and backtracked as non-branched.
	inline void branchBt(longTermMemory< longCellPtr >*& ltm, 
											  std::queue<mblPoint>& branchPoints,
											  mbllist* brp,
											  positionType i,
											  positionType j,
											  positionType k,
											  positionType l);



	// Calculate the length of branching structure
	inline lengthType calcStructureLength(const positionType rightStart,
						const lengthType rightWindow,
						const positionType leftStart) const {
		return lengthType(rightStart + rightWindow - leftStart);
	};

	// A get the next branch point, used in branchBt()
	inline void getNextBranchPoint(mbllist*& brp,
						std::queue<mblPoint>& branchPoints,
						alignPos& begin, alignPos& end,
						alignPos& leftSide, alignPos& rightSide, mblPoint& mbl,
						scoreType& branchSimilarityScore,
						scoreType& branchEnergyScore, bool branch, bool right);
					
	inline void outputAlignment(positionType start1, positionType end1,
					positionType start2, positionType end2, alignmentScore best_score);

	// A clean up function for the fullBt function
	inline void fullBtCleanup(stack<positionType>*& leftPos_I,
							   stack<positionType>*& leftPos_K,
							   stack<positionType>*& leftBasepair_I,
							   stack<positionType>*& leftBasepair_K);

	//Iteratively backtracks from the end points through all seeds
	inline alignmentScore iterativeBacktrack(
								positionType begin_1, positionType end_1,
								positionType begin_2, positionType end_2);

	inline alignmentScore iterativeBacktrackBranch(
				positionType begin_1, positionType end_1,
				positionType begin_2, positionType end_2,
				std::queue<mblPoint>& branchPoints);

	inline alignmentScore iterativeBacktrackNoBranch(
				positionType begin_1, positionType end_1,
				positionType begin_2, positionType end_2,
				std::queue<mblPoint>& branchPoints);

	inline alignmentScore stemSegmentBacktrack(positionType begin_1,
					positionType begin_2, std::queue<mblPoint>& branchPoints);

	// Realign and backtrack a stem segment
	inline alignmentScore noBranchBt(alignPos& begin, alignPos& end,
										 scoreType startSimilarityScore,
										 scoreType startEnergyScore,
	                      		 stack<positionType>*& leftPos_I,
								 		 stack<positionType>*& leftBasepair_I,
								 		 stack<positionType>*& leftPos_K,
								 		 stack<positionType>*& leftBasepair_K,
								 		 stack<positionType>*& rightPos_I,
								 		 stack<positionType>*& rightBasepair_I,
								 		 stack<positionType>*& rightPos_K,
								 		 stack<positionType>*& rightBasepair_K);

	// Trace_state finds the previous state and position in
	// the dynamicprograming matrix
	inline stateType trace_state(positionType& i, positionType& k,
										  lengthType& Wi, lengthType& Wk,
										  stateType old_state,
										  longTermMemory< longCellState >*& localLtm,
	                   		     stack<positionType>*& leftPos_I,
										  stack<positionType>*& leftBasepair_I,
										  stack<positionType>*& leftPos_K,
										  stack<positionType>*& leftBasepair_K,
										  stack<positionType>*& rightPos_I,
										  stack<positionType>*& rightBasepair_I,
										  stack<positionType>*& rightPos_K,
										  stack<positionType>*& rightBasepair_K);

	inline void stacksPop(stack<positionType>*& rightPos_I,
						stack<positionType>*& rightPos_K,
						stack<positionType>*& rightBasepair_I,
						stack<positionType>*& rightBasepair_K,
						bool& right,
						positionType& position_i,
						positionType& position_k,
						stack< stack<positionType>* >* rightStackPos_I,
						stack< stack<positionType>* >* rightStackPos_K,
						stack< stack<positionType>* >* rightStackBasepair_I,
						stack< stack<positionType>* >* rightStackBasepair_K,
						stack< bool >* rightStackRigh,
						stack< positionType >* rightStackPosition_i,
						stack< positionType >* rightStackPosition_k);
						
	inline void stacksPush(stack<positionType>* rightPos_I,
						stack<positionType>* rightPos_K,
						stack<positionType>* rightBasepair_I,
						stack<positionType>* rightBasepair_K,
						bool right,
						positionType position_i,
						positionType position_k,
						stack< stack<positionType>* >* rightStackPos_I,
						stack< stack<positionType>* >* rightStackPos_K,
						stack< stack<positionType>* >* rightStackBasepair_I,
						stack< stack<positionType>* >* rightStackBasepair_K,
						stack< bool >* rightStackRight,
						stack< positionType >* rightStackPosition_i,
						stack< positionType >* rightStackPosition_k);
	
	// Empties the right stacks into the left stacks, and deletes the right
	// stacks
	inline void emptyRight(stack<positionType>*& leftPos_I,
								  stack<positionType>*& leftBasepair_I,
								  stack<positionType>*& leftPos_K,
								  stack<positionType>*& leftBasepair_K,
								  stack<positionType>*& rightPos_I,
								  stack<positionType>*& rightBasepair_I,
								  stack<positionType>*& rightPos_K,
								  stack<positionType>*& rightBasepair_K);

	// A monster function for calculating and printing the -backtrack_info
	inline void branch_info(longCellState*& cCell,
	                        const positionType& i, const positionType& k,
									const lengthType& Wi, const lengthType& Wk,
									const stateType oldState,
									longTermMemory< longCellState >*& localLtm);

	inline int makeSequences(char*& seq_1, char*& seq_2, char*& struc, 
									 int*& org_pos_1, int*& org_bp_1, int*& ali_bp_1,
									 stack<positionType>*& leftPos_1,
									 stack<positionType>*& leftBasepair_1,
									 int*& org_pos_2, int*& org_bp_2, int*& ali_bp_2,
									 stack<positionType>*& leftPos_2,
									 stack<positionType>*& leftBasepair_2);

	inline float identity(char sequence1[], char sequence2[], int len, float& similar, float& total);
	
	// Set the flow control arrays and the explanation array used by the
	// branch_info
	inline void set_explain();

	// Print the sorry no global alignment found message.
	inline void noGlobal();
	
											
	// Prints the coordinates and score of a stem
	inline void printStem(alignPos& b, alignPos& e, scoreType ss, scoreType es);

	inline lengthType size_k_dimension(const lengthType del) const 
		{return 4*del+3;}

	inline positionType calc_k_offset(const positionType begin_I,
										  const positionType begin_K,
										  const lengthType del) const
		{return begin_I - begin_K + 2*del+1;}

	//*************************************
	// The global variables

	const sequence* const seq_1; // Stores the sequences
	const sequence* const seq_2;

	arguments& arg;

	scorematrix& s_matrix; // Contains the score matrixs

	output* out; //Object containing the functions for outputing

	constraints* cons;
	constraints* backtrackConstraints;

	const lengthType delta; // Parameters
	const lengthType lambda;
	const lengthType min_loop;
	const bool global;
	const bool globalPrune;
	const bool plot_score;
	const bool back;
	const lengthType chunk_size;
	const std::string id;           // Alignment ID
	const std::string score_matrix; // The name of the score_matrix
	const bool nobranch;       // Nobranch (true) or branch (false) version
	const bool flip; // Is the order of the sequences switched
	const lengthType n_threads;

	const scoreType mblHelix;
	const scoreType mbl;
	const scoreType mblNuc;
	const scoreType gap_bonus;       // Affine gap bonus

	const positionType len_1; // The length of the sequences
	const positionType len_2;

	positionType end_1;
	positionType end_2;
	
	bool change_plot;
	lengthType stack_size;
	
	stack<positionType>* leftPos_I;
	stack<positionType>* leftPos_K;
	stack<positionType>* leftBasepair_I;
	stack<positionType>* leftBasepair_K;
	long branch_count;

	// The four loop lengths.
	// These are used during recalculation of the stm score from the ltm score
	// This only happens when the -backtrack_info option is used.
	lengthType len1;
	lengthType len2;
	lengthType len3;
	lengthType len4;

	// True if the non-GC cost score must be corrected
	// At the start of new branch (or the start of the backtrack)
	// all bWibWk and stem states have the correct non-GC cost added
	// but after the first base-pair (in either sequence) the cost must be
	// subtracted from score of alignments with theses states.
	// If this is not done correctly the score reported by -backtrack_info
	// will not be correct, everything else should be unaffected.
	bool gcCorrect_I;	
	bool gcCorrect_K;	
	
	bool gcCorrect[flow_size];
	bool hbi_loop[flow_size];
	std::string state_explain[flow_size];
	void (backtrack::*p2end[flow_size])(scoreType& score, const cell*& cCell);

	const tables< false, false > array;
	
};

//******************************************
// The constructor                         *
//******************************************


inline backtrack::backtrack(scoreType bestSimilarityScore, 
							scoreType bestEnergyScore,
									 positionType begin_1, 
									 positionType end_1,
									 positionType begin_2,
									 positionType end_2,
									 const sequence* const one,
									 const sequence* const two,
									 arguments& argu,
									 scorematrix& score,
									 output*& outp,
									 constraints* constr,
									 constraints* btConstr)
: seq_1(one),
  seq_2(two),
  arg(argu),
  s_matrix(score),
  out(outp),
  cons(constr),
  backtrackConstraints(btConstr),
  delta(arg.ltOpt("-max_diff")),
  lambda(arg.ltOpt("-max_length")),
  min_loop(arg.ltOpt("-min_loop")),
  global(arg.boolOpt("-global")),
  globalPrune(arg.boolOpt("-use_global_pruning")),
  plot_score(arg.boolOpt("-plot_score")),
  back(arg.boolOpt("-backtrack_info")),
  chunk_size(arg.ltOpt("-chunk_size")),
  id(arg.stringOpt("-ID")),
  score_matrix(arg.stringOpt("-score_matrix")),
  nobranch(arg.boolOpt("-nobranch")),
  flip(arg.boolOpt("switch")),
  n_threads(arg.ltOpt("-number_of_processors")),
  mblHelix(2*s_matrix.getMblAffine()),
  mbl(2*s_matrix.getMbl()),
  mblNuc(s_matrix.getMblNuc()),
  gap_bonus(s_matrix.getGap()),
  len_1(lengthType(end_1 - begin_1 +1)),
  len_2(lengthType(end_2 - begin_2 +1)),
  change_plot(false),
  stack_size(2*(len_2+1)+delta),
  branch_count(0)
   {




	//=======================================
	// Initialise objects and variabels

	arg.setBool("mblrealign", true);

	set_explain();

	// plot_score should be false unless this is global alignment
	if ( !global ) {
		if (plot_score) {
			arg.setBool("-plot_score",false);
			change_plot = true;
		}
		out->head();
	}
	else {
	
		if ( !plot_score ) {
			// Output the alignment head
			out->head();
		}

	}

	// A stem has a left side and a right side and maybe a loop connecting them.
	// The loop belongs in the left side.
	// For each position in the stem the nucleotide is stored in a 
	// <left, right>Pos stack. For each nucleotide the <left, right>Basepair
	// stack stores the position with which the nucleotide basepairs (-1 if the
	// nucleotide is not basepaired).
	
	// When the end of a stem is reached the right stacks are emptied into the
	// left stacks.
	// In the case of mbl the right stacks are put on a stack for later emptying.
	// These stacks are the rightStackPos and rightStackBasepair.

	leftPos_I = new stack<positionType>(stack_size);
	leftPos_K = new stack<positionType>(stack_size);
	leftBasepair_I = new stack<positionType>(stack_size);
	leftBasepair_K = new stack<positionType>(stack_size);
	
	if (leftPos_I == 0 || leftPos_K == 0 ||
	    leftBasepair_I == 0 || leftBasepair_K == 0) {
			std::string error = "Out of memory";
			throw exception(error, false);
	}


	alignmentScore best_score = {bestSimilarityScore, bestEnergyScore};
	best_score = doIterativeAlignment(begin_1, end_1, begin_2, end_2, best_score);

	outputAlignment(begin_1, end_1, begin_2, end_2, best_score);
	
	fullBtCleanup(leftPos_I, leftPos_K, leftBasepair_I, leftBasepair_K);
}

inline backtrack::alignmentScore backtrack::doIterativeAlignment(
					positionType begin_1, positionType end_1,
					positionType begin_2, positionType end_2, alignmentScore best_score) {


	alignmentScore alignScore = {big_neg, big_neg};
	try {
		alignScore = iterativeBacktrack(begin_1, end_1, begin_2, end_2);
		
		if (global) {
			best_score = alignScore;
		}
	}
	catch ( exception& exc ) {

		out->saveOutputError(arg.boolOpt("-plot_score"), exc.getMessage());

		fullBtCleanup(leftPos_I, leftPos_K, leftBasepair_I, leftBasepair_K);
		throw;
	}
	catch ( ... ) {

		out->saveOutputError(arg.boolOpt("-plot_score"));

		fullBtCleanup(leftPos_I, leftPos_K, leftBasepair_I, leftBasepair_K);
		throw;
	}

	return best_score;
}

inline backtrack::alignmentScore backtrack::iterativeBacktrack(
				positionType begin_1, positionType end_1,
				positionType begin_2, positionType end_2) {

	// This stack/data structure holds the information about the stem segments.
	std::queue<mblPoint> branchPoints;
	alignmentScore realign_score = {big_neg, big_neg};

	if (!nobranch) {
		realign_score = iterativeBacktrackBranch(begin_1, end_1, begin_2, end_2,
										branchPoints);
	}
	else {
		realign_score = iterativeBacktrackNoBranch(begin_1, end_1, begin_2, end_2,
										branchPoints);
	}

	return realign_score;
}

inline backtrack::alignmentScore backtrack::iterativeBacktrackBranch(
				positionType begin_1, positionType end_1,
				positionType begin_2, positionType end_2,
				std::queue<mblPoint>& branchPoints) {

	alignmentScore realign_score = {big_neg, big_neg};

	results reres;
	longTermMemory< longCellPtr >* ltm = 0;
	longTermMemory< longCellState >* startCoord = 0;
	stack_ssl<mbllist, int>* mbl_stack = 0;

	try {
//std::cout << "MBL REALIGN" << " " << begin_1 << " -> " << end_1 << " " << begin_2 << " -> " << end_2 << std::endl;
		realign_score = mblRealign< cell >(reres, ltm, begin_1, end_1, begin_2, end_2,
								        startCoord, mbl_stack);
//std::cout << "MBL REALIGN DONE ---------------------------------------------------" << std::endl;
	}
	catch ( ... ) {
		if (ltm != 0) {delete ltm; ltm = 0;}
		if (startCoord != 0) {delete startCoord; startCoord = 0;}
		if (mbl_stack != 0) {delete mbl_stack; mbl_stack = 0;}
			
		throw;
	}
	
	//====================================================================
	//
	//STEP 2: Backtrack using the branch point information to get a list
	// 		 of stem segments.
	//

	// If no alignment was found the score == big_neg and there is nothing
	// to backtrack

//std::cout << "BRANCH BT" << std::endl;
	if (realign_score.similarityScore != big_neg || realign_score.energyScore != big_neg) {
		mbllist* brp = reres.getPointer();
		branchBt(ltm, branchPoints, brp, begin_1, end_1, begin_2, end_2);
	}
//std::cout << "BRANCH BT DONE ................................................." << std::endl;
	// The reres, mbl_stack and the ltm is no longer needed.
	delete mbl_stack;
	delete ltm;

	mbl_stack = 0;
	ltm = 0;

	if (realign_score.similarityScore == big_neg || realign_score.energyScore == big_neg) {
		return realign_score;
	}

	//=========================================
	// Step 3
	//
	// Each stem segment is realigned and backtrack
	// A stem-segment may be traced into a seed constraint which must
	// then be fully backtrack.
	if (backtrackConstraints != 0) {
		cons = backtrackConstraints;
	}
//std::cout << "STEM SEGMENT BT" << std::endl;
	stemSegmentBacktrack(begin_1, begin_2, branchPoints);
//std::cout << "STEM SEGMENT BT DONE ==============================================" << std::endl;
	return realign_score;
}


inline backtrack::alignmentScore backtrack::iterativeBacktrackNoBranch(
				positionType begin_1, positionType end_1,
				positionType begin_2, positionType end_2,
				std::queue<mblPoint>& branchPoints) {

		//====================================================================
		//
		// Handle the non-branched case.
		// It is given that there is no branch points and hence no realignment
		// is needed to find them. This takes care of STEP 1 and 2.

	alignmentScore realign_score = {big_neg, big_neg};

	alignPos begin = {0, 0, 0, 0};
	alignPos end = {begin_1, begin_2, lengthType(end_1 - begin_1), lengthType(end_2 - begin_2)};
	
	mblPoint mbl = {begin, end, big_neg, false, false};

	branchPoints.push(mbl);

	if ( back ) {
		printStem(begin, end, big_neg, big_neg);
	}

	realign_score = stemSegmentBacktrack(begin_1, begin_2, branchPoints);

	return realign_score;
}

template< class stmClass, class ltmClass, class startClass >
inline backtrack::alignmentScore backtrack::mblRealign(results& reres,
						 	   longTermMemory< ltmClass >*& localLtm,
								const positionType& i,
								const positionType& j,
								const positionType& k,
								const positionType& l,
								longTermMemory< startClass >*& startCoord,
								stack_ssl<mbllist, int>*& mbl_stack) {


	const positionType i_dimension = j-i+1; //len_1;
	const positionType k_dimension = size_k_dimension(delta);

	// Allocate the mbl_stack
	// The mbl_stack is accessed through pointers which violates the
	// encapsulation for the mbl_stack object.

	mbl_stack = new stack_ssl<mbllist, int>();
	if (mbl_stack == 0) {
		std::string error = "Could not allocate the mbl stack. Most likely cause: Out of memory";
		throw exception(error, false);
	}

	// Do the realignment

	try {
		if (global) {
			foldfunc< ltmClass, startClass, mblcell, true, false, true >
			        (reres, localLtm, i_dimension, k_dimension, len_1, len_2, true,
						i, j, k, l, startCoord, mbl_stack);
		}
		else {
			foldfunc< ltmClass, startClass, mblcell, false, false, true >
		   	     (reres, localLtm, i_dimension, k_dimension, len_1, len_2, true,
						i, j, k, l, startCoord, mbl_stack);
		}
	}
	catch ( ... ) {
	
		if (mbl_stack != 0) {delete mbl_stack; mbl_stack = 0;}
		
		throw;
	}
	
	if (global && plot_score) {
		out->head();
		arg.setBool("-plot_score", false);
		change_plot = true;
	}

	if ( global ) {

		if (reres.getScore() == big_neg) {
			noGlobal();
			alignmentScore bad = {big_neg, big_neg};
			return bad;
		}

	}

	alignmentScore good = {reres.getSimilarityScore(), reres.getEnergyScore()};
	return good;
}

inline void backtrack::noGlobal() {
	std::cerr << "No good global alignment was found between sequences ";
	if ( ! flip ) {
		std::cerr  << seq_1->getName() << " and " << seq_2->getName() << "." << std::endl;
	}
	else {
		std::cerr  << seq_2->getName() << " and " << seq_1->getName() << "." << std::endl;
	}
	out->parameters();
	out->errorNoGlobal();
	out->saveOutput(arg.boolOpt("-plot_score"));
}

inline void backtrack::printStem(alignPos& b, alignPos& e, scoreType ss, scoreType es) {
	
	scoreType s = combineSimilarityEnergy(ss, es);
	out->backtrackPrintStem(e.i, e.k, e.Wi, e.Wk, b.i, b.k, b.Wi, b.Wk, s);
}


inline void backtrack::branchBt(longTermMemory< longCellPtr >*& ltm, 
										  std::queue<mblPoint>& branchPoints,
										  mbllist* brp,
										  positionType i,
										  positionType j,
										  positionType k,
										  positionType l) {

	//=====================================================================
	//
	// MBL backtrack
	//
	// The mbl backtrack aims at locating all branch points.
	// Each elemet in ltm has a pointer to an element (or 0) in mbl_stack.
	// This element is the mbl last passed by this alignment. The pointer 
	// from the best_score element leeds to the first branch point. The two
	// parts of this mbl each has a new pointer. These pointers are followed
	// to the next mbl. etc.

	stack_ssl<alignPos, int> rightSides;

	// The middel branch point coordinates and the alignment score
	scoreType branch_similarity_score = big_neg;
	scoreType branch_energy_score = big_neg;
	bool branch = false;
	bool right = false;
	
	alignPos end = {i, k, lengthType(j-i), lengthType(l-k)};

	alignPos begin = {0, 0, 0, 0};
	
	alignPos leftSide = {0, 0, 0, 0};
	alignPos rightSide = {0, 0, 0, 0};

	if (brp != 0) {
		brp->get(leftSide.i, leftSide.k, leftSide.Wi, leftSide.Wk, 
					rightSide.i, rightSide.k, rightSide.Wi, rightSide.Wk,
					branch_similarity_score, branch_energy_score);

		branch = true;

		lengthType leftWi = calcStructureLength(rightSide.i, rightSide.Wi, leftSide.i);
		lengthType leftWk = calcStructureLength(rightSide.k, rightSide.Wk, leftSide.k);

		begin.set(leftSide.i, leftSide.k, leftWi, leftWk);
	}
	
	mblPoint mbl = {begin, end, branch_similarity_score, branch_energy_score, branch, right};
	
	branchPoints.push(mbl);

	if ( back ) {
		printStem(begin, end, branch_similarity_score, branch_energy_score);
	}
	
	while (brp != 0) {

		if ( back ) {
			lengthType leftWi = calcStructureLength(rightSide.i, rightSide.Wi, leftSide.i);
			lengthType leftWk = calcStructureLength(rightSide.k, rightSide.Wk, leftSide.k);

			scoreType score = combineSimilarityEnergy(branch_similarity_score, branch_energy_score);
			out->backtrackPrintMBL(leftSide.i, leftSide.k, leftWi, leftWk,
								leftSide.i, leftSide.k, leftSide.Wi, leftSide.Wk,
								rightSide.i, rightSide.k, rightSide.Wi, rightSide.Wk,
								score);
		}

		if (rightSide.Wi != 0 && rightSide.Wk != 0) {
			rightSides.push(rightSide);
		
		
			end.set(leftSide.i, leftSide.k,
					  lengthType(rightSide.i -1 - leftSide.i),
					  lengthType(rightSide.k -1 - leftSide.k));

			// Check: Does the left side contain a branch point.
			brp = ltm->getPos(leftSide.i, leftSide.k, leftSide.Wi, leftSide.Wk)
			         ->getPointer();
		}
		else {
			brp = 0;
			end.set(leftSide.i, leftSide.k, leftSide.Wi, leftSide.Wk);
		}

		if (brp == 0) {

			// The left side does not contain a branch point (its a hairpin)
			// Put the left side into the branchPoints list and start processing
			// the stored right sides.
			
			// A hairpin has no begin coordinate.
			begin.set(0, 0, 0, 0);
			branch_similarity_score = big_neg;
			branch_energy_score = big_neg;
			
			mbl.set(begin, end, branch_similarity_score, branch_energy_score, false, false);
			
			// Putting it the left side away.
			branchPoints.push(mbl);

			if ( back ) {
				printStem(begin, end, branch_similarity_score, branch_energy_score);
			}

			while (rightSides.size() > 0) {

				// Process the stored right sides.

				end = rightSides.pop();
				
				brp = ltm->getPos(end.i, end.k, end.Wi, end.Wk)->getPointer();
				
				if (brp == 0) {
				
					// The store right side is a hairpin. Put it away and get the
					// next right side.

					mbl.set(begin, end, branch_similarity_score, branch_energy_score, false, true);				

					if ( back ) {printStem(begin, end, branch_similarity_score, branch_energy_score);}

					branchPoints.push(mbl);

				}
				else {

					// This is not a hairpin. Get the coordinates and give it the
					// full processing.
					getNextBranchPoint(brp, branchPoints, begin, end,
											leftSide, rightSide, mbl,
											branch_similarity_score, branch_energy_score, true, false);
					
					break;
				}
			}

		}
		else {
			getNextBranchPoint(brp, branchPoints, begin, end, leftSide, rightSide,
									mbl, branch_similarity_score, branch_energy_score, true, false);
		}
	}
		
}

inline void backtrack::getNextBranchPoint(mbllist*& brp,
					std::queue<mblPoint>& branchPoints,
					alignPos& begin, alignPos& end,
					alignPos& leftSide, alignPos& rightSide, mblPoint& mbl,
					scoreType& branchSimilarityScore,
					scoreType& branchEnergyScore, bool branch, bool right) {
					
		
	// Get the coordinates of the next branch point.

	brp->get(leftSide.i, leftSide.k, leftSide.Wi, leftSide.Wk,
				rightSide.i, rightSide.k, rightSide.Wi, rightSide.Wk, branchSimilarityScore, branchEnergyScore);
					
	lengthType leftWi = calcStructureLength(rightSide.i, rightSide.Wi, leftSide.i);
	lengthType leftWk = calcStructureLength(rightSide.k, rightSide.Wk, leftSide.k);

	begin.set(leftSide.i, leftSide.k, leftWi, leftWk);

	mbl.set(begin, end, branchSimilarityScore, branchEnergyScore, branch, right);				
		
	if ( back ) {
		printStem(begin, end, branchSimilarityScore, branchEnergyScore);
	}

	branchPoints.push(mbl);
}

inline void backtrack::outputAlignment(positionType start1, positionType end1,
					positionType start2, positionType end2, alignmentScore best_score) {

	//==================================================
	//
	// The backtrack is done.
	// Output the result
	//
	if ( !nobranch ) {arg.setBool("-nobranch", false);}
	if (change_plot) {arg.setBool("-plot_score", true);}
	

	if ( (best_score.similarityScore == big_neg ||
		  best_score.energyScore == big_neg) && nobranch ) {
		noGlobal();
	}
	else {
//		if ( nobranch ) { best_score = alignScore;}

		char* sequence_1 = new char[2*lambda]; // Stores the sequence
		char* sequence_2 = new char[2*lambda];
		char* structure = new char[2*lambda];  // Stores the structure of the first sequence
		int* org_pos_1 = new int[2*lambda];   // Position in original sequence
		int* org_pos_2 = new int[2*lambda];
		int* org_bp_1 = new int[2*lambda];    // Base-pairing in original coordinates
		int* org_bp_2 = new int[2*lambda];
		int* ali_bp_1 = new int[2*lambda];    // Base-pairing in the new coordinates
		int* ali_bp_2 = new int[2*lambda];

		for(int p=0; p < 2*lambda; p++) {structure[p] = '!';}

		// Make sequences and the structure
		int total_length = makeSequences(sequence_1, sequence_2, structure,
													org_pos_1, org_bp_1, ali_bp_1,
													leftPos_I, leftBasepair_I,
													org_pos_2, org_bp_2, ali_bp_2,
													leftPos_K, leftBasepair_K);

		float similar;
		float total;
		float sim = 100*identity(sequence_1, sequence_2, total_length, similar, 
		                         total);

		scoreType score = combineSimilarityEnergy(best_score.similarityScore, best_score.energyScore);

		out->out(sim, similar, total, sequence_1, sequence_2, structure, 
		        total_length, start1, end1, start2, end2,
				  score, org_pos_1, org_bp_1, ali_bp_1, org_pos_2, org_bp_2,
				  ali_bp_2);
				  
		// Clean up
		delete[] sequence_1;
		delete[] sequence_2;
		delete[] structure;
		delete[] org_pos_1;
		delete[] org_pos_2;
		delete[] org_bp_1;
		delete[] org_bp_2;
		delete[] ali_bp_1;
		delete[] ali_bp_2;
	}

}

inline backtrack::alignmentScore backtrack::stemSegmentBacktrack(positionType begin_1,
					positionType begin_2, std::queue<mblPoint>& branchPoints) {


	
	if ( !nobranch ) {
		arg.setBool("-nobranch", true);
	}

	arg.setBool("mblrealign", false);
	arg.setBool("realigning", true);
	
	
	positionType position_i = begin_1;
	positionType position_k = begin_2;
	
	alignmentScore alignScore = {big_neg, big_neg};

	stack<positionType>* rightPos_I;
	stack<positionType>* rightPos_K;
	stack<positionType>* rightBasepair_I;
	stack<positionType>* rightBasepair_K;

	stack< stack<positionType>* >* rightStackPos_I = new stack<stack<positionType>* >(stack_size);
	stack< stack<positionType>* >* rightStackPos_K = new stack<stack<positionType>* >(stack_size);
	stack< stack<positionType>* >* rightStackBasepair_I = new stack<stack<positionType>* >(stack_size);
	stack< stack<positionType>* >* rightStackBasepair_K = new stack<stack<positionType>* >(stack_size);
	stack< bool >* rightStackRight = new stack< bool >(stack_size);
	stack< positionType >* rightStackPosition_i = new stack< positionType >(stack_size);
	stack< positionType >* rightStackPosition_k = new stack< positionType >(stack_size);
	
	while (!branchPoints.empty()) {
	
		mblPoint mbl = branchPoints.front();
		branchPoints.pop();
		
		alignPos begin = mbl.begin;
		alignPos end   = mbl.end;

		scoreType branch_similarity_score = mbl.similarityScore;
		scoreType branch_energy_score = mbl.energyScore;
		bool right = mbl.right;
		
		while (end.i > position_i+1 || end.k > position_k+1) {

			// This is the end of a right side of a branch point.
				
			// Get the next right stacks from the stacks of stacks
			stacksPop(rightPos_I, rightPos_K, rightBasepair_I, rightBasepair_K,
						right, position_i, position_k,
						rightStackPos_I, rightStackPos_K,
						rightStackBasepair_I, rightStackBasepair_K, rightStackRight,
						rightStackPosition_i, rightStackPosition_k);
						
			emptyRight(leftPos_I, leftBasepair_I, leftPos_K, leftBasepair_K,
						  rightPos_I, rightBasepair_I, rightPos_K, rightBasepair_K);

			if (back) {
				out->backtrackEnd();
			}
		}			

		rightPos_I = new stack<positionType>(stack_size);
		rightPos_K = new stack<positionType>(stack_size);
		rightBasepair_I = new stack<positionType>(stack_size);
		rightBasepair_K = new stack<positionType>(stack_size);

		if (rightPos_I == 0 || rightPos_K == 0 || 
		    rightBasepair_I == 0 || rightBasepair_K == 0) {
			 	std::string error = "Out of memory";
				throw exception(error, false);
		}

		alignScore = noBranchBt(begin, end, branch_similarity_score, 
								branch_energy_score,
					   				leftPos_I, leftBasepair_I,
										leftPos_K, leftBasepair_K,
					   				rightPos_I, rightBasepair_I,
										rightPos_K, rightBasepair_K);

		if (begin.i == 0 && begin.k == 0 && begin.Wi == 0 && begin.Wk == 0) {

			// This stem does not start at a branch point
	
			// Empty the right stacks into the left stacks, and delete the right
			// stacks
			emptyRight(leftPos_I, leftBasepair_I, leftPos_K, leftBasepair_K,
						  rightPos_I, rightBasepair_I, rightPos_K, rightBasepair_K);

			position_i = end.i + end.Wi;
			position_k = end.k + end.Wk;

			if (back) {
				out->backtrackEnd();
			}
		}
		else {
			
			// Branch point. Store the right stacks in a stack of stacks
			stacksPush(rightPos_I, rightPos_K, rightBasepair_I, rightBasepair_K,
					right, end.i+end.Wi, end.k+end.Wk,
					rightStackPos_I, rightStackPos_K,
					rightStackBasepair_I, rightStackBasepair_K, rightStackRight,
					rightStackPosition_i, rightStackPosition_k);
			
			position_i = begin.i;
			position_k = begin.k;

			if ( back ) {
				branch_count++;
				out->backtrackStart(branch_count);
			}
		}
	}

	while ( rightStackBasepair_I->getsize() > 0 ) {

		// Make sure that the right basepair stacks are all empty

		// Get the next right stacks from the stacks of stacks
		bool right; // just there so stacksPop() can be used
		stacksPop(rightPos_I, rightPos_K, rightBasepair_I, rightBasepair_K,
					right, position_i, position_k,
					rightStackPos_I, rightStackPos_K,
					rightStackBasepair_I, rightStackBasepair_K, rightStackRight,
					rightStackPosition_i, rightStackPosition_k);

		// Empty the right stacks into the left stacks and delete the right
		// stacks
		emptyRight(leftPos_I, leftBasepair_I, leftPos_K, leftBasepair_K,
					  rightPos_I, rightBasepair_I, rightPos_K, rightBasepair_K);
	}
	
	if ( rightStackBasepair_K->getsize() > 0 ) {
		// This should never be true.
		std::string error = 
		   "Program error: The right basepair stack for sequence K is not empty";
		throw exception(error, false);
	}
	
	delete rightStackPos_I;
	delete rightStackPos_K;
	delete rightStackBasepair_I;
	delete rightStackBasepair_K;
	delete rightStackRight;
	delete rightStackPosition_i;
	delete rightStackPosition_k;

	if (nobranch) {
//		realign_score = alignScore;
	}
	return alignScore;


}

inline void backtrack::stacksPush(stack<positionType>* rightPos_I,
						stack<positionType>* rightPos_K,
						stack<positionType>* rightBasepair_I,
						stack<positionType>* rightBasepair_K,
						bool right,
						positionType position_i,
						positionType position_k,
						stack< stack<positionType>* >* rightStackPos_I,
						stack< stack<positionType>* >* rightStackPos_K,
						stack< stack<positionType>* >* rightStackBasepair_I,
						stack< stack<positionType>* >* rightStackBasepair_K,
						stack< bool >* rightStackRight,
						stack< positionType >* rightStackPosition_i,
						stack< positionType >* rightStackPosition_k) {


	rightStackPos_I->push(rightPos_I);
	rightStackPos_K->push(rightPos_K);
	rightStackBasepair_I->push(rightBasepair_I);
	rightStackBasepair_K->push(rightBasepair_K);
	rightStackRight->push(right);
	rightStackPosition_i->push(position_i);
	rightStackPosition_k->push(position_k);

}

inline void backtrack::stacksPop(stack<positionType>*& rightPos_I,
						stack<positionType>*& rightPos_K,
						stack<positionType>*& rightBasepair_I,
						stack<positionType>*& rightBasepair_K,
						bool& right,
						positionType& position_i,
						positionType& position_k,
						stack< stack<positionType>* >* rightStackPos_I,
						stack< stack<positionType>* >* rightStackPos_K,
						stack< stack<positionType>* >* rightStackBasepair_I,
						stack< stack<positionType>* >* rightStackBasepair_K,
						stack< bool >* rightStackRight,
						stack< positionType >* rightStackPosition_i,
						stack< positionType >* rightStackPosition_k) {


	rightPos_I = rightStackPos_I->pop();
	rightPos_K = rightStackPos_K->pop();
	rightBasepair_I = rightStackBasepair_I->pop();
	rightBasepair_K = rightStackBasepair_K->pop();
	right = rightStackRight->pop();
	position_i = rightStackPosition_i->pop();
	position_k = rightStackPosition_k->pop();

}

inline void backtrack::fullBtCleanup(stack<positionType>*& leftPos_I,
							   stack<positionType>*& leftPos_K,
							   stack<positionType>*& leftBasepair_I,
							   stack<positionType>*& leftBasepair_K) {
	delete leftPos_I;
	delete leftPos_K;
	delete leftBasepair_I;
	delete leftBasepair_K;

	leftPos_I = 0;
	leftPos_K = 0;
	leftBasepair_I = 0;
	leftBasepair_K = 0;
	
}

inline void backtrack::emptyRight(stack<positionType>*& leftPos_I,
								 			 stack<positionType>*& leftBasepair_I,
											 stack<positionType>*& leftPos_K,
								 			 stack<positionType>*& leftBasepair_K,
								 			 stack<positionType>*& rightPos_I,
								 			 stack<positionType>*& rightBasepair_I,
								 			 stack<positionType>*& rightPos_K,
								 			 stack<positionType>*& rightBasepair_K) {

	// Empty the right stacks into the left stacks
	while(rightPos_I->getsize() > 0) {
		leftPos_I->push(rightPos_I->pop());
		leftPos_K->push(rightPos_K->pop());
		leftBasepair_I->push(rightBasepair_I->pop());
		leftBasepair_K->push(rightBasepair_K->pop());
	}

	// Delete the right stacks
	delete rightPos_I;
	delete rightPos_K;
	delete rightBasepair_I;
	delete rightBasepair_K;

	rightPos_I = 0;
	rightPos_K = 0;
	rightBasepair_I = 0;
	rightBasepair_K = 0;
}


template< class ltmClass, class startClass, class stmClass, 
			 bool global_align, bool realigning, bool mblrealign >
inline backtrack::alignmentScore backtrack::foldfunc(results& reres,
										 	 	 longTermMemory< ltmClass >*& localLtm,
										 		 const positionType& i_dimension,
											 	 const positionType& k_dimension,
										 		 const positionType& sublen_I,
										 		 const positionType& sublen_K,
										 		 const bool& mem_global,
										 		 const positionType& i,
										 		 const positionType& j,
												 const positionType& k,
												 const positionType& l,
										 		 longTermMemory< startClass >*& startCoord,
										 		 stack_ssl<mbllist, int>*& mbl_stack) {


	//*********************************************************************
	//
	// Call the fold algorithm
	//
	// The reres and localLtm objects must be delete by the caller. Also in
	// case of exceptions.


	// Allocate result object
	std::string error = "Could not allocate the backtrack result object. Most likely cause: Out of memory.";
	reres.store(big_neg, double(big_neg), big_neg, big_neg, noState, 0, 0, 0, 0, 0);

	// Allocate ltm
	error = "Could not allocate long term memory (during backtrack) Most likely cause: Out of memory.";
	const positionType k_offset = calc_k_offset(i,k,delta);

	localLtm = new longTermMemory< ltmClass >(k_dimension, i_dimension, sublen_I,
												 sublen_K, mem_global, n_threads, k_offset);
	if ( localLtm == 0) {throw exception(error, false);}
	
	if (startCoord != 0 ) {
		arg.setBool("hpstart", false);
	}
	else {
		arg.setBool("hpstart", true);
	}
	

//std::cout << "FOLDING: " <<i << " " << j << " " << k << " " << l << " " << global_align<< " " << realigning << " " << mblrealign<< std::endl;

	try {
		if (n_threads == 1) {
			if (globalPrune) {
				fold<stmClass, ltmClass, startClass, global_align, true, 
						realigning, mblrealign, false>
					RNA(i, j, k, l, reres, seq_1, seq_2, arg, s_matrix,
						localLtm, out, false, cons, startCoord, mbl_stack);
			}
			else {
				fold<stmClass, ltmClass, startClass, global_align, false, 
						realigning, mblrealign, false>
					RNA(i, j, k, l, reres, seq_1, seq_2, arg, s_matrix,
						localLtm, out, false, cons, startCoord, mbl_stack);
			}
		}
		else {
			if (globalPrune) {
				fold<stmClass, ltmClass, startClass, global_align, true, 
						realigning, mblrealign, true>
					RNA(i, j, k, l, reres, seq_1, seq_2, arg, s_matrix,
						localLtm, out, false, cons, startCoord, mbl_stack);
			}
			else {
				fold<stmClass, ltmClass, startClass, global_align, false, 
						realigning, mblrealign, true>
					RNA(i, j, k, l, reres, seq_1, seq_2, arg, s_matrix,
						localLtm, out, false, cons, startCoord, mbl_stack);
			}
		}
	}
	catch ( std::bad_alloc& bad ) {

		std::string error = "Memory error ";
		error += bad.what();

		out->saveOutputError(arg.boolOpt("-plot_score"), error);
	
		throw exception(error, false);
	}
	catch ( exception& exc ) {

		out->saveOutputError(arg.boolOpt("-plot_score"), exc.getMessage());

		throw;
	}
	catch ( ... ) {

		out->saveOutputError(arg.boolOpt("-plot_score"));

		std::string error = "Unknown fold error. Could not align sequences.";
		throw exception(error, false);
	}
//std::cout << "FOLDING DONE" << std::endl;

	alignmentScore good = {reres.getSimilarityScore(), reres.getEnergyScore()};
	return good;
}

inline backtrack::alignmentScore backtrack::noBranchBt(alignPos& begin, alignPos& end,
											 scoreType startSimilarityScore,
											 scoreType startEnergyScore,
											 stack<positionType>*& leftPos_I,
											 stack<positionType>*& leftBasepair_I,
											 stack<positionType>*& leftPos_K,
											 stack<positionType>*& leftBasepair_K,
											 stack<positionType>*& rightPos_I,
											 stack<positionType>*& rightBasepair_I,
											 stack<positionType>*& rightPos_K,
											 stack<positionType>*& rightBasepair_K) {

	//=================================================================
	//
	// Realign and backtrack a stem segment.
	//
	// Currently the full stem segment is realigned and backtracked, but
	// as the segment is known to non-branch it should be easy to make
	// a qubic version of a linear space model known from sequence alignment

	results reres;
	longTermMemory< longCellState >* localLtm;
	stack_ssl<mbllist, int>* mbl_stack = 0;
	
	const lengthType sublen_I = end.Wi - begin.Wi +1;
	const lengthType sublen_K = end.Wk - begin.Wk +1;

	const positionType i_dimension = sublen_I;
	const positionType k_dimension = size_k_dimension(delta);

	const positionType k_offset = calc_k_offset(end.i, end.k, delta);

	gcCorrect_I = true;
	gcCorrect_K = true;

//	scoreType alignSimilarityScore = startSimilarityScore;
//	scoreType alignEnergyScore = startEnergyScore;

	// Allocate the start coordinates object
	arg.setBool("hpstart", true);

	constraints* save_cons = cons;
	
//	bool checkSeed = true;
	
	longTermMemory< longCellState >* startCoord = 0;
	if (begin.i != 0 || begin.k != 0) {

		startCoord = new longTermMemory< longCellState >(k_dimension, i_dimension,
											 sublen_I, sublen_K, true, n_threads, k_offset,
											 begin.i);
		if ( startCoord == 0) {
			std::string error = "Could not allocate start coordinates object (during backtrack) Most likely cause: Out of memory.";

			throw exception(error, false);
		}
		startCoord->setPositions(end.i, end.k);

		longCellState* startpos = startCoord->putPos(begin.i, begin.k, 
													begin.Wi, begin.Wk);
		startpos->set(startSimilarityScore, startEnergyScore, mblIK);

		arg.setBool("hpstart", false);
		cons = 0;
//		checkSeed = false;
	}	
	//===========================================
	//
	// Do the stem segment realignment

	alignmentScore alignScore;
	try {
		if (global) {
			alignScore = foldfunc< longCellState, longCellState, cell, true, true, false >
									(reres, localLtm, i_dimension, k_dimension,
									 sublen_I, sublen_K, true,
									 end.i, positionType(end.i+end.Wi),
									 end.k, positionType(end.k+end.Wk),
									 startCoord, mbl_stack);
		}
		else {
			alignScore = foldfunc< longCellState, longCellState, cell, false, true, false >
									(reres, localLtm, i_dimension, k_dimension,
									 sublen_I, sublen_K, true,
									 end.i, positionType(end.i+end.Wi),
									 end.k, positionType(end.k+end.Wk),
									 startCoord, mbl_stack);
		}
	}
	catch ( ... ) {
		if (localLtm != 0) {delete localLtm; localLtm = 0;}
		throw;
	}

	cons = save_cons;
	
	stateType state = reres.getState();

	if (startCoord != 0) {
		delete startCoord;
		startCoord = 0;
	}
	
	if (state == noState) {

		// This may be due to a branched alignment and nobrached alignments can
		// give slightly different results. The stem segment realignment is also
		// nobranched while the other alignments may be branched.

		std::string error = "Could not realign stem segment";
		
		delete localLtm;
		
		throw exception(error, false);
	}

	state = noState;

	//===========================================
	//
	// Backtrack the segment
	//
	
	if ( (!nobranch) || alignScore.similarityScore != big_neg ||
		alignScore.energyScore != big_neg) {

		positionType i = end.i;
		positionType k = end.k;
		lengthType  Wi = end.Wi;
		lengthType  Wk = end.Wk;

		while (state != hp_init) {

			if (array.isThisASeedConstraint(state)) {
				
				delete localLtm;
				localLtm = 0;
			
				if ( !nobranch ) {arg.setBool("-nobranch", false);}

				alignScore = iterativeBacktrack(i, i+Wi, k, k+Wk);

				return alignScore;
			}

			state = trace_state(i, k, Wi, Wk, state, localLtm,
									  leftPos_I, leftBasepair_I, leftPos_K, leftBasepair_K,
						   		  rightPos_I, rightBasepair_I, rightPos_K, rightBasepair_K);


			if ((state == noState) || (state >= flow_size)) {

				std::string error = "Program error! Illegal state found during backtrack. State: ";
				error += int(state);
				delete localLtm;
				throw exception(error, false);
			}

		}
	}

	delete localLtm;
	localLtm = 0;

	return alignScore;
}

inline stateType backtrack::trace_state(positionType& i,
                                        positionType& k,
													 lengthType& Wi,
													 lengthType& Wk,
													 stateType old_state,
													 longTermMemory< longCellState >*& localLtm,
													 stack<positionType>*& leftPos_I,
													 stack<positionType>*& leftBasepair_I,
													 stack<positionType>*& leftPos_K,
													 stack<positionType>*& leftBasepair_K,
													 stack<positionType>*& rightPos_I,
													 stack<positionType>*& rightBasepair_I,
													 stack<positionType>*& rightPos_K,
													 stack<positionType>*& rightBasepair_K) {

	// This function does most of the backtrack

	// He is Bob, eager for fun. When he wears a smile, everybody run! 

	longCellState* branch = localLtm->getPos(i, k, Wi, Wk);

	if (branch == 0) {
		std::cerr << "Found empty cell (branch) during traceback (trace_state)"<< std::endl;
		helper::pcerr("Previous cell:", i,k,Wi,Wk);
		std::cerr << "state: " << int(old_state) << " " << int(array.convertSeed2state(old_state)) << std::endl;
		return noState;
	}
	

	if (back) {
		// Outputs backtrack information if wanted
		branch_info(branch, i, k, Wi, Wk, array.convertSeed2state(old_state),
					localLtm);
	}
	// You punks owe me ten grand! Leo needs a new pair of shoes! 

	const stateType state = branch->getState();
	const stateType noSeedState = array.convertSeed2state(state);

	positionType ii=i, kk=k; // The next position
	lengthType Wii=Wi, Wkk=Wk; // More indicies to use

	switch (noSeedState) {
		case noState:
			std::cerr << "Program error. Warning state noState found" << std::endl;
			ii = i+1; Wii = Wi-2; kk=k+1; Wkk=Wk-2; // The new indicies
			leftPos_I->push(i);
			leftPos_K->push(k);
			leftBasepair_I->push(i+Wi);
			leftBasepair_K->push(k+Wk);
			rightPos_I->push(i+Wi);
			rightPos_K->push(k+Wk);
			rightBasepair_I->push(i);
			rightBasepair_K->push(k);
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return stem_IK;

		// Align ikWiWk. With a potential base-pair
		case hp_pb_IK:
		case bi_pb_I_bk_pb_K:
		case bWi_pb_I_bWk_pb_K:
		case il_pb_I_il_pb_K:
		case mbl_il_pb_I_mbl_il_pb_K:
			ii = i+1; Wii = Wi-2; kk=k+1; Wkk=Wk-2; // The new indicies
			leftPos_I->push(i);
			rightPos_I->push(i+Wi);
			leftPos_K->push(k);
			rightPos_K->push(k+Wk);
			if ((old_state == stem_IK) ||
			    (old_state == stem_I_stem_gap_kWk) || (old_state == stem_no_mbl_I_stem_gap_kWk) ||
				 (old_state == stem_gap_iWi_stem_K) || (old_state == stem_no_mbl_gap_iWi_stem_K))
			{
				leftBasepair_K->push(k+Wk);
				rightBasepair_K->push(k);
				leftBasepair_I->push(i+Wi);
				rightBasepair_I->push(i);
			}
			else {
				leftBasepair_K->push(-1);
				rightBasepair_K->push(-1);
				leftBasepair_I->push(-1);
				rightBasepair_I->push(-1);
			}
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return state;
		
		// Align ikWiWk with basepair
		case stem_IK:
			ii = i+1; Wii = Wi-2; kk=k+1; Wkk=Wk-2; // The new indicies
			leftPos_I->push(i);
			leftPos_K->push(k);
			leftBasepair_I->push(i+Wi);
			leftBasepair_K->push(k+Wk);
			rightPos_I->push(i+Wi);
			rightPos_K->push(k+Wk);
			rightBasepair_I->push(i);
			rightBasepair_K->push(k);
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return state;

		// Align iWi to gaps. Basepair iWi.
		case stem_I_stem_gap_kWk:
		case stem_no_mbl_I_stem_gap_kWk:
			ii = i+1; Wii = Wi-2; kk=k; Wkk=Wk; // The new indicies
			leftPos_I->push(i);
			rightPos_I->push(i+Wi);
			leftPos_K->push(0);
			rightPos_K->push(0);
			leftBasepair_K->push(-1);
			rightBasepair_K->push(-1);
			leftBasepair_I->push(i+Wi);
			rightBasepair_I->push(i);
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return state;

		// Align kWk to gaps. Basepair kWk
		case stem_gap_iWi_stem_K:
		case stem_no_mbl_gap_iWi_stem_K:
			ii = i; Wii = Wi; kk=k+1; Wkk=Wk-2; // The new indicies
			leftPos_I->push(0);
			rightPos_I->push(0);
			leftPos_K->push(k);
			rightPos_K->push(k+Wk);
			leftBasepair_K->push(k+Wk);
			rightBasepair_K->push(k);
			leftBasepair_I->push(-1);
			rightBasepair_I->push(-1);
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return state;

		// Align i and k
		case hp_init:
		case hp_init_align_ik:
		case hp_init_gap_Wi_ik:
		case hp_init_gap_Wk_ik:
		case bi_I_bk_K:
		case bi_I_bk_gap_Wk:
		case bi_gap_Wi_bk_K:
		case il_I_il_K_ik:
		case il_I_il_gap_Wk_ik:
		case il_gap_Wi_il_K_ik:
		case mbl_il_I_mbl_il_K_ik:
		case mbl_il_I_mbl_il_gap_Wk_ik:
		case mbl_il_gap_Wi_mbl_il_K_ik:
			ii = i+1; Wii = Wi-1; kk=k+1; Wkk=Wk-1; // The new indicies
			leftPos_I->push(i);
			leftPos_K->push(k);
			leftBasepair_I->push(-1);
			leftBasepair_K->push(-1);
			len1--; len3--;
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return state;

		// Align Wi and Wk
		case hp_init_align_WiWk:
		case hp_init_gap_I_WiWk:
		case hp_init_gap_K_WiWk:
		case bWi_I_bWk_K:
		case bWi_I_bWk_gap_k:
		case bWi_gap_i_bWk_K:
		case il_I_il_K_WiWk:
		case il_I_il_gap_k_WiWk:
		case il_gap_i_il_K_WiWk:
		case mbl_bWi_I_mbl_bWk_K:
		case mbl_il_I_mbl_il_K_WiWk:
		case mbl_il_I_mbl_il_gap_k_WiWk:
		case mbl_il_gap_i_mbl_il_K_WiWk:
			ii = i; Wii = Wi-1; kk=k; Wkk=Wk-1; // The new indicies
			rightPos_I->push(i+Wi);
			rightPos_K->push(k+Wk);
			rightBasepair_I->push(-1);
			rightBasepair_K->push(-1);
			len2--; len4--;
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return state;

		// Align i to a gap
		case hp_init_gap_K_i:
		case hp_init_gap_kWk_i:
		case hp_init_gap_kWi_i:
		case bi_I_bk_gap_k:
		case bi_I_bk_gap_kWk:
		case bi_gap_Wi_bk_gap_k:
		case il_I_il_gap_k_i:
		case il_I_il_gap_kWk_i:
		case il_gap_Wi_il_gap_k_i:
		case mbl_il_I_mbl_il_gap_k_i:
		case mbl_il_I_mbl_il_gap_kWk_i:
		case mbl_il_gap_Wi_mbl_il_gap_k_i:
			ii = i+1; Wii = Wi-1; kk=k; Wkk=Wk; // The new indicies
			leftPos_I->push(i);
			leftPos_K->push(0);
			leftBasepair_I->push(-1);
			leftBasepair_K->push(-1);
			len1--;
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return state;

		// Align Wi to a gap
		case hp_init_gap_Wk_Wi:
		case hp_init_gap_kWk_Wi:
		case hp_init_gap_iWk_Wi:
		case bWi_I_bWk_gap_kWk:
		case bWi_I_bWk_gap_Wk:
		case bWi_gap_i_bWk_gap_Wk:
		case il_I_il_gap_Wk_Wi:
		case il_I_il_gap_kWk_Wi:
		case il_gap_i_il_gap_Wk_Wi:
		case mbl_bWi_I_mbl_bWk_gap_Wk:
		case mbl_il_I_mbl_il_gap_kWk_Wi:
		case mbl_il_I_mbl_il_gap_Wk_Wi:
		case mbl_il_gap_i_mbl_il_gap_Wk_Wi:
			ii = i; Wii = Wi-1; kk=k; Wkk=Wk; // The new indicies
			rightPos_I->push(i+Wi);
			rightPos_K->push(0);
			rightBasepair_I->push(-1);
			rightBasepair_K->push(-1);
			len2--;
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return state;

		// Align k to a gap
		case hp_init_gap_I_k:
		case hp_init_gap_iWi_k:
		case hp_init_gap_iWk_k:
		case bi_gap_i_bk_K:
		case bi_gap_i_bk_gap_Wk:
		case bi_gap_iWi_bk_K:
		case il_gap_iWi_il_K_k:
		case il_gap_i_il_K_k:
		case il_gap_i_il_gap_Wk_k:
		case mbl_il_gap_iWi_mbl_il_K_k:
		case mbl_il_gap_i_mbl_il_K_k:
		case mbl_il_gap_i_mbl_il_gap_Wk_k:
			ii = i; Wii = Wi; kk=k+1; Wkk=Wk-1; // The new indicies
			leftPos_I->push(0);
			leftPos_K->push(k);
			leftBasepair_I->push(-1);
			leftBasepair_K->push(-1);
			len3--;
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return state;

		// Align Wk to a gap
		case hp_init_gap_Wi_Wk:
		case hp_init_gap_iWi_Wk:
		case hp_init_gap_kWi_Wk:
		case bWi_gap_iWi_bWk_K:
		case bWi_gap_Wi_bWk_K:
		case bWi_gap_Wi_bWk_gap_k:
		case il_gap_iWi_il_K_Wk:
		case il_gap_Wi_il_K_Wk:
		case il_gap_Wi_il_gap_k_Wk:
		case mbl_bWi_gap_Wi_mbl_bWk_K:
		case mbl_il_gap_Wi_mbl_il_K_Wk:
		case mbl_il_gap_Wi_mbl_il_gap_k_Wk:
		case mbl_il_gap_iWi_mbl_il_K_Wk:
			ii = i; Wii = Wi; kk=k; Wkk=Wk-1; // The new indicies
			rightPos_I->push(0);
			rightPos_K->push(k+Wk);
			rightBasepair_I->push(-1);
			rightBasepair_K->push(-1);
			len4--;
			i = ii; Wi = Wii; k = kk; Wk = Wkk;
			return state;

		// Branch point
		case mblIK:
			if (state == mblIK) {
				return hp_init;
			}
			return state;

		// The default state should never be reached
		default:
			stateType dummy = int(state);
			std::string error = "Program error. Unknown state ";
			error += int(dummy);
			error += " (";
			error += int(i);
			error += ",";
			error += int(Wi);
			error += ") (";
			error += int(k);
			error += ",";
			error += int(Wk);
			error + ")";
			throw exception(error, false);
			
	}
}

inline int backtrack::makeSequences(char*& sequence_1, char*& sequence_2, 
										char*& struc, 
               					int*& org_pos_1, int*& org_bp_1, int*& ali_bp_1, 
									   stack<positionType>*& leftPos_1,
				   					stack<positionType>*& leftBasepair_1,
               					int*& org_pos_2, int*& org_bp_2, int*& ali_bp_2,
				   					stack<positionType>*& leftPos_2,
				   					stack<positionType>*& leftBasepair_2)
{
	// Get a sequence and a structure from the backtrack stacks.
	// Returns the length of the sequence
	int counter = 0;
	lengthType length = 2*lambda;
	stack<positionType> bp_1(length);
	stack<positionType> bp_2(length);

	while (leftPos_1->getfifosize() > 0) {

		// Get the position
		int pos_1 = leftPos_1->fifo();
		int pos_2 = leftPos_2->fifo();
		int basepair_1 = leftBasepair_1->fifo();
		int basepair_2 = leftBasepair_2->fifo();

		// Build the sequences
		if (pos_1 > 0) { // Handles nucleotides (not gaps)
			sequence_1[counter] = toupper(seq_1->getLetter(pos_1));
			org_pos_1[counter] = pos_1;
		}
		else { // Handle gaps
			sequence_1[counter] = s_matrix.getLetter(0);
			org_pos_1[counter] = 0;
		}

		if (pos_2 > 0) { // Handles nucleotides (not gaps)
			sequence_2[counter] = toupper(seq_2->getLetter(pos_2));
			org_pos_2[counter] = pos_2;
		}
		else { // Handle gaps
			sequence_2[counter] = s_matrix.getLetter(0);
			org_pos_2[counter] = 0;
		}

		// Get the base-pairing information
		
		if ((basepair_1 > -1) && (basepair_2 > -1)) {
			if (pos_1 < basepair_1) {
				struc[counter] = open_bp;
				bp_1.push(counter);
				bp_2.push(counter);
			}
			else {
				struc[counter] = close_bp;
				int left_1 = bp_1.pop();
				int left_2 = bp_2.pop();
				ali_bp_1[counter] = left_1;
				ali_bp_2[counter] = left_2;
				ali_bp_1[left_1] = counter;
				ali_bp_2[left_2] = counter;
			}
			org_bp_1[counter] = basepair_1;
			org_bp_2[counter] = basepair_2;
		}
		else if ((basepair_1 > -1) || (basepair_2 > -1)) {
		
			if (basepair_1 > -1) {
				org_bp_1[counter] = basepair_1;
				org_bp_2[counter] = 0;
				ali_bp_2[counter] = -1;
				if (pos_1 < basepair_1) {
					struc[counter] = open_one_seq;
					bp_1.push(counter);
				}
				else {
					struc[counter] = close_one_seq;
					int left_1 = bp_1.pop();
					ali_bp_1[counter] = left_1;
					ali_bp_1[left_1] = counter;
				}
			}
			else {
				org_bp_1[counter] = 0;
				org_bp_2[counter] = basepair_2;
				ali_bp_1[counter] = -1;
				if (pos_2 < basepair_2) {
					struc[counter] = open_one_seq;
					bp_2.push(counter);
				}
				else {
					struc[counter] = close_one_seq;
					int left_2 = bp_2.pop();
					ali_bp_2[counter] = left_2;
					ali_bp_2[left_2] = counter;
				}
			}
		}		
		else {
			struc[counter] = unpaired;
			org_bp_1[counter] = 0;
			org_bp_2[counter] = 0;
			ali_bp_1[counter] = -1;
			ali_bp_2[counter] = -1;
		}
		counter++;
	}
	sequence_1[counter] = '\0';
	sequence_2[counter] = '\0';
	struc[counter] = '\0';
	return counter;
}

inline float backtrack::identity(char sequence1[], char sequence2[], int len, float& similar, float& total) {
	total = 0; // The total number of nucleotides. Two aligned gaps are not counted
	similar = 0; // The number of similar nucleotides. Gaps do not count as similar

	for(int i=0; i < len; i++) {
		if ((sequence1[i] == '\0') || (sequence2[i] == '\0')) {
			if (sequence1[i] != sequence2[i]) {
				std::cerr << "Program error. The two sequences in the identity function do not have the same length" << std::endl;
				throw -1;
			}
			return similar/total;
		}
		if (sequence1[i] == sequence2[i]) {
			// The nucleotides are similar
			if (sequence1[i] != '-') {
				// They are not gaps
				total++;
				similar++;
			}
		}
		else {
			total++;
		}
	}
	return similar/total;
}


inline void backtrack::branch_info(longCellState*& cCell, const positionType& i, const positionType& k, const lengthType& Wi, const lengthType& Wk, const stateType oldState, longTermMemory< longCellState >*& localLtm) {

	// A simple output function gone bad.
	
	// Handles the output for the -backtrack_info option


	scoreType similarityScore = cCell->getSimilarityScore();
	scoreType energyScore = cCell->getEnergyScore();
	const stateType state = array.convertSeed2state(cCell->getState());
	
	// The stm score is the score where hairpin bulge and internal loops are not
	// treated as mbl loops.
	scoreType stm_score = energyScore;

	//*******************************************************************
	// Calculate the short term memory score

	if ( (! hbi_loop[oldState]) &&  hbi_loop[state]) {

		// The current state is the first of a series of non-mbl loop nucleotides
		// The four loop lengths are calculated.
		// They are used to calculate the score as if it was a non-mbl loop score 
		// than a mbl score loop scores.
		// The lengths are found be backtracking.
		
		stateType tmpState = state;

		positionType local_i = i;
		positionType local_k = k;
		lengthType local_Wi = Wi;
		lengthType local_Wk = Wk;
		len1 = len2 = len3 = len4 = 0;

		while ( hbi_loop[tmpState] ) {
			
			switch (tmpState) {

				// i & k
				case hp_init:
				case hp_init_align_ik:
				case hp_init_gap_Wi_ik:
				case hp_init_gap_Wk_ik:
				case bi_I_bk_K:
				case bi_I_bk_gap_Wk:
				case bi_gap_Wi_bk_K:
				case il_I_il_K_ik:
				case il_I_il_gap_Wk_ik:
				case il_gap_Wi_il_K_ik:
					local_i++;
					local_k++;
					local_Wi--;
					local_Wk--;
					len1++;
					len3++;
					break;

				// Align Wi and Wk
				case hp_init_align_WiWk:
				case hp_init_gap_I_WiWk:
				case hp_init_gap_K_WiWk:
				case bWi_I_bWk_K:
				case bWi_I_bWk_gap_k:
				case bWi_gap_i_bWk_K:
				case il_I_il_K_WiWk:
				case il_I_il_gap_k_WiWk:
				case il_gap_i_il_K_WiWk:
					local_Wi--;
					local_Wk--;
					len2++;
					len4++;
					break;

				// Align i to a gap
				case hp_init_gap_K_i:
				case hp_init_gap_kWk_i:
				case hp_init_gap_kWi_i:
				case bi_I_bk_gap_k:
				case bi_I_bk_gap_kWk:
				case bi_gap_Wi_bk_gap_k:
				case il_I_il_gap_k_i:
				case il_I_il_gap_kWk_i:
				case il_gap_Wi_il_gap_k_i:
					local_i++;
					local_Wi--;
					len1++;
					break;

				// Align Wi to a gap
				case hp_init_gap_Wk_Wi:
				case hp_init_gap_kWk_Wi:
				case hp_init_gap_iWk_Wi:
				case bWi_I_bWk_gap_kWk:
				case bWi_I_bWk_gap_Wk:
				case bWi_gap_i_bWk_gap_Wk:
				case il_I_il_gap_Wk_Wi:
				case il_I_il_gap_kWk_Wi:
				case il_gap_i_il_gap_Wk_Wi:
					local_Wi--;
					len2++;
					break;

				// Align k to a gap
				case hp_init_gap_I_k:
				case hp_init_gap_iWi_k:
				case hp_init_gap_iWk_k:
				case bi_gap_i_bk_K:
				case bi_gap_i_bk_gap_Wk:
				case bi_gap_iWi_bk_K:
				case il_gap_iWi_il_K_k:
				case il_gap_i_il_K_k:
				case il_gap_i_il_gap_Wk_k:
					local_k++;
					local_Wk--;
					len3++;
					break;

				// Align Wk to a gap
				case hp_init_gap_Wi_Wk:
				case hp_init_gap_iWi_Wk:
				case hp_init_gap_kWi_Wk:
				case bWi_gap_iWi_bWk_K:
				case bWi_gap_Wi_bWk_K:
				case bWi_gap_Wi_bWk_gap_k:
				case il_gap_iWi_il_K_Wk:
				case il_gap_Wi_il_K_Wk:
				case il_gap_Wi_il_gap_k_Wk:
					local_Wk--;
					len4++;
					break;

				default:
					std::string error = "Program error! Unknown state found during local backtrack. State: ";
					error += tmpState;
					throw exception(error, false);
			}

			if (tmpState == hp_init) {break;}
			longCellState* branch = localLtm->getPos(local_i, local_k, local_Wi, local_Wk);
			tmpState = branch->getState();
		}
	}
	
	if ( hbi_loop[state] ) {

		switch (state) {
		
			case hp_init:
			case hp_init_align_ik:
			case hp_init_align_WiWk:
			case hp_init_gap_I_k:
			case hp_init_gap_I_WiWk:
			case hp_init_gap_K_i:
			case hp_init_gap_K_WiWk:
			case hp_init_gap_Wi_ik:
			case hp_init_gap_Wi_Wk:
			case hp_init_gap_Wk_ik:
			case hp_init_gap_Wk_Wi:
			case hp_init_gap_iWi_k:
			case hp_init_gap_iWi_Wk:
			case hp_init_gap_kWk_i:
			case hp_init_gap_kWk_Wi:
			case hp_init_gap_iWk_k:
			case hp_init_gap_iWk_Wi:
			case hp_init_gap_kWi_i:
			case hp_init_gap_kWi_Wk:
				stm_score += s_matrix.getHpLength(len1, len3);
				stm_score -= mblNuc * (len1 + len2 + len3 + len4);
				break;
				
			case bi_I_bk_K:
			case bi_I_bk_gap_k:
			case bi_I_bk_gap_kWk:
			case bi_I_bk_gap_Wk:
			case bi_gap_i_bk_K:
			case bi_gap_i_bk_gap_Wk:
			case bi_gap_iWi_bk_K:
			case bi_gap_Wi_bk_K:
			case bi_gap_Wi_bk_gap_k:
				stm_score += s_matrix.getBulgeLength(len1, len3);
				stm_score -= mblNuc * (len1 + len2 + len3 + len4);
				break;

			case bWi_I_bWk_K:
			case bWi_I_bWk_gap_k:
			case bWi_I_bWk_gap_kWk:
			case bWi_I_bWk_gap_Wk:
			case bWi_gap_i_bWk_K:
			case bWi_gap_i_bWk_gap_Wk:
			case bWi_gap_iWi_bWk_K:
			case bWi_gap_Wi_bWk_K:
			case bWi_gap_Wi_bWk_gap_k:
				stm_score += s_matrix.getBulgeLength(len2, len4);
				stm_score -= mblNuc * (len1 + len2 + len3 + len4);
				break;

			case il_I_il_K_ik:
			case il_I_il_K_WiWk:
			case il_I_il_gap_k_i:
			case il_I_il_gap_k_WiWk:
			case il_I_il_gap_Wk_ik:
			case il_I_il_gap_Wk_Wi:
			case il_I_il_gap_kWk_i:
			case il_I_il_gap_kWk_Wi:
			case il_gap_iWi_il_K_k:
			case il_gap_iWi_il_K_Wk:
			case il_gap_i_il_K_k:
			case il_gap_i_il_K_WiWk:
			case il_gap_i_il_gap_Wk_k:
			case il_gap_i_il_gap_Wk_Wi:
			case il_gap_Wi_il_K_ik:
			case il_gap_Wi_il_K_Wk:
			case il_gap_Wi_il_gap_k_i:
			case il_gap_Wi_il_gap_k_Wk:
				stm_score += s_matrix.getIntLoopLength(len1, len2);
				stm_score += s_matrix.getIntLoopLength(len3, len4);
				stm_score -= mblNuc * (len1 + len2 + len3 + len4);
				break;

			default:
				std::string error = "Program error! Illegal state found during backtrack rescoring. State: ";
				error += char(int(state));
				throw exception(error, false);
		}
	}			
	//*******************************************************************
	// The cases of mbl and long bulges the longTerm Memory score must be 
	// corrected for the nonGC close cost.

	// When a new stem is started the nonGC stem end score is correct. The
	// gcCorrect_I and gcCorrect_K are set to true.
	// The first time a base-pair is passed the score is correct but all other
	// times it must be correct. After this the gcCorrect is set to false


	if ( gcCorrect[state] ) {
	
		// The owls are not what they seem.
		
		if ( (gcCorrect_I) && (gcCorrect_K) ) {

			// The score is correct.
			
			// If this is bulge state then noting happens
			
			if ( array.get_right_branch(state) ) {

				// It is a stem case change either both or one of the sides
				// are no longer correct

				if (stem_IK == state) {gcCorrect_I = false; gcCorrect_K = false;}
				else if  (stem_I_stem_gap_kWk == state) {gcCorrect_I = false;}
				else {gcCorrect_K = false;}
				
			}
		}
		else if ( !gcCorrect_I && gcCorrect_K ) {
		
			energyScore -= s_matrix.getNonGCEnd(seq_1->getPos(i), seq_1->getPos(i+Wi));
			gcCorrect_K = false;
		}
		else if ( !gcCorrect_K && gcCorrect_I ) {

			energyScore -= s_matrix.getNonGCEnd(seq_2->getPos(k), seq_2->getPos(k+Wk));
			gcCorrect_I = false;
		}
		else {
			
			lengthType len2 = 0;
			lengthType len4 = 0;

			if ( ! array.get_right_branch(state) ) {

				// sometimes my arms bend back...

				// The bulge states must be backtracked to a stem state

				stateType tmpState = state;

				while ( !array.get_right_branch(tmpState) ) {
				
					switch (tmpState) {
						case bWi_I_bWk_K:
						case bWi_I_bWk_gap_k:
						case bWi_gap_i_bWk_K:
							len2++;
							len4++;
							break;
						case bWi_I_bWk_gap_kWk:
						case bWi_I_bWk_gap_Wk:
						case bWi_gap_i_bWk_gap_Wk:
							len2++;
							break;
						case bWi_gap_iWi_bWk_K:
						case bWi_gap_Wi_bWk_K:
						case bWi_gap_Wi_bWk_gap_k:
							len4++;
							break;
						default:
							std::string error = "Program error! Illegal state found during branch_info backtrack: ";
							error += int(tmpState);
							throw exception(error, false);
					}
					lengthType nWi = Wi - len2;
					lengthType nWk = Wk - len4;
					tmpState = localLtm->getPos(i, k, nWi, nWk)->getState();
					
				}
			}

			energyScore -= s_matrix.getNonGCEnd(seq_1->getPos(i), seq_1->getPos(i+Wi-len2));
			energyScore -= s_matrix.getNonGCEnd(seq_2->getPos(k), seq_2->getPos(k+Wk-len4));
		
		}
	}

	// End of longTerm memorys GC correction
	//******************************************************************

	if ( ! hbi_loop[state] ) {len1 = len2 = len3 = len4 = 0;}

	scoreType score = combineSimilarityEnergy(similarityScore, energyScore);
	stm_score  = combineSimilarityEnergy(similarityScore, stm_score);
	
	out->backtrack(score, stm_score, i, k, Wi, Wk, state, len1, len2, len3, len4, state_explain[state]);
}


//************************************************************
// An initialization function                                *
//************************************************************

inline void backtrack::set_explain() {

	// This function sets the human readable explanations for the states.
	// If the input sequences are switched then some of the explanations must
	// also be switched. This is also handled by this function.

	// Other array initializations below


	std::string error_msg = "Program error. Illegal state.";

	helper::init_array(state_explain, flow_size, error_msg);
	helper::init_array(gcCorrect, flow_size, false);
	helper::init_array(hbi_loop, flow_size, false);

	state_explain[noState]                       = "Program error! This is the default no state.";
	state_explain[hp_init]                       = "Initial Hairpin ik";
	state_explain[hp_init_align_ik]              = "Hairpin ik";
	state_explain[hp_init_align_WiWk]            = "Hairpin jl";
	state_explain[hp_init_gap_I_k]               = "Hairpin k gap i";
	state_explain[hp_init_gap_I_WiWk]            = "Hairpin jl gap i";
	state_explain[hp_init_gap_K_i]               = "Hairpin i gap k";
	state_explain[hp_init_gap_K_WiWk]            = "Hairpin jl gap k";
	state_explain[hp_init_gap_Wi_ik]             = "Hairpin ik gap j";
	state_explain[hp_init_gap_Wi_Wk]             = "Hairpin l gap j";
	state_explain[hp_init_gap_Wk_ik]             = "Hairpin ik gap l";
	state_explain[hp_init_gap_Wk_Wi]             = "Hairpin j gap l";
	state_explain[hp_init_gap_iWi_k]             = "Hairpin k gap i & j";
	state_explain[hp_init_gap_iWi_Wk]            = "Hairpin l gap i & j";
	state_explain[hp_init_gap_kWk_i]             = "Hairpin i gap k & l";
	state_explain[hp_init_gap_kWk_Wi]            = "Hairpin j gap k & l";
	state_explain[hp_init_gap_iWk_k]             = "Hairpin k gap i & l";
	state_explain[hp_init_gap_iWk_Wi]            = "Haiprin j gap i & l";
	state_explain[hp_init_gap_kWi_i]             = "Hairpin i gap k & j";
	state_explain[hp_init_gap_kWi_Wk]            = "Hairpin l gap k & j";
	state_explain[hp_pb_IK]                      = "Hairpin -> stem ik jl";
	state_explain[stem_IK]                       = "Basepair ik jl";
	state_explain[stem_I_stem_gap_kWk]           = "Basepair i j gap k & l";
	state_explain[stem_no_mbl_I_stem_gap_kWk]    = "Basepair i j gap k & l k & l must bp later";
	state_explain[stem_gap_iWi_stem_K]           = "Basepair k l gap i & j";
	state_explain[stem_no_mbl_gap_iWi_stem_K]    = "Basepair k l gap i & j i & j must bp later";
	state_explain[bi_I_bk_K]                     = "Bulge ik";
	state_explain[bi_I_bk_gap_k]                 = "Bulge i gap k";
	state_explain[bi_I_bk_gap_kWk]               = "Bulge i gap k & l";
	state_explain[bi_I_bk_gap_Wk]                = "Bulge ik gap l";
	state_explain[bi_pb_I_bk_pb_K]               = "Bulge -> stem ik jl";
	state_explain[bi_gap_i_bk_K]                 = "Bulge k gap i";
	state_explain[bi_gap_i_bk_gap_Wk]            = "Bulge k gap i & l";
	state_explain[bi_gap_iWi_bk_K]               = "Bulge k gap i & j";
	state_explain[bi_gap_Wi_bk_K]                = "Bulge ik gap j ";
	state_explain[bi_gap_Wi_bk_gap_k]            = "Bulge i gap k & j";
	state_explain[bWi_I_bWk_K]                   = "Bulge jl";
	state_explain[bWi_I_bWk_gap_k]               = "Bulge jl gap k";
	state_explain[bWi_I_bWk_gap_kWk]             = "Bulge j gap k & l";
	state_explain[bWi_I_bWk_gap_Wk]              = "Bulge j gap l";
	state_explain[bWi_pb_I_bWk_pb_K]             = "Bulge -> stem ij kl";
	state_explain[bWi_gap_i_bWk_K]               = "Bulge il gap i";
	state_explain[bWi_gap_i_bWk_gap_Wk]          = "Bulge i gap i";
	state_explain[bWi_gap_iWi_bWk_K]             = "Bulge l gap i & j";
	state_explain[bWi_gap_Wi_bWk_K]              = "Bulge l gap j";
	state_explain[bWi_gap_Wi_bWk_gap_k]          = "Bulge l gap k & j";
	state_explain[il_I_il_K_ik]                  = "Internal loop ik";
	state_explain[il_I_il_K_WiWk]                = "Internal loop jl";
	state_explain[il_I_il_gap_k_i]               = "Internal loop i gap k";
	state_explain[il_I_il_gap_k_WiWk]            = "Internal loop jl gap k";
	state_explain[il_I_il_gap_Wk_ik]             = "Internal loop ik gap l";
	state_explain[il_I_il_gap_Wk_Wi]             = "Internal loop j gap l";
	state_explain[il_I_il_gap_kWk_i]             = "Internal loop i gap k & l";
	state_explain[il_I_il_gap_kWk_Wi]            = "Internal loop j gap k & l";
	state_explain[il_pb_I_il_pb_K]               = "Internal loop -> stem ik jl";
	state_explain[il_gap_iWi_il_K_k]             = "Internal loop k gap i & j";
	state_explain[il_gap_iWi_il_K_Wk]            = "Internal loop l gap i & j";
	state_explain[il_gap_i_il_K_k]               = "Internal loop k gap i";
	state_explain[il_gap_i_il_K_WiWk]            = "Internal loop jl gap i";
	state_explain[il_gap_i_il_gap_Wk_k]          = "Internal loop k gap i & l";
	state_explain[il_gap_i_il_gap_Wk_Wi]         = "Internal loop j gap i & l";
	state_explain[il_gap_Wi_il_K_ik]             = "Internal loop ik gap j";
	state_explain[il_gap_Wi_il_K_Wk]             = "Internal loop l gap j";
	state_explain[il_gap_Wi_il_gap_k_i]          = "Internal loop i gap k & j";
	state_explain[il_gap_Wi_il_gap_k_Wk]         = "Internal loop l gap k & j";
	state_explain[mblIK]                         = "Bifurcation";
	state_explain[mbl_bWi_I_mbl_bWk_K]           = "External loop jl Bifurcation possible";
	state_explain[mbl_bWi_I_mbl_bWk_gap_Wk]      = "External loop j gap l Bifurcation possible";
	state_explain[mbl_bWi_gap_Wi_mbl_bWk_K]      = "External loop l gap j Bifurcation possible";
	state_explain[mbl_il_I_mbl_il_K_ik]          = "External loop ik";
	state_explain[mbl_il_I_mbl_il_K_WiWk]        = "External loop jl";
	state_explain[mbl_il_I_mbl_il_gap_k_i]       = "External loop i gap k";
	state_explain[mbl_il_I_mbl_il_gap_k_WiWk]    = "External loop jl gap k";
	state_explain[mbl_il_I_mbl_il_gap_kWk_i]     = "External loop i gap k & l";
	state_explain[mbl_il_I_mbl_il_gap_kWk_Wi]    = "External loop j gap k & l";
	state_explain[mbl_il_I_mbl_il_gap_Wk_ik]     = "External loop ik gap l";
	state_explain[mbl_il_I_mbl_il_gap_Wk_Wi]     = "External loop j gap l";
	state_explain[mbl_il_pb_I_mbl_il_pb_K]       = "External loop -> stem ik jl";
	state_explain[mbl_il_gap_Wi_mbl_il_K_ik]     = "External loop ik gap j";
	state_explain[mbl_il_gap_Wi_mbl_il_K_Wk]     = "External loop l gap j";
	state_explain[mbl_il_gap_Wi_mbl_il_gap_k_i]  = "External loop i gap k & j";
	state_explain[mbl_il_gap_Wi_mbl_il_gap_k_Wk] = "External loop l gap k & j";
	state_explain[mbl_il_gap_iWi_mbl_il_K_k]     = "External loop k gap i & j";
	state_explain[mbl_il_gap_iWi_mbl_il_K_Wk]    = "External loop l gap i & j";
	state_explain[mbl_il_gap_i_mbl_il_K_k]       = "External loop k gap i";
	state_explain[mbl_il_gap_i_mbl_il_K_WiWk]    = "External loop jl gap i";
	state_explain[mbl_il_gap_i_mbl_il_gap_Wk_k]  = "External loop k gap i";
	state_explain[mbl_il_gap_i_mbl_il_gap_Wk_Wi] = "External loop j gap i";

	// If the sequences are swithch, then some of the explanations must also
	// be switched
	if (flip) {
		helper::swap(state_explain[hp_init_gap_I_k], state_explain[hp_init_gap_K_i]);
		helper::swap(state_explain[hp_init_gap_Wi_ik], state_explain[hp_init_gap_Wk_ik]);
		helper::swap(state_explain[hp_init_gap_Wi_Wk], state_explain[hp_init_gap_Wk_Wi]);
		helper::swap(state_explain[hp_init_gap_iWi_k], state_explain[hp_init_gap_kWk_i]);
		helper::swap(state_explain[hp_init_gap_iWi_Wk], state_explain[hp_init_gap_kWk_Wi]);
		helper::swap(state_explain[hp_init_gap_iWk_k], state_explain[hp_init_gap_kWi_i]);
		helper::swap(state_explain[hp_init_gap_iWk_Wi], state_explain[hp_init_gap_kWi_Wk]);
		helper::swap(state_explain[stem_I_stem_gap_kWk], state_explain[stem_gap_iWi_stem_K]);
		helper::swap(state_explain[stem_no_mbl_I_stem_gap_kWk], state_explain[stem_no_mbl_gap_iWi_stem_K]);
		helper::swap(state_explain[bi_I_bk_gap_k], state_explain[bi_gap_i_bk_K]);
		helper::swap(state_explain[bi_I_bk_gap_kWk], state_explain[bi_gap_iWi_bk_K]);
		helper::swap(state_explain[bi_I_bk_gap_Wk], state_explain[bi_gap_Wi_bk_K]);
		helper::swap(state_explain[bi_gap_i_bk_gap_Wk], state_explain[bi_gap_Wi_bk_gap_k]);
		helper::swap(state_explain[bWi_I_bWk_gap_k], state_explain[bWi_gap_i_bWk_K]);
		helper::swap(state_explain[bWi_I_bWk_gap_kWk], state_explain[bWi_gap_iWi_bWk_K]);
		helper::swap(state_explain[bWi_I_bWk_gap_Wk], state_explain[bWi_gap_Wi_bWk_K]);
		helper::swap(state_explain[bWi_gap_i_bWk_gap_Wk], state_explain[bWi_gap_Wi_bWk_gap_k]);
		helper::swap(state_explain[il_I_il_gap_k_i], state_explain[il_gap_i_il_K_k]);
		helper::swap(state_explain[il_I_il_gap_k_WiWk], state_explain[il_gap_i_il_K_WiWk]);
		helper::swap(state_explain[il_I_il_gap_Wk_ik], state_explain[il_gap_Wi_il_K_ik]);
		helper::swap(state_explain[il_I_il_gap_Wk_Wi], state_explain[il_gap_Wi_il_K_Wk]);
		helper::swap(state_explain[il_I_il_gap_kWk_i], state_explain[il_gap_iWi_il_K_k]);
		helper::swap(state_explain[il_I_il_gap_kWk_Wi], state_explain[il_gap_iWi_il_K_Wk]);
		helper::swap(state_explain[il_gap_i_il_gap_Wk_k], state_explain[il_gap_Wi_il_gap_k_i]);
		helper::swap(state_explain[il_gap_i_il_gap_Wk_Wi], state_explain[il_gap_Wi_il_gap_k_Wk]);
		helper::swap(state_explain[mbl_bWi_I_mbl_bWk_gap_Wk], state_explain[mbl_bWi_gap_Wi_mbl_bWk_K]);
		helper::swap(state_explain[mbl_il_I_mbl_il_gap_k_i], state_explain[mbl_il_gap_i_mbl_il_K_k]);
		helper::swap(state_explain[mbl_il_I_mbl_il_gap_k_WiWk], state_explain[mbl_il_gap_i_mbl_il_K_WiWk]);
		helper::swap(state_explain[mbl_il_I_mbl_il_gap_kWk_i], state_explain[mbl_il_gap_iWi_mbl_il_K_k]);
		helper::swap(state_explain[mbl_il_I_mbl_il_gap_kWk_Wi], state_explain[mbl_il_gap_iWi_mbl_il_K_Wk]);
		helper::swap(state_explain[mbl_il_I_mbl_il_gap_Wk_ik], state_explain[mbl_il_gap_Wi_mbl_il_K_ik]);
		helper::swap(state_explain[mbl_il_I_mbl_il_gap_Wk_Wi], state_explain[mbl_il_gap_Wi_mbl_il_K_Wk]);
		helper::swap(state_explain[mbl_il_gap_Wi_mbl_il_gap_k_i], state_explain[mbl_il_gap_i_mbl_il_gap_Wk_k]);
		helper::swap(state_explain[mbl_il_gap_Wi_mbl_il_gap_k_Wk], state_explain[mbl_il_gap_i_mbl_il_gap_Wk_Wi]);
	}


	gcCorrect_I = true;
	gcCorrect_K = true;
	
	gcCorrect[stem_IK] = true; 
	gcCorrect[stem_I_stem_gap_kWk] = true;
	gcCorrect[stem_gap_iWi_stem_K] = true;
	gcCorrect[bWi_I_bWk_K] = true;
	gcCorrect[bWi_I_bWk_gap_k] = true;
	gcCorrect[bWi_I_bWk_gap_kWk] = true;
	gcCorrect[bWi_I_bWk_gap_Wk] = true;
	gcCorrect[bWi_gap_iWi_bWk_K] = true;
	gcCorrect[bWi_gap_i_bWk_K] = true;
	gcCorrect[bWi_gap_i_bWk_gap_Wk] = true;
	gcCorrect[bWi_gap_Wi_bWk_K] = true;
	gcCorrect[bWi_gap_Wi_bWk_gap_k] = true;

	hbi_loop[hp_init] = true;
	hbi_loop[hp_init_align_ik] = true;
	hbi_loop[hp_init_align_WiWk] = true;
	hbi_loop[hp_init_gap_I_k] = true;
	hbi_loop[hp_init_gap_I_WiWk] = true;
	hbi_loop[hp_init_gap_K_i] = true;
	hbi_loop[hp_init_gap_K_WiWk] = true;
	hbi_loop[hp_init_gap_Wi_ik] = true;
	hbi_loop[hp_init_gap_Wi_Wk] = true;
	hbi_loop[hp_init_gap_Wk_ik] = true;
	hbi_loop[hp_init_gap_Wk_Wi] = true;
	hbi_loop[hp_init_gap_iWi_k] = true;
	hbi_loop[hp_init_gap_iWi_Wk] = true;
	hbi_loop[hp_init_gap_kWk_i] = true;
	hbi_loop[hp_init_gap_kWk_Wi] = true;
	hbi_loop[hp_init_gap_iWk_k] = true;
	hbi_loop[hp_init_gap_iWk_Wi] = true;
	hbi_loop[hp_init_gap_kWi_i] = true;
	hbi_loop[hp_init_gap_kWi_Wk] = true;
	hbi_loop[bi_I_bk_K] = true;
	hbi_loop[bi_I_bk_gap_k] = true;
	hbi_loop[bi_I_bk_gap_kWk] = true;
	hbi_loop[bi_I_bk_gap_Wk] = true;
	hbi_loop[bi_gap_i_bk_K] = true;
	hbi_loop[bi_gap_i_bk_gap_Wk] = true;
	hbi_loop[bi_gap_iWi_bk_K] = true;
	hbi_loop[bi_gap_Wi_bk_K] = true;
	hbi_loop[bi_gap_Wi_bk_gap_k] = true;
	hbi_loop[bWi_I_bWk_K] = true;
	hbi_loop[bWi_I_bWk_gap_k] = true;
	hbi_loop[bWi_I_bWk_gap_kWk] = true;
	hbi_loop[bWi_I_bWk_gap_Wk] = true;
	hbi_loop[bWi_gap_i_bWk_K] = true;
	hbi_loop[bWi_gap_i_bWk_gap_Wk] = true;
	hbi_loop[bWi_gap_iWi_bWk_K] = true;
	hbi_loop[bWi_gap_Wi_bWk_K] = true;
	hbi_loop[bWi_gap_Wi_bWk_gap_k] = true;
	hbi_loop[il_I_il_K_ik] = true;
	hbi_loop[il_I_il_K_WiWk] = true;
	hbi_loop[il_I_il_gap_k_i] = true;
	hbi_loop[il_I_il_gap_k_WiWk] = true;
	hbi_loop[il_I_il_gap_Wk_ik] = true;
	hbi_loop[il_I_il_gap_Wk_Wi] = true;
	hbi_loop[il_I_il_gap_kWk_i] = true;
	hbi_loop[il_I_il_gap_kWk_Wi] = true;
	hbi_loop[il_gap_iWi_il_K_k] = true;
	hbi_loop[il_gap_iWi_il_K_Wk] = true;
	hbi_loop[il_gap_i_il_K_k] = true;
	hbi_loop[il_gap_i_il_K_WiWk] = true;
	hbi_loop[il_gap_i_il_gap_Wk_k] = true;
	hbi_loop[il_gap_i_il_gap_Wk_Wi] = true;
	hbi_loop[il_gap_Wi_il_K_ik] = true;
	hbi_loop[il_gap_Wi_il_K_Wk] = true;
	hbi_loop[il_gap_Wi_il_gap_k_i] = true;
	hbi_loop[il_gap_Wi_il_gap_k_Wk] = true;


}
#endif /*BACKTRACK*/
