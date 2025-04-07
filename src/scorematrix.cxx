#ifndef SCOREMATRIX
#define SCOREMATRIX

#include "foldalign.hxx"
#include "readfile.cxx"
#include "helper.cxx"
#include "stack_ssl.cxx"
#include "exception.cxx"
#include "matrix4d.cxx"
#include "matrix2d.cxx"
#include "prune.cxx"

#include <string>
#include <iostream>
#include <math.h>

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

//*********************************************************
// This class holds the scoremtricies needed by foldalign *
// Implemented by Jakob Hull Havgaard 2004                *
// hull@bioinf.ku.dk                                     *
//*********************************************************

class scorematrix {
public:
	//**************************************************************************
	// The class has 4 constructors.
	// a default which sets up the default values
	// one which updates the default values with those read from a file
	// the last two are the assignment and copy constructors

	// A default constructor
	scorematrix(lengthType delta, const bool global = false, const bool noPrune = false) {
		score_matrix(delta, global, noPrune);
	};

	// One that reads all or some of the matrices from a file
	scorematrix(std::string& filename, lengthType delta, const bool global = false,
		    const bool noPrune = false);

	//==========================================================================
	// This function must be called before the object can be used.
	//Expand the size of the tables.

	inline void checkSize(const int size);


	// Change the parameters which are dependent on the gc content of the
	// sequences
	void setGCparameters(const float gcSeq1, const float gcSeq2,
			     const positionType seqLength1,
			     const positionType seqLength2,
			     const bool global = false);

	void updatePruning(const positionType seqLength1, const positionType seqLength2, lengthType delta) {
	  prune_table->setNewSeqLengths(seqLength1, seqLength2, getGap(), delta);
	}

	// and an assignment constructor
	scorematrix& operator=(const scorematrix& s){
		if (this != &s) {
			clean();
			// Clean can not take care of asymTable so it is done here
			delete[] asymTable;
			copy(s);
		}
		return *this;
	};

	// The copy constructor
	scorematrix(const scorematrix& s) {copy(s);};

	// The destructor.
	~scorematrix() {
		clean();
		// asymTable must be deleted separately.
		delete[] asymTable;
	};

	//**************************************************************************
	// Functions for accessing the information

	//==========================================================================
	// The stacking bonus
	// last stack on first
	scoreType getStack(const int i, const int j, const int ii, const int jj,const int k, const int l, const int kk, const int ll) const {return (stack_matrix->get(i,j,ii,jj) + stack_matrix->get(k,l,kk,ll));};
	scoreType getStack(const int i, const int j, const int ii, const int jj) const {return stack_matrix->get(i,j,ii,jj);};

	//==========================================================================
	// The affine gap penalty
	scoreType getGap() const {return elongation_cost;};
	// The gap open penalty
	scoreType getGapOpen() const {return gap_open;};
	scoreType getBpGapOpen() const {return bp_gap_open;};
	scoreType getBpAffineGap() const {return bp_elongation_cost;};
	//==========================================================================
	// The cost of base pair substitutions
	scoreType getScore(const int i, const int j, const int k, const int l) const {return s_matrix->get(i,j,k,l);};
	// The single strand substitution cost
	scoreType getInit(const int i, const int k) const {return i_matrix->get(i,k);};

	//==========================================================================
	// The name of the matrix: default or filename
	std::string getName() const {return matrixname;};

	//==========================================================================
	// True if the bases basepairs
	bool getBasepair(const int i, const int k) const {return basepair->get(i,k);};

	//==========================================================================
	// The size of the alfabet
	int alfaSize() const {return alfa_size;};
	// The letter at position pos in the alfabet
	char getLetter(const int pos) const {return alfabet[pos];};
	// The alfabet position of the letter letter
	inline int alfa(char letter) const;

	//==========================================================================
	// Returns the multibranch loop opening score
	scoreType getMbl() const {return mbl;}
	// Returns the affine multibranch loop cost
	scoreType getMblAffine() const {return mblAffine;}
	// Returns the cost for adding a nucleotide to a mbl
	scoreType getMblNuc() const {return mblNuc;}

	//==========================================================================
	// Returns the endNonGC value
	scoreType getNonGCEnd(const int i, const int j,
	                      const int k, const int l) const {
		return endNonGC->get(i, j) + endNonGC->get(k, l);
	};

	scoreType getNonGCEnd(const int i, const int j) const {return endNonGC->get(i, j);}

	//==========================================================================
	// These functions are for handling the loop close scores

	scoreType getHpClose(const int i, const int j,
	                     const int ii, const int jj,
								const int k, const int l,
								const int kk, const int ll) const {
		return hp_close_matrix->get(i,j,ii,jj) + hp_close_matrix->get(k,l,kk,ll);
	};
	scoreType getHpClose(const int i, const int j,
	                     const int ii, const int jj) const {
		return hp_close_matrix->get(i,j,ii,jj);
	};
	scoreType getIntLoopOpen(const int bi, const int bj,
	                         const int si, const int sj) const {
		return internal_loop_matrix->get(bi,bj,si,sj);
	};
	scoreType getIntLoopClose(const int bi, const int bj,
	                          const int si, const int sj) const {
		return internal_loop_matrix->get(bi,bj,si,sj);
	};
	scoreType getHpLength(const int i, const int j) const {
		return hpLength[i]+hpLength[j];
	};
	scoreType getIntLoopLength(const int i, const int j) const;
	scoreType getBulgeLength(const int i, const int j) const {
		return bulgeLength[i]+bulgeLength[j];
	};
	scoreType getBulgeLength(const int i) const {
		return bulgeLength[i];
	};

	//==========================================================================
	// Return the pruning minimum score given a length

	inline scoreType getPruneScore(const scoreType Wi, const scoreType Wk) const {return prune_table->get(Wi, Wk);};
	inline lengthType getPruneTableLength() const {return prune_table->getLength();}

	//==========================================================================
	// Return the log (base 2)

	inline double getLog(const lengthType Wi) const {return logTable[Wi];}

private:

	//======================================================================
	// The matrix's name
	std::string matrixname; // Default or filename

	//======================================================================
	// The energy matrices and tables
	int loopTableLength; // The length of the loop cost table
	matrix4d<scoreType>* stack_matrix;  // The stacking matrix
	matrix4d<scoreType>* hp_close_matrix;  // The hair-pin open matrix
	matrix4d<scoreType>* internal_loop_matrix;  // The internal loop open/close matrix
	matrix2d<scoreType>* endNonGC;	// The extra cost for terminal non GC base-pair. Zero for GC-pairs.
	scoreType* hpLength;		 // The hairpin loop length cost
	scoreType* bulgeLength; // The bulge loop length cost
	scoreType* ilLength;    // The internal loop length cost
	scoreType ilLong;         // The scale factor for long internal loops
	scoreType bulgeLong;      // The scale factor for long bulge loops
	scoreType hpLong;         // The scale factor for long hairpin loops
	scoreType asym;				// The per base cost for asymetric internal loops
	scoreType masym;				// The per base cost for asymetric internal loops
	scoreType* asymTable;	   // A table of asymetry cost.
	lengthType asymTableSize;   // The size of the asymTable is 2*asymTalbeSize+1
	scoreType mbl;				// Multibranchloop opening cost
	scoreType mblAffine;	   // The affine multibranch cost
	scoreType mblNuc;         // The cost of adding a nucleotide to the mbl


	//======================================================================
	// The substitution and base pair matrices
	matrix4d< scoreType >* s_matrix;  // The substitution matrix
	matrix2d< scoreType >* i_matrix;    // The scorematrix for initial alignment
	matrix2d< bool >* basepair;   // The basepairing matrix


	//======================================================================
	// The alfabet
	char* alfabet;     // The alfabet
	int alfa_size;     // The size of the alfabet
	int tmp_alfa_size; // A tmp variable alfa_size. Used when the size of the matrices is changed
	bool alfadie;      // Used to make sure the matricies in a file is ordered correctly

	//======================================================================
	// Gaps
	scoreType gap_open;      // The gap open score
	scoreType elongation;    // The gap elongation score (the one read from file)
	scoreType bp_gap_open;      // The gap open score for basepairs
	scoreType bp_elongation;    // The gap elongation score for basepairs

	scoreType elongation_cost; // Holds the full elongation cost
	scoreType bp_elongation_cost; // Holds the full elongation cost

	//======================================================================
	// Prunning
	prune* prune_table;	// Alignments with scores below these scores are pruned

	//======================================================================
	// The log table
	double* logTable; 		// Table of logrithms. Base 2.

	//======================================================================
	// Define the defualts
	inline void score_matrix(lengthType delta, const bool global, const bool noPrune);
	inline void setNameGaps(std::string name,
					scoreType gapOpen, scoreType gapElongation,
					scoreType bpGapOpen, scoreType bpGapElongation);

	inline void setImatrix(
			const scoreType aa,
			const scoreType ac,
			const scoreType ag,
			const scoreType au,
			const scoreType cc,
			const scoreType cg,
			const scoreType cu,
			const scoreType gg,
			const scoreType gu,
			const scoreType uu,
			const scoreType gap);

	inline void setS_matrix(
			const scoreType auau,
			const scoreType aucg,
			const scoreType augc,
			const scoreType augu,
			const scoreType auua,
			const scoreType auug,
			const scoreType cgcg,
			const scoreType cggc,
			const scoreType cggu,
			const scoreType cgua,
			const scoreType cgug,
			const scoreType gcgc,
			const scoreType gcgu,
			const scoreType gcua,
			const scoreType gcug,
			const scoreType gugu,
			const scoreType guua,
			const scoreType guug,
			const scoreType uaua,
			const scoreType uaug,
			const scoreType ugug);

	inline void submatrixScale0_33();
	inline void submatrixScale0_5();
	inline void submatrixScale2();
	inline void submatrixScale3();

	//======================================================================
	// Reading functions
	inline void setAlfabet(readfile*& file);     // Read a new alfabet from file
	inline void setBasepair(readfile*& file);    // Read a new basepairings matrix from file
	inline void setInitMatrix(readfile*& file); // Read a new init matrix from file
	inline void setLoopTable(readfile*&file);    // Read a new loop cost table
	inline void setMisc(readfile*& file);        // Read the misc values

	//======================================================================
	// Helper functions
	inline void store(std::string name, scoreType value);        // A helper function to setMisc
	inline void parseLine(std::string line, int*& i_letters, int& len); // Parse a annotation line
	template<class data>
	inline void read4dMatrixFromFile(const int& size, readfile*& file, matrix4d<data>*& matrix);
	inline int getLetterIndex(std::string line, int& prev);

	//======================================================================
	//Misc calculations
	inline scoreType calcLongScore(const int factor, const lengthType i, const lengthType length, const scoreType zero);
	inline void calcScore(const int factor, const int size, scoreType*& array);
	inline void calcAsym(const int size);
	inline void makeLogTable(lengthType end);

	//======================================================================
	// Helper functions for coping and destruction
	// Deletes every thing except the asymetry table
	inline void clean();
	// Copies everything.
	inline void copy(const scorematrix& s);


	// Used to store table values read from a file
	struct prunes {
		scoreType score;
		lengthType index;
	};
	inline void assignFromStack(scoreType*& store, stack_ssl<prunes, lengthType>& lineStack);
};

//********************************************
// Read a scorematrix file

inline scorematrix::scorematrix(std::string& filename, lengthType delta, const bool global, const bool noPrune) {

	score_matrix(delta, global, noPrune); // Sets all parameters to the default values
	matrixname = filename;
	std::string line;

        // Opening the new file
	readfile* file;
	try {
		file = new readfile(filename);
	}
	catch ( ... ) {
		std::string error = "It is not possible to open the scorematrix file: ";
		error += filename;
		throw exception(error, false);
	}

	while (file->get_line(line)) {

		if (!line.compare("Alfabet:")) {setAlfabet(file);}
		else if (!line.compare("Base-pair:")) {setBasepair(file);}
		else if (!line.compare("Base-pair substitution:")) {
			read4dMatrixFromFile(alfaSize(), file, s_matrix);
		}
		else if (!line.compare("Single strand substitution:")) {setInitMatrix(file);}
		else if (!line.compare("Stacking:")) {
			read4dMatrixFromFile(alfaSize(), file, stack_matrix);
		}
		else if (!line.compare("Hairpin Close:")) {
			read4dMatrixFromFile(alfaSize(), file, hp_close_matrix);
		}
		else if (!line.compare("Internal loop:")) {
			read4dMatrixFromFile(alfaSize(), file, internal_loop_matrix);
		}
		else if (!line.compare("5' Dangle:")) {
			std::cerr << "Note! The 5' dangle matrix is no longer used. Any values will be ignored." << std::endl;
			std::cerr << "Skipping to the next empty line" << std::endl;
			file->skip_to_empty();
		}
		else if (!line.compare("3' Dangle:")) {
			std::cerr << "Note! The 3' dangle matrix is no longer used. Any values will be ignored." << std::endl;
			std::cerr << "Skipping to the next empty line" << std::endl;
			file->skip_to_empty();
		}
		else if (!line.compare("Loop length costs:")) {setLoopTable(file);}
		else if (!line.compare("Miscellaneous:")) {setMisc(file);}
		else if (!line.compare("Pruning:")) {
			prune_table->buildLocalPruneTableFromFile(file);
		}
		else if (!line.compare("")) {}
		else {
			std::string error = "The score matrix do not have the right format\n";
			error += line;
			throw exception(error, false);
		}
	}

	delete file;

	alfa_size = tmp_alfa_size;

	// Setup the gap opening cost
	for(int i=0; i<alfa_size; i++) {
		i_matrix->set(i,0, gap_open);
		i_matrix->set(0,i, gap_open);
	}
	i_matrix->set(0,0, gap_open);

	calcAsym(loopTableLength);
	prune_table->reSizeTable(loopTableLength);

	elongation_cost = elongation - gap_open;
	bp_elongation_cost = bp_elongation - bp_gap_open;
}

inline void scorematrix::setGCparameters(const float gcSeq1, const float gcSeq2,
				const positionType seqLength1,
				const positionType seqLength2,
				const bool global) {


//	const float gcHigh = gcSeq1 > gcSeq2 ? gcSeq1 : gcSeq2;
//	const float gcLow = gcSeq1 <= gcSeq2 ? gcSeq1 : gcSeq2;

	if (global) {

return;
	}
	else {
return;
/*	  delete prune_table;
		if (gcHigh < 0.3) {
			if (gcLow < 0.2) {
				std::string name = matrixname + "_gc_0.2_0.1";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-166, 1, loopTableLength, seqLength1, seqLength2, getGap(), delta);
				submatrixScale0_33();
			}
			else {
				std::string name = matrixname + "_gc_0.2_0.2";
				setNameGaps(name, -36, -18, -73, -36);
				prune_table->updatePruneTable(-227, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_33();
			}
		}
		else if (gcHigh < 0.4) {
			if (gcLow < 0.1) {
				std::string name = matrixname + "_gc_0.3_0.0";
				setNameGaps(name, -55, -28, -110, -55);
				prune_table->updatePruneTable(-247, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_5();
			}
			else if (gcLow < 0.2) {
				std::string name = matrixname + "_gc_0.3_0.1";
				setNameGaps(name, -55, -28, -110, -55);
				prune_table->updatePruneTable(-178, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_33();
			}
			else if (gcLow < 0.3) {
				std::string name = matrixname + "_gc_0.3_0.2";
				setNameGaps(name, -55, -28, -110, -55);
				prune_table->updatePruneTable(-236, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_5();
			}
			else {
				std::string name = matrixname + "_gc_0.3_0.3";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-341, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
		}
		else if (gcHigh < 0.5) {
			if (gcLow < 0.1) {
				std::string name = matrixname + "_gc_0.4_0.0";
				setNameGaps(name, -110, -55, -220, -110);
				prune_table->updatePruneTable(-247, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_5();
			}
			else if (gcLow < 0.2) {
				std::string name = matrixname + "_gc_0.4_0.1";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-253, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_5();
			}
			else if (gcLow < 0.3) {
				std::string name = matrixname + "_gc_0.4_0.2";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-200, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
			else if (gcLow < 0.4) {
				std::string name = matrixname + "_gc_0.4_0.3";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-272, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
			else {
				std::string name = matrixname + "_gc_0.4_0.4";
				setNameGaps(name, -110, -55, -220, -110);
				prune_table->updatePruneTable(-611, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
		}
		else if (gcHigh < 0.6) {
			if (gcLow < 0.1) {
				std::string name = matrixname + "_gc_0.5_0.0";
				setNameGaps(name, -55, -28, -110, -55);
				prune_table->updatePruneTable(-191, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_33();
			}
			else if (gcLow < 0.2) {
				std::string name = matrixname + "_gc_0.5_0.1";
				setNameGaps(name, -110, -55, -220, -110);
				prune_table->updatePruneTable(-299, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_5();
			}
			else if (gcLow < 0.3) {
				std::string name = matrixname + "_gc_0.5_0.2";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-233, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
			else if (gcLow < 0.4) {
				std::string name = matrixname + "_gc_0.5_0.3";
				setNameGaps(name, -330, -165, -660, -330);
				prune_table->updatePruneTable(-612, 1, loopTableLength, seqLength1, seqLength2, getGap());
//				submatrixScale0_5();
			}
			else if (gcLow < 0.5) {
				std::string name = matrixname + "_gc_0.5_0.4";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-510, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_5();
			}
			else {
				std::string name = matrixname + "_gc_0.5_0.5";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-1050, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale2();
			}
		}
		else if (gcHigh < 0.7) {
			if (gcLow < 0.1) {
				std::string name = matrixname + "_gc_0.6_0.0";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-132, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
			else if (gcLow < 0.2) {
				std::string name = matrixname + "_gc_0.6_0.1";
				setNameGaps(name, -110, -55, -220, -110);
				prune_table->updatePruneTable(-280, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_5();
			}
			else if (gcLow < 0.3) {
				std::string name = matrixname + "_gc_0.6_0.2";
				setNameGaps(name, -330, -165, -660, -330);
				prune_table->updatePruneTable(-418, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale2();
			}
			else if (gcLow < 0.4) {
				std::string name = matrixname + "_gc_0.6_0.3";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-201, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
			else if (gcLow < 0.5) {
				std::string name = matrixname + "_gc_0.6_0.4";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-811, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale3();
			}
			else if (gcLow < 0.6) {
				std::string name = matrixname + "_gc_0.6_0.5";
				setNameGaps(name, -330, -165, -660, -330);
				prune_table->updatePruneTable(-334, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
			else {
				std::string name = matrixname + "_gc_0.6_0.6";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-286, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
		}
		else {
			if (gcLow < 0.3) {
				std::string name = matrixname + "_gc_0.7_0.2";
				setNameGaps(name, -330, -165, -660, -330);
				prune_table->updatePruneTable(-224, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale2();
			}
			else if (gcLow < 0.4) {
				std::string name = matrixname + "_gc_0.7_0.3";
				setNameGaps(name, -330, -165, -660, -330);
				prune_table->updatePruneTable(-205, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
			else if (gcLow < 0.5) {
				std::string name = matrixname + "_gc_0.7_0.4";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-251, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale0_5();
			}
			else if (gcLow < 0.6) {
				std::string name = matrixname + "_gc_0.7_0.5";
				setNameGaps(name, -220, -110, -440, -220);
				prune_table->updatePruneTable(-147, 1, loopTableLength, seqLength1, seqLength2, getGap());
			}
			else {
				std::string name = matrixname + "_gc_0.7_0.6";
				setNameGaps(name, -330, -165, -660, -330);
				prune_table->updatePruneTable(-202, 1, loopTableLength, seqLength1, seqLength2, getGap());
				submatrixScale3();
			}
		}*/
	}
}

inline void scorematrix::submatrixScale0_33() {
	setImatrix(6, -7, -6, -6, 4, -8, -5, 3, -7, 4, gap_open);
	setS_matrix(14,  4,   6,  0, 3, -5,
	            17,  7,  -1,  6, 1,
	            18,  2,   3, -2,
	            11, -2,  -6,
	            16,  1,
	            11);
}

inline void scorematrix::submatrixScale0_5() {
	setImatrix(10, -11, -9, -10, 6, -13, -8, 5, -10, 7, gap_open);
	setS_matrix(21,  6,  10,  0, 4, -7,
	            27, 10,  -2, 10, 2,
	            28,  4,   5, -3,
	            17, -3, -10,
	            24,  1,
	            16);
}

inline void scorematrix::submatrixScale2() {
	setImatrix(38, -44, -36, -38, 22, -50, -30, 18, -40, 26, gap_open);
	setS_matrix(84, 22, 38, 0, 16, -28,
		        106, 40, -8, 38, 8,
		        112, 14, 18, -12,
		        68, -12, -38,
		        96, 4,
		        64);
}

inline void scorematrix::submatrixScale3() {
	setImatrix(57, -66, -54, -57, 33, -75, -45, 27, -60, 39, gap_open);
	setS_matrix(126,  33,  57,   0, 24, -42,
		        159,  60, -12,  57, 12,
		        168,  21,  27, -18,
		        102, -18, -57,
		        144,   6,
		         96);
}

//*********************************************
// Assignment operator for the scorematrix object

inline int scorematrix::alfa(char letter) const {
	char upper = toupper(letter);
	if(upper=='T') {upper = 'U';}
	for(int i =0; i<alfa_size; i++) {
		if (alfabet[i]==upper) {return i;}
	}
	std::cerr << "Unknown nucleotide: " << letter << " found." << std::endl;
	return alfa_size -1;
}

inline scoreType scorematrix::getIntLoopLength(const int i, const int j) const {

	scoreType res = ilLength[i+j];

	// The Ninno equation (in the case where the cost is independent of the length)
	// A value of 5 is ten times the values set in the miscloop file of dynalign
	res += asymTable[i-j+asymTableSize];

	return res;
}

//******************************************************
// The private functions

inline void scorematrix::clean() {
	// Note assymTable is not deleted by this function.

	delete s_matrix;
	delete stack_matrix;
	delete hp_close_matrix;
	delete internal_loop_matrix;
	delete i_matrix;
	delete basepair;
	delete endNonGC;
	delete prune_table;
	delete[] alfabet;
	delete[] hpLength;
	delete[] bulgeLength;
	delete[] ilLength;
	delete[] logTable;
}

inline void scorematrix::copy(const scorematrix& s) {
	matrixname = s.matrixname;

	alfa_size = s.alfaSize();
	tmp_alfa_size = s.tmp_alfa_size;
	alfabet = new char[alfa_size];
	helper::copyArray(alfabet, s.alfabet, alfa_size);
	alfadie = s.alfadie;

	gap_open = s.gap_open;
	elongation = s.elongation;
	elongation_cost = s.elongation_cost;
	bp_gap_open = s.bp_gap_open;
	bp_elongation = s.bp_elongation;
	bp_elongation_cost = s.bp_elongation_cost;
	asym = s.asym;
	masym = s.masym;
	mbl = s.mbl;
	mblAffine = s.mblAffine;
	mblNuc = s.mblNuc;

	s_matrix = new matrix4d<scoreType>(*s.s_matrix);
	stack_matrix = new matrix4d<scoreType>(*s.stack_matrix);
	hp_close_matrix = new matrix4d<scoreType>(*s.hp_close_matrix);
	internal_loop_matrix = new matrix4d<scoreType>(*s.internal_loop_matrix);

	i_matrix = new matrix2d<scoreType>(*s.i_matrix);
	basepair = new matrix2d<bool>(*s.basepair);
	endNonGC = new matrix2d<scoreType>(*s.endNonGC);

	prune_table = new prune(*s.prune_table);

	loopTableLength = s.loopTableLength;

	hpLong = s.hpLong;
	hpLength = new scoreType[loopTableLength];
	helper::copyArray(hpLength, s.hpLength, loopTableLength);

	bulgeLong = s.bulgeLong;
	bulgeLength = new scoreType[loopTableLength];
	helper::copyArray(bulgeLength, s.bulgeLength, loopTableLength);

	ilLong = s.ilLong;
	ilLength = new scoreType[loopTableLength];
	helper::copyArray(ilLength, s.ilLength, loopTableLength);

	// A dummy pointer. It will be deleted by the calcAsym function
	asymTableSize = s.asymTableSize;
	asymTable = new scoreType[2];
	calcAsym(asymTableSize);

	asym = s.asym;
	masym = s.masym;
	mbl = s.mbl;
	mblAffine = s.mblAffine;
	mblNuc = s.mblNuc;

	logTable = new double[loopTableLength];
	helper::copyArray(logTable, s.logTable, loopTableLength);

}

inline scoreType scorematrix::calcLongScore(const int factor, const lengthType i, const lengthType length, const scoreType zero) {
	return static_cast<scoreType>(log((static_cast<float>(factor*i)/static_cast<float>(factor*(length))))) + zero;
}

// A simple way of calculate the cost of very long loops.
// A better function would be based on the log of i scaled with the factor and the top value
inline void scorematrix::calcScore(const int factor, const int size, scoreType*& array) {
	scoreType* tmp = new scoreType[size+1];
	helper::copyArray(tmp, array, loopTableLength);
	if (factor != 0) {
		for(int i=loopTableLength; i<=size; i++) {
			tmp[i] = calcLongScore(factor, i, loopTableLength-1, array[loopTableLength-1]);
		}
	}
	else {
		for(int i=loopTableLength; i<=size; i++) {
			tmp[i] = 0;
		}
	}
	delete[] array;
	array = tmp;
}

inline void scorematrix::calcAsym(const int size) {

	asymTableSize = size;
	const int real_size = 2*size +1;
	delete[] asymTable;

	asymTable = new scoreType[real_size];
	for (int i = 0; i < real_size; i++) {
		asymTable[i] = asym*(i - size);

		// The final cost has to be negative
		if (asymTable[i] > 0) {asymTable[i] *= -1;}

		// Limits the maximum asymmetric cost
		if (asymTable[i] < masym) {asymTable[i] = masym;}

	}
}

inline void scorematrix::makeLogTable(lengthType end) {
	// Make sure the zero index is not 0 or very small
	if (end > 0) {logTable[0] = 1;}

	// The length of a sequence with window size Wi is Wi +1 hence the i+1
	for(lengthType i = 1; i < end; i++) {
		logTable[i] = log(i+1)/log(2);
	}
}

inline void scorematrix::checkSize(const int size) {
	if (size >= loopTableLength) {
		calcScore(hpLong, size, hpLength);
		calcScore(bulgeLong, size, bulgeLength);
		calcScore(ilLong, size, ilLength);
		calcAsym(size);
		helper::expandArray(logTable, loopTableLength, size+1);
		makeLogTable(size);
		loopTableLength = size+1;
	}
	if (size >= prune_table->getLength() -1) {
		prune_table->reSizeTable(size);
	}
}

//********************
// Functions for reading matrixes

inline void scorematrix::setAlfabet(readfile*& file) {
	if (alfadie) {
		std::string error = "The alfabet has to be defined before anything else in the score matrix file.";
		throw exception(error, false);
	}
	delete[] alfabet;
	int prev = 0;
	int pos, end_pos;
	std::string line;
	file->get_line(line);
	pos = line.find_last_not_of(" ");
	end_pos = line.find_last_of(" ", pos);

	tmp_alfa_size = atoi(line.substr(pos, (end_pos - pos)).c_str());
	alfabet = new char[tmp_alfa_size];
	for(int i=0; i<tmp_alfa_size; i++) {
		pos = line.find_first_not_of(" ", prev);
		alfabet[i] = line[pos];
		prev=pos+1;
	}
}

inline void scorematrix::setBasepair(readfile*& file) {
	alfadie=true;
	int as = alfaSize();
	delete basepair;
	as = tmp_alfa_size;
	basepair = new matrix2d<bool>(as);
	basepair->init(false);

	std::string line;
	int len = as-1;
	int prev = 0;
	int pos, end_pos;
	int* i_letters = new int[len];
	int k_letter;

	file->get_line(line);
	parseLine(line, i_letters, len);
	for(int i=0; i<len; i++) {
		file->get_line(line);

		int length = line.length();
		k_letter = getLetterIndex(line, length);

		prev = 0;
		for(int j=0; j<len; j++) {
			pos = line.find_first_not_of(" ", prev);
			end_pos = line.find(" ", pos);
			if (!line.substr(pos, (end_pos - pos)).compare("1")) {
				basepair->set(i_letters[j], k_letter, true);
			}
			else {
				basepair->set(i_letters[j], k_letter, false);
			}
			prev = end_pos+1;
		}
	}
	delete[] i_letters;
}

inline void scorematrix::setInitMatrix(readfile*& file) {
	alfadie=true;

	delete i_matrix;
	i_matrix = new matrix2d<scoreType>(alfaSize());
	i_matrix->init(scoreType(0));

	std::string line;
	int len = alfaSize() -2; // -1 for gaps and -1 for ambiguous nucleotides
	int prev = 0;
	int* i_letters = new int[len];
	int k_letter;

	file->get_line(line);
	parseLine(line, i_letters, len);
	for(int i=0; i<len; i++) {
		file->get_line(line);

		int length = line.length();
		k_letter = getLetterIndex(line, length);

		prev = 0;
		for(int j=0; j<len; j++) {
			i_matrix->set(i_letters[j], k_letter, helper::getValue(prev, line));
		}
	}
	delete[] i_letters;
}

template< class data>
inline void scorematrix::read4dMatrixFromFile(const int& size, readfile*& file, matrix4d<data>*& matrix) {

	alfadie=true;

	delete matrix;
	matrix = new matrix4d<data>(size);
	matrix->init(0);
	std::string line;
	int len = (size-2)*(size-2); // -1 for gaps and -1 for ambiguous nucleotides

	int prev = 0;
	int* i_letters = new int[len];
	int* j_letters = new int[len];
	int k_letter;
	int l_letter;

	file->get_line(line);
	parseLine(line, i_letters, len);

	file->get_line(line);
	parseLine(line, j_letters, len);
	for(int i=0; i<len; i++) {
		file->get_line(line);
		prev=line.length();
		l_letter = getLetterIndex(line, prev);
		k_letter = getLetterIndex(line, prev);

		prev = 0;
		for(int j=0; j<len; j++) {
			matrix->set(i_letters[j], j_letters[j], k_letter, l_letter, helper::getValue(prev, line));
		}
	}
	delete[] i_letters;
	delete[] j_letters;
}

inline int scorematrix::getLetterIndex(std::string line, int& prev) {
		int pos = line.find_last_not_of(" ", prev);
		char letter = line[pos];
		prev = pos -1;
		return alfa(letter);
}


inline void scorematrix::setLoopTable(readfile*& file) {

	alfadie=true; // It is no longer allowed to change the alfabet.
	std::string line;

	stack_ssl<prunes, lengthType> hpStack = stack_ssl<prunes, lengthType>();
	stack_ssl<prunes, lengthType> blStack = stack_ssl<prunes, lengthType>();
	stack_ssl<prunes, lengthType> ilStack = stack_ssl<prunes, lengthType>();

	// Add the zero values
	prunes zero = {0, 0};
	hpStack.push(zero);
	blStack.push(zero);
	ilStack.push(zero);

	lengthType last_index = 0;

	while(file->get_line_failEmpty(line)) {
		// Read the values stored in the file
		int prev = 0;
		lengthType index = helper::getValue(prev, line);
		if (index != last_index+1) {std::cerr << "Warning: Length " << index << " in the Loop length costs table is shorter than the privious length " << last_index << ". Parts of the table may be corrupt or missing. Ignoring this value." << std::endl; continue;}

		prunes hp = {scoreType(helper::getValue(prev, line)), index};
		prunes bl = {scoreType(helper::getValue(prev, line)), index};
		prunes il = {scoreType(helper::getValue(prev, line)), index};

		hpStack.push(hp);
		blStack.push(bl);
		ilStack.push(il);

		last_index = index;
	}

	// Shut down the old arrays an make new ones.
	delete[] hpLength;
	delete[] bulgeLength;
	delete[] ilLength;
	loopTableLength = last_index+1;
	hpLength = new scoreType[loopTableLength];
	bulgeLength = new scoreType[loopTableLength];
	ilLength = new scoreType[loopTableLength];

	assignFromStack(hpLength, hpStack);
	assignFromStack(bulgeLength, blStack);
	assignFromStack(ilLength, ilStack);

}


inline void scorematrix::assignFromStack(scoreType*& store, stack_ssl<prunes, lengthType>& lineStack) {

	const lengthType size = lineStack.size();
	for(lengthType i = 0; i < size; i++) {
		prunes prune = lineStack.pop();
		store[prune.index] = prune.score;
	}
//	store[0] = store[1];
}

inline void scorematrix::setMisc(readfile*& file) {
	std::string line;
	std::string name;
	scoreType value;
	while (file->get_line_failEmpty(line)) {
		int prev = 0;
		name = helper::findName(prev, line);
		value = helper::getValue(prev,line);
		store(name, value);
	}
}

inline void scorematrix::store(std::string name, scoreType value) {
	     if (!name.compare("Gap_open:")) {gap_open=value;}
	else if (!name.compare("Elongation_bonus:")) {elongation=value;}
	else if (!name.compare("Stem_gap_open:")) {bp_gap_open=value;}
	else if (!name.compare("Stem_gap_elongation_bonus:")) {bp_elongation=value;}
	else if (!name.compare("Multibranchloop:")) {mbl=value;}
	else if (!name.compare("Multibranchloop_helix:")) {mblAffine=value;}
	else if (!name.compare("Multibranchloop_nucleotide:")) {mblNuc=value;}
	else if (!name.compare("Multibranchloop_non_GC_stem_end:")) {
		endNonGC->init(value);
		endNonGC->set(2,3,0);
		endNonGC->set(3,2,0);
	}
	else if (!name.compare("Asymmetric_cost:")) {asym=value;}
	else if (!name.compare("Asymmetric_cost_limit:")) {masym=value;}
	else if (!name.compare("Long_hairpin_loop_factor:")) {ilLong=value;}
	else if (!name.compare("Long_bulge_loop_factor:")) {bulgeLong=value;}
	else if (!name.compare("Long_Internal_loop_factor:")) {hpLong=value;}
	else if (!name.compare("Linear_prunings_coefficient:")) {prune_table->setPruneCoefficient(value);}
	else {std::cerr << "Warning: Unknown parameter " << name << " with value " << value << " found" << std::endl;}
}

//*****************************************************
// These function reads and parses an input line

inline void scorematrix::parseLine(std::string line, int*& i_letters, int& len) {
	char letter;
	int pos=0;
	int length = line.length();
	for(int i=0; i<len; i++) {
		while ((line[pos] == ' ') || (line[pos] == '\t')) {pos++;}
		if (pos > length) {
			std::string error = "Could not read score matrix. Somethings wrong with this line\n";
			error += line;
			throw exception(error, false);
		}

		letter = line[pos];
		i_letters[i] = alfa(letter);
		pos++;
	}
}

inline void scorematrix::score_matrix(lengthType delta, const bool global, const bool noPrune) {
	//*****************************************************
	// This function sets the default values.

	// The alphabet can stil be changed.
	alfadie=false;
	alfa_size = tmp_alfa_size = 6;
	alfabet = new char[alfa_size];
	alfabet[0] = '-';
	alfabet[1] = 'A';
	alfabet[2] = 'C';
	alfabet[3] = 'G';
	alfabet[4] = 'U';
	alfabet[5] = 'N';

	// Make the matrices. They will be filled below
	loopTableLength = 31;
	int i_size=alfaSize();
	s_matrix = new matrix4d<scoreType>(i_size);
	stack_matrix = new matrix4d<scoreType>(i_size);
	hp_close_matrix = new matrix4d<scoreType>(i_size);
	internal_loop_matrix = new matrix4d<scoreType>(i_size);
	i_matrix = new matrix2d<scoreType>(i_size);
	basepair = new matrix2d<bool>(i_size);
	endNonGC = new matrix2d<scoreType>(i_size);

	hpLength = new scoreType[loopTableLength];
	bulgeLength = new scoreType[loopTableLength];
	ilLength = new scoreType[loopTableLength];

	// Build log table
	logTable = new double[loopTableLength];
	makeLogTable(loopTableLength);

	if (global) {

		setNameGaps("default_global", -50, -25, -100, -50);
		prune_table = new prune(-200, 1, loopTableLength, 0, 0, delta, elongation, global, noPrune);

		// The gaps will be set after the i_matrix has been initilized
		s_matrix->init(scoreType(0));
		setS_matrix(21,  5,  9,  0,  4, -6,
		            26, 10, -2, 10,  2,
					28,  4,  4, -2,
					17, -2, -9,
					24,  1,
					16);

		i_matrix->init(scoreType(0));
		setImatrix(10, -10, -8, -9, 6, -12, -7, 5, -10, 6, gap_open);
		//i_matrix->print_matrix();
	}
	else {
		setNameGaps("default_local", -110, -55, -220, -110);
		prune_table = new prune(-400, 1, loopTableLength, 0, 0, delta, elongation, global, noPrune);

		// The gaps will be set after the i_matrix has been initilized
		s_matrix->init(scoreType(0));
		setS_matrix(42, 11, 19, 0, 8, -14,
		            53, 20, -4, 19, 4,
		            56, 7, 9, -6,
		            34, -6, -19,
		            48, 2,
		            32);
		i_matrix->init(scoreType(0));
		setImatrix(19, -22, -18, -19, 11, -25, -15, 9, -20, 13, gap_open);
		// Setting the gaps
//		setGap4(s_matrix);
	}

	// The non GC ends set. Ambigouos treated as non GC
	scoreType default_nonGC_score = -5;
	endNonGC->init(default_nonGC_score);
	endNonGC->set(2,3,0);
	endNonGC->set(3,2,0);

	// Setting the misc values
	ilLong = -11;
	bulgeLong = -11;
	hpLong = -11;
	asym = -5;
	masym = -30;
	mbl=0;
	mblAffine=-4;
	mblNuc=-1;

	// A dummy memory cell. It will be deleted and resized calcAsym
	asymTable = new scoreType[2];
	calcAsym(loopTableLength);

	// Fill the matrices
	basepair->init(false);
	basepair->set(1,4,true);
	basepair->set(2,3,true);
	basepair->set(3,2,true);
	basepair->set(3,4,true);
	basepair->set(4,1,true);
	basepair->set(4,3,true);

	stack_matrix->init(scoreType(0));
	stack_matrix->set(1,4,1,4,9);
	stack_matrix->set(2,3,1,4,21);
	stack_matrix->set(3,2,1,4,24);
	stack_matrix->set(3,4,1,4,13);
	stack_matrix->set(4,1,1,4,13);
	stack_matrix->set(4,3,1,4,10);
	stack_matrix->set(1,4,2,3,22);
	stack_matrix->set(2,3,2,3,33);
	stack_matrix->set(3,2,2,3,34);
	stack_matrix->set(3,4,2,3,25);
	stack_matrix->set(4,1,2,3,24);
	stack_matrix->set(4,3,2,3,15);
	stack_matrix->set(1,4,3,2,21);
	stack_matrix->set(2,3,3,2,24);
	stack_matrix->set(3,2,3,2,33);
	stack_matrix->set(3,4,3,2,21);
	stack_matrix->set(4,1,3,2,21);
	stack_matrix->set(4,3,3,2,14);
	stack_matrix->set(1,4,3,4,6);
	stack_matrix->set(2,3,3,4,14);
	stack_matrix->set(3,2,3,4,15);
	stack_matrix->set(3,4,3,4,5);
	stack_matrix->set(4,1,3,4,10);
	stack_matrix->set(4,3,3,4,-3);
	stack_matrix->set(1,4,4,1,11);
	stack_matrix->set(2,3,4,1,21);
	stack_matrix->set(3,2,4,1,22);
	stack_matrix->set(3,4,4,1,14);
	stack_matrix->set(4,1,4,1,9);
	stack_matrix->set(4,3,4,1,6);
	stack_matrix->set(1,4,4,3,14);
	stack_matrix->set(2,3,4,3,21);
	stack_matrix->set(3,2,4,3,25);
	stack_matrix->set(3,4,4,3,-13);
	stack_matrix->set(4,1,4,3,13);
	stack_matrix->set(4,3,4,3,5);

	hp_close_matrix->init(scoreType(0));
	hp_close_matrix->set(1,4,1,1,3);
	hp_close_matrix->set(2,3,1,1,15);
	hp_close_matrix->set(3,2,1,1,11);
	hp_close_matrix->set(3,4,1,1,-2);
	hp_close_matrix->set(4,1,1,1,5);
	hp_close_matrix->set(4,3,1,1,5);
	hp_close_matrix->set(1,4,1,2,5);
	hp_close_matrix->set(2,3,1,2,15);
	hp_close_matrix->set(3,2,1,2,15);
	hp_close_matrix->set(3,4,1,2,5);
	hp_close_matrix->set(4,1,1,2,3);
	hp_close_matrix->set(4,3,1,2,3);
	hp_close_matrix->set(1,4,1,3,3);
	hp_close_matrix->set(2,3,1,3,14);
	hp_close_matrix->set(3,2,1,3,13);
	hp_close_matrix->set(3,4,1,3,3);
	hp_close_matrix->set(4,1,1,3,6);
	hp_close_matrix->set(4,3,1,3,6);
	hp_close_matrix->set(1,4,1,4,3);
	hp_close_matrix->set(2,3,1,4,18);
	hp_close_matrix->set(3,2,1,4,21);
	hp_close_matrix->set(3,4,1,4,3);
	hp_close_matrix->set(4,1,1,4,5);
	hp_close_matrix->set(4,3,1,4,5);
	hp_close_matrix->set(1,4,2,1,1);
	hp_close_matrix->set(2,3,2,1,10);
	hp_close_matrix->set(3,2,2,1,11);
	hp_close_matrix->set(3,4,2,1,1);
	hp_close_matrix->set(4,1,2,1,2);
	hp_close_matrix->set(4,3,2,1,2);
	hp_close_matrix->set(1,4,2,2,2);
	hp_close_matrix->set(2,3,2,2,9);
	hp_close_matrix->set(3,2,2,2,7);
	hp_close_matrix->set(3,4,2,2,2);
	hp_close_matrix->set(4,1,2,2,1);
	hp_close_matrix->set(4,3,2,2,1);
	hp_close_matrix->set(1,4,2,3,15);
	hp_close_matrix->set(2,3,2,3,29);
	hp_close_matrix->set(3,2,2,3,24);
	hp_close_matrix->set(3,4,2,3,15);
	hp_close_matrix->set(4,1,2,3,12);
	hp_close_matrix->set(4,3,2,3,17);
	hp_close_matrix->set(1,4,2,4,2);
	hp_close_matrix->set(2,3,2,4,8);
	hp_close_matrix->set(3,2,2,4,5);
	hp_close_matrix->set(3,4,2,4,2);
	hp_close_matrix->set(1,4,3,1,11);
	hp_close_matrix->set(2,3,3,1,22);
	hp_close_matrix->set(3,2,3,1,24);
	hp_close_matrix->set(3,4,3,1,9);
	hp_close_matrix->set(4,1,3,1,14);
	hp_close_matrix->set(4,3,3,1,8);
	hp_close_matrix->set(1,4,3,2,12);
	hp_close_matrix->set(2,3,3,2,20);
	hp_close_matrix->set(3,2,3,2,29);
	hp_close_matrix->set(3,4,3,2,11);
	hp_close_matrix->set(4,1,3,2,12);
	hp_close_matrix->set(4,3,3,2,12);
	hp_close_matrix->set(1,4,3,3,2);
	hp_close_matrix->set(2,3,3,3,16);
	hp_close_matrix->set(3,2,3,3,14);
	hp_close_matrix->set(3,4,3,3,3);
	hp_close_matrix->set(4,1,3,3,7);
	hp_close_matrix->set(4,3,3,3,3);
	hp_close_matrix->set(1,4,3,4,-2);
	hp_close_matrix->set(2,3,3,4,11);
	hp_close_matrix->set(3,2,3,4,12);
	hp_close_matrix->set(4,1,3,4,2);
	hp_close_matrix->set(4,3,3,4,7);
	hp_close_matrix->set(1,4,4,1,3);
	hp_close_matrix->set(2,3,4,1,17);
	hp_close_matrix->set(3,2,4,1,19);
	hp_close_matrix->set(3,4,4,1,3);
	hp_close_matrix->set(4,1,4,1,3);
	hp_close_matrix->set(4,3,4,1,6);
	hp_close_matrix->set(1,4,4,2,3);
	hp_close_matrix->set(2,3,4,2,14);
	hp_close_matrix->set(3,2,4,2,10);
	hp_close_matrix->set(3,4,4,2,3);
	hp_close_matrix->set(4,1,4,2,1);
	hp_close_matrix->set(4,3,4,2,1);
	hp_close_matrix->set(1,4,4,3,6);
	hp_close_matrix->set(2,3,4,3,18);
	hp_close_matrix->set(3,2,4,3,22);
	hp_close_matrix->set(3,4,4,3,4);
	hp_close_matrix->set(4,1,4,3,5);
	hp_close_matrix->set(4,3,4,3,6);
	hp_close_matrix->set(1,4,4,4,11);
	hp_close_matrix->set(2,3,4,4,20);
	hp_close_matrix->set(3,2,4,4,15);
	hp_close_matrix->set(3,4,4,4,11);
	hp_close_matrix->set(4,1,4,4,8);
	hp_close_matrix->set(4,3,4,4,8);


	internal_loop_matrix->init(scoreType(0));
	internal_loop_matrix->set(1,4,1,1,-7);
	internal_loop_matrix->set(3,4,1,1,-7);
	internal_loop_matrix->set(4,1,1,1,-7);
	internal_loop_matrix->set(4,3,1,1,-7);
	internal_loop_matrix->set(1,4,1,2,-7);
	internal_loop_matrix->set(3,4,1,2,-7);
	internal_loop_matrix->set(4,1,1,2,-7);
	internal_loop_matrix->set(4,3,1,2,-7);
	internal_loop_matrix->set(1,4,1,3,4);
	internal_loop_matrix->set(2,3,1,3,11);
	internal_loop_matrix->set(3,2,1,3,11);
	internal_loop_matrix->set(3,4,1,3,4);
	internal_loop_matrix->set(4,1,1,3,4);
	internal_loop_matrix->set(4,3,1,3,4);
	internal_loop_matrix->set(1,4,1,4,-7);
	internal_loop_matrix->set(3,4,1,4,-7);
	internal_loop_matrix->set(4,1,1,4,-7);
	internal_loop_matrix->set(4,3,1,4,-7);
	internal_loop_matrix->set(1,4,2,1,-7);
	internal_loop_matrix->set(3,4,2,1,-7);
	internal_loop_matrix->set(4,1,2,1,-7);
	internal_loop_matrix->set(4,3,2,1,-7);
	internal_loop_matrix->set(1,4,2,2,-7);
	internal_loop_matrix->set(3,4,2,2,-7);
	internal_loop_matrix->set(4,1,2,2,-7);
	internal_loop_matrix->set(4,3,2,2,-7);
	internal_loop_matrix->set(1,4,2,3,-7);
	internal_loop_matrix->set(3,4,2,3,-7);
	internal_loop_matrix->set(4,1,2,3,-7);
	internal_loop_matrix->set(4,3,2,3,-7);
	internal_loop_matrix->set(1,4,2,4,-7);
	internal_loop_matrix->set(3,4,2,4,-7);
	internal_loop_matrix->set(4,1,2,4,-7);
	internal_loop_matrix->set(4,3,2,4,-7);
	internal_loop_matrix->set(1,4,3,1,4);
	internal_loop_matrix->set(2,3,3,1,11);
	internal_loop_matrix->set(3,2,3,1,11);
	internal_loop_matrix->set(3,4,3,1,4);
	internal_loop_matrix->set(4,1,3,1,4);
	internal_loop_matrix->set(4,3,3,1,4);
	internal_loop_matrix->set(1,4,3,2,-7);
	internal_loop_matrix->set(3,4,3,2,-7);
	internal_loop_matrix->set(4,1,3,2,-7);
	internal_loop_matrix->set(4,3,3,2,-7);
	internal_loop_matrix->set(1,4,3,3,-7);
	internal_loop_matrix->set(3,4,3,3,-7);
	internal_loop_matrix->set(4,1,3,3,-7);
	internal_loop_matrix->set(4,3,3,3,-7);
	internal_loop_matrix->set(1,4,3,4,-7);
	internal_loop_matrix->set(3,4,3,4,-7);
	internal_loop_matrix->set(4,1,3,4,-7);
	internal_loop_matrix->set(4,3,3,4,-7);
	internal_loop_matrix->set(1,4,4,1,-7);
	internal_loop_matrix->set(3,4,4,1,-7);
	internal_loop_matrix->set(4,1,4,1,-7);
	internal_loop_matrix->set(4,3,4,1,-7);
	internal_loop_matrix->set(1,4,4,2,-7);
	internal_loop_matrix->set(3,4,4,2,-7);
	internal_loop_matrix->set(4,1,4,2,-7);
	internal_loop_matrix->set(4,3,4,2,-7);
	internal_loop_matrix->set(1,4,4,3,-7);
	internal_loop_matrix->set(3,4,4,3,-7);
	internal_loop_matrix->set(4,1,4,3,-7);
	internal_loop_matrix->set(4,3,4,3,-7);
	internal_loop_matrix->set(2,3,4,4,7);
	internal_loop_matrix->set(3,2,4,4,7);


	// Setting up the loop length cost tables
	int index=0;
	index =  0; hpLength[index] =   0; bulgeLength[index] =   0; ilLength[index] =   0;
	index =  1; hpLength[index] = -57; bulgeLength[index] = -38; ilLength[index] = -17;
	index =  2; hpLength[index] = -57; bulgeLength[index] = -28; ilLength[index] = -17;
	index =  3; hpLength[index] = -57; bulgeLength[index] = -32; ilLength[index] = -17;
	index =  4; hpLength[index] = -56; bulgeLength[index] = -36; ilLength[index] = -17;
	index =  5; hpLength[index] = -56; bulgeLength[index] = -40; ilLength[index] = -18;
	index =  6; hpLength[index] = -54; bulgeLength[index] = -44; ilLength[index] = -20;
	index =  7; hpLength[index] = -59; bulgeLength[index] = -46; ilLength[index] = -22;
	index =  8; hpLength[index] = -56; bulgeLength[index] = -47; ilLength[index] = -23;
	index =  9; hpLength[index] = -64; bulgeLength[index] = -48; ilLength[index] = -24;
	index = 10; hpLength[index] = -65; bulgeLength[index] = -49; ilLength[index] = -25;
	index = 11; hpLength[index] = -66; bulgeLength[index] = -50; ilLength[index] = -26;
	index = 12; hpLength[index] = -67; bulgeLength[index] = -51; ilLength[index] = -27;
	index = 13; hpLength[index] = -68; bulgeLength[index] = -52; ilLength[index] = -28;
	index = 14; hpLength[index] = -69; bulgeLength[index] = -53; ilLength[index] = -29;
	index = 15; hpLength[index] = -69; bulgeLength[index] = -54; ilLength[index] = -30;
	index = 16; hpLength[index] = -70; bulgeLength[index] = -54; ilLength[index] = -30;
	index = 17; hpLength[index] = -71; bulgeLength[index] = -55; ilLength[index] = -31;
	index = 18; hpLength[index] = -71; bulgeLength[index] = -55; ilLength[index] = -31;
	index = 19; hpLength[index] = -72; bulgeLength[index] = -56; ilLength[index] = -32;
	index = 20; hpLength[index] = -72; bulgeLength[index] = -57; ilLength[index] = -33;
	index = 21; hpLength[index] = -73; bulgeLength[index] = -57; ilLength[index] = -33;
	index = 22; hpLength[index] = -73; bulgeLength[index] = -58; ilLength[index] = -34;
	index = 23; hpLength[index] = -74; bulgeLength[index] = -58; ilLength[index] = -34;
	index = 24; hpLength[index] = -74; bulgeLength[index] = -58; ilLength[index] = -34;
	index = 25; hpLength[index] = -75; bulgeLength[index] = -59; ilLength[index] = -35;
	index = 26; hpLength[index] = -75; bulgeLength[index] = -59; ilLength[index] = -35;
	index = 27; hpLength[index] = -75; bulgeLength[index] = -60; ilLength[index] = -36;
	index = 28; hpLength[index] = -76; bulgeLength[index] = -60; ilLength[index] = -36;
	index = 29; hpLength[index] = -76; bulgeLength[index] = -60; ilLength[index] = -36;
	index = 30; hpLength[index] = -77; bulgeLength[index] = -61; ilLength[index] = -37;

}

inline void scorematrix::setNameGaps(std::string name,
				scoreType gapOpen, scoreType gapElongation,
				scoreType bpGapOpen, scoreType bpGapElongation) {

	matrixname = name;

	gap_open = gapOpen;
	elongation = gapElongation;
	for(int i = 1; i < 5; i++) {
		i_matrix->set(0, i, gap_open);
		i_matrix->set(i, 0, gap_open);
	}

	elongation_cost = elongation - gap_open;

	bp_gap_open = bpGapOpen;
	bp_elongation = bpGapElongation;

	bp_elongation_cost = bp_elongation - bp_gap_open;

//	prune_table = updatePruneTable(pruneStart, linearPrune, loopTableLength, seqLength1, seqLength2);
}

inline void scorematrix::setImatrix(
			const scoreType aa,
			const scoreType ac,
			const scoreType ag,
			const scoreType au,
			const scoreType cc,
			const scoreType cg,
			const scoreType cu,
			const scoreType gg,
			const scoreType gu,
			const scoreType uu,
			const scoreType gap) {
		i_matrix->set(1,1,aa);
		i_matrix->set(2,1,ac);
		i_matrix->set(3,1,ag);
		i_matrix->set(4,1,au);
		i_matrix->set(1,2,ac);
		i_matrix->set(2,2,cc);
		i_matrix->set(3,2,cg);
		i_matrix->set(4,2,cu);
		i_matrix->set(1,3,ag);
		i_matrix->set(2,3,cg);
		i_matrix->set(3,3,gg);
		i_matrix->set(4,3,gu);
		i_matrix->set(1,4,au);
		i_matrix->set(2,4,cu);
		i_matrix->set(3,4,gu);
		i_matrix->set(4,4,uu);
		i_matrix->set(1,0,gap_open);
		i_matrix->set(2,0,gap_open);
		i_matrix->set(3,0,gap_open);
		i_matrix->set(4,0,gap_open);
		i_matrix->set(5,0,gap_open);
		i_matrix->set(0,0,gap_open);
		i_matrix->set(0,1,gap_open);
		i_matrix->set(0,2,gap_open);
		i_matrix->set(0,3,gap_open);
		i_matrix->set(0,4,gap_open);
		i_matrix->set(0,5,gap_open);
}

inline void scorematrix::setS_matrix(
			const scoreType auau,
			const scoreType aucg,
			const scoreType augc,
			const scoreType augu,
			const scoreType auua,
			const scoreType auug,
			const scoreType cgcg,
			const scoreType cggc,
			const scoreType cggu,
			const scoreType cgua,
			const scoreType cgug,
			const scoreType gcgc,
			const scoreType gcgu,
			const scoreType gcua,
			const scoreType gcug,
			const scoreType gugu,
			const scoreType guua,
			const scoreType guug,
			const scoreType uaua,
			const scoreType uaug,
			const scoreType ugug) {
	s_matrix->set(1,4,1,4,auau);
	s_matrix->set(2,3,1,4,aucg);
	s_matrix->set(3,2,1,4,augc);
	s_matrix->set(3,4,1,4,augu);
	s_matrix->set(4,1,1,4,auua);
	s_matrix->set(4,3,1,4,auug);
	s_matrix->set(1,4,2,3,aucg);
	s_matrix->set(2,3,2,3,cgcg);
	s_matrix->set(3,2,2,3,cggc);
	s_matrix->set(3,4,2,3,cggu);
	s_matrix->set(4,1,2,3,cgua);
	s_matrix->set(4,3,2,3,cgug);
	s_matrix->set(1,4,3,2,augc);
	s_matrix->set(2,3,3,2,cggc);
	s_matrix->set(3,2,3,2,gcgc);
	s_matrix->set(3,4,3,2,gcgu);
	s_matrix->set(4,1,3,2,gcua);
	s_matrix->set(4,3,3,2,gcug);
	s_matrix->set(1,4,3,4,augu);
	s_matrix->set(2,3,3,4,cggu);
	s_matrix->set(3,2,3,4,gcgu);
	s_matrix->set(3,4,3,4,gugu);
	s_matrix->set(4,1,3,4,guua);
	s_matrix->set(4,3,3,4,guug);
	s_matrix->set(1,4,4,1,auua);
	s_matrix->set(2,3,4,1,cgua);
	s_matrix->set(3,2,4,1,gcua);
	s_matrix->set(3,4,4,1,guua);
	s_matrix->set(4,1,4,1,uaua);
	s_matrix->set(4,3,4,1,uaug);
	s_matrix->set(1,4,4,3,auug);
	s_matrix->set(2,3,4,3,cgug);
	s_matrix->set(3,2,4,3,gcug);
	s_matrix->set(3,4,4,3,guug);
	s_matrix->set(4,1,4,3,uaug);
	s_matrix->set(4,3,4,3,ugug);
}

#endif /* SCOREMATRIX */
