#ifndef WRITECOL
#define WRITECOL
#include "arguments.cxx"
#include "sequence.cxx"
#include "helper.cxx"
#include "stack.cxx"
#include "foldalign.hxx"
#include "output.hxx"
#include "output.interface.cxx"

#include <iostream>
#include <string>
#include <iomanip>

/******************************************************************************
*                                                                             *
*   Copyright 2004 -2007 Jakob Hull Havgaard, hull@bioinf.ku.dk               *
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


// Written by Jakob Hull Havgaard, 2004, hull@bioinf.kvl.dk.
class writecol : public outputInterface {
public:
	inline writecol(arguments& argu, 
					  const sequence* const one, 
					  const sequence* const two, scorematrix& score) 
		: arg(argu),
		  seq_1(one),
		  seq_2(two),
		  s_matrix(score),
		  flip(arg.boolOpt("switch")) {};
	
	inline void head();
	
	inline void localscorehead();
	
	inline void parameters();

	inline void plotscoreSep() {writeStarLine();};
	
	inline void out(const float sim, const float similar, const float total,
								  char  sequence_1[], char sequence_2[],
								  char  structure[], const int total_length,
								  int align_start_1, int align_end_1,
								  int align_start_2, int align_end_2, int align_score,
								  int* org_pos_1, int* org_bp_1, int* ali_bp_1,
								  int* org_pos_2, int* org_bp_2, int* ali_bp_2);


	inline void backtrack(const scoreType& score, const scoreType& stm_score, 
	                      const positionType& i, const positionType& k,
								 const lengthType& Wi, const lengthType& Wk,
								 const stateType& state,
								 const lengthType& len1, const lengthType& len2,
								 const lengthType& len3, const lengthType& len4,
								 const std::string& explain) const;

	inline void backtrackPrintStem(const positionType& ei,
									       const positionType& ek,
											 const lengthType& eWi,
											 const lengthType& eWk,
									       const positionType& bi,
									       const positionType& bk,
											 const lengthType& bWi,
											 const lengthType& bWk,
											 const scoreType& score) const;

	inline void backtrackPrintMBL(const positionType& oi,
									      const positionType& ok,
											const lengthType& oWi,
											const lengthType& oWk,
									      const positionType& li,
									      const positionType& lk,
											const lengthType& lWi,
											const lengthType& lWk,
									      const positionType& ri,
									      const positionType& rk,
											const lengthType& rWi,
											const lengthType& rWk,
											const scoreType& score) const;

	inline void backtrackStart(const long& count) const {
		std::cout << "; BACKTRACK           Bifurcation " << count << std::endl;
	}

	inline void backtrackEnd() const {
		std::cout << "; BACKTRACK           Branch end" << std::endl;
	}
	
	inline void errorNoGlobal() const {
		std::cout << "; NOTE                No global alignment was found. This can either be due to" << std::endl;
		std::cout << "; NOTE                the pruning, or because no structural alignment exists." << std::endl;
		std::cout << "; NOTE                Option: -no_prune can used to turn off pruning." << std::endl;
	}
	inline void errorNoLocal() const {
		std::cout << "; NOTE                No structural alignment was found." << std::endl;
	}

	inline void saveOutput(const bool& plot_score) const {
	
		if ( plot_score ) {writeStarLine();}
		writeEqualLine();

		writeMinusLine();
		writeStarLine();
		writeMinusLine();
		writeStarLine();
	}

	inline void saveOutputError(const bool& plot_score,
											   const std::string& error) const {

		std::cout << "; ERROR               ";
		if (error.compare("")) {std::cout << error << std::endl;}
		else {std::cout << "Unknown" << std::endl;}
		saveOutput(plot_score);
	}

	inline std::string getLShead() const {return "LS ";}

private:

	inline void entry_head(const sequence* const seqCurrent, const int score, const int start, const int end, const int length, char seq[], int org_pos[], int org_bp[], int ali_bp[]);

	inline int makeSequences(char*& seq_1, char*& seq_2, char*& struc, 
									 int*& org_pos_1, int*& org_bp_1, int*& ali_bp_1,
									 stack<positionType>*& leftPos_1,
									 stack<positionType>*& leftBasepair_1,
									 int*& org_pos_2, int*& org_bp_2, int*& ali_bp_2,
									 stack<positionType>*& leftPos_2,
									 stack<positionType>*& leftBasepair_2);

	inline void outputAlign(std::string name_1, std::string name_2, char sequen_1[], char sequen_2[], char struc[], int length);

	inline void printLine(char line[], int start, int end, std::string head, std::string title);
	
	inline float identity(char sequence1[], char sequence2[], int len, float& similar, float& total);
	
	inline void pf(const std::string& field, const bool newline = false) {
		std::cout << std::left << std::setw(front_length) << field << std::right;
		if (newline) {std::cout << std::endl;}
	}

	inline void writeStarLine() const {
		std::cout <<"; ******************************************************************************" << std::endl;
	}

	inline void writeEqualLine() const {
		std::cout << "; ==============================================================================" << std::endl;
	}

	inline void writeMinusLine() const {
		std::cout << "; ------------------------------------------------------------------------------" << std::endl;
	}


	arguments& arg;
	const sequence* const seq_1; // Stores the sequences
	const sequence* const seq_2;
	scorematrix& s_matrix; // Contains the score matrixs
	
	const bool flip;
};

inline void writecol::head() {
	pf("; FOLDALIGN");
	std::cout << arg.stringOpt("version") << std::endl;
	helper::print_array(reference, reference_size, "; REFERENCE", front_length);
	pf("; ALIGNMENT_ID");
	std::cout << arg.stringOpt("-ID") << std::endl;
	pf("; ALIGNING");
	std::cout << seq_1->getName() << " against " << seq_2->getName() << std::endl;

	parameters();
	
	pf("; SEQUENCE_LENGTH");
	std::cout << seq_1->getName() << " = " << seq_1->getLength() << std::endl;
	pf("; SEQUENCE_LENGTH");
	std::cout << seq_2->getName() << " = " << seq_2->getLength() << std::endl;
//	std::cout << seq_2->getLength() << std::endl;
	pf("; GC_CONTENT_SEQ_1");
	std::cout << std::setprecision(2) << seq_1->getGCcontent() << std::endl;
	pf("; GC_CONTENT_SEQ_2");
	std::cout << std::setprecision(2) << seq_2->getGCcontent() << std::endl;

}

inline void writecol::localscorehead() {
	head();
	pf("; SEQUENCE_1_COMMENT");
	std::cout << seq_1->getComment() << std::endl;
	pf("; SEQUENCE_2_COMMENT");
	std::cout << seq_2->getComment() << std::endl;
/*	pf("; LENGTH_SEQUENCE_1");
	std::cout << seq_1->getLength() << std::endl;
	pf("; LENGTH_SEQUENCE_2");
	std::cout << seq_2->getLength() << std::endl;
	pf("; GC_CONTENT_SEQ_1");
	std::cout << std::setprecision(2) << seq_1->getGCcontent() << std::endl;
	pf("; GC_CONTENT_SEQ_2");
	std::cout << std::setprecision(2) << seq_2->getGCcontent() << std::endl;
*/
	if (seq_2->getFilename().compare(seq_1->getFilename())) {
		pf("; FILENAMES");
		std::cout << seq_1->getFilename() << " and ";
		std::cout << seq_2->getFilename() << std::endl;
	}
	else {
		pf("; FILENAME");
		std::cout << seq_1->getFilename() << std::endl;
	}

	parameters();

	pf("; TYPE");
	std::cout << "Foldalign_local_scores" << std::endl;
	pf("; COL 1");
	std::cout << "label" << std::endl;
	pf("; COL 2");
	std::cout << "Alignment_start_position_sequence_1" << std::endl;
	pf("; COL 3");
	std::cout << "Alignment_end_position_sequence_1" << std::endl;
	pf("; COL 4");
	std::cout << "Alignment_start_position_sequence_2" << std::endl;
	pf("; COL 5");
	std::cout << "Alignment_end_position_sequence_2" << std::endl;
	pf("; COL 6");
	std::cout << "Alignment_score" << std::endl;
	std::cout <<"; ------------------------------------------------------------------------------" << std::endl;
}

inline void writecol::out(const float sim, const float similar, const float total,
								  char  sequence_1[], char sequence_2[],
								  char  structure[], const int total_length,
								  int align_start_1, int align_end_1,
								  int align_start_2, int align_end_2, int align_score,
								  int* org_pos_1, int* org_bp_1, int* ali_bp_1,
								  int* org_pos_2, int* org_bp_2, int* ali_bp_2
								 ) {

	std::string local_name_1 = seq_1->getName().substr(0,length_name_field);
	std::string local_name_2 = seq_2->getName().substr(0,length_name_field);
	const int length_comment = 80 - front_length;

	// Print the comments
	pf("; SEQUENCE_1_COMMENT");
	std::cout << seq_1->getComment().substr(0,length_comment) << std::endl;
	pf("; SEQUENCE_2_COMMENT");
	std::cout << seq_2->getComment().substr(0,length_comment) << std::endl;

	pf("; ALIGN", true);

	// Print the score and idendity
	pf("; ALIGN");
	std::cout << "Score: " << align_score << std::endl;
	pf("; ALIGN");

	int pres = 2;
	if (sim >= 100) {pres=3;}
	std::cout << "Identity: " << std::setprecision(pres) << sim << std::setprecision(6) << " % ( " << similar << " / " << total << " )" << std::endl;

	pf("; ALIGN", true);

	// Print the start coordinates
	pf("; ALIGN");
	std::cout << std::left << std::setw(length_name_field) << local_name_1;
	std::cout << " Begin " << align_start_1 << std::endl;
	pf("; ALIGN");
	std::cout << std::left << std::setw(length_name_field) << local_name_2;
	std::cout << " Begin " << align_start_2 << std::endl;

	pf("; ALIGN", true);

	// Print the alignment
	outputAlign(local_name_1, local_name_2, sequence_1, sequence_2, structure, total_length);

	pf("; ALIGN", true);

	// Print the end positions
	pf("; ALIGN");
	std::cout << std::left << std::setw(length_name_field) << local_name_1;
	std::cout << " End " << align_end_1 << std::endl;
	pf("; ALIGN");
	std::cout << std::left << std::setw(length_name_field) << local_name_2;
	std::cout << " End " << align_end_2 << std::endl;

	std::cout <<"; ==============================================================================" << std::endl;

	// If more than just the summary is wanted print the rest.
	entry_head(seq_1, align_score, align_start_1, align_end_1, total_length, sequence_1, org_pos_1, org_bp_1, ali_bp_1);
	entry_head(seq_2, align_score, align_start_2, align_end_2, total_length, sequence_2, org_pos_2, org_bp_2, ali_bp_2);

}

inline void writecol::parameters() {
	pf("; PARAMETER");
	std::cout <<"max_length=" << arg.ltOpt("-max_length") << std::endl;
	pf("; PARAMETER");
	std::cout <<"max_diff=" << arg.ltOpt("-max_diff") << std::endl;
	pf("; PARAMETER");
	std::cout <<"min_loop=" << arg.ltOpt("-min_loop") << std::endl;
	pf("; PARAMETER");
	std::cout <<"score_matrix=";
	if (arg.stringOpt("-score_matrix").compare("<default>")) {
		std::cout << arg.stringOpt("-score_matrix") << std::endl;
	}
	else {
		std::cout << "<" << s_matrix.getName() << ">" << std::endl;
	}
//	pf("; PARAMETER");
//	std::cout << "seed_constraints=" << arg.stringOpt("-seed_constraints") << std::endl;
	pf("; PARAMETER");
	std::cout <<"nobranching=";
	if (arg.boolOpt("-nobranch")) {std::cout << "<true>" << std::endl;}
	else {std::cout << "<false>" << std::endl;}
	pf("; PARAMETER");
	std::cout <<"global=";
	if (arg.boolOpt("-global")) {std::cout << "<true>" << std::endl;}
	else {std::cout << "<false>" << std::endl;}
	pf("; PARAMETER");
	std::cout <<"use_global_pruning=";
	if (arg.boolOpt("-use_global_pruning")) {std::cout << "<true>" << std::endl;}
	else {std::cout << "<false>" << std::endl;}

	positionType i = arg.ptOpt("-i");
	positionType j = arg.ptOpt("-j");
	positionType k = arg.ptOpt("-k");
	positionType l = arg.ptOpt("-l");
	if (flip) {
		helper::swap(i,k);
		helper::swap(j,l);
	}

	pf("; PARAMETER");
	std::cout <<"i=" << i << std::endl;
	pf("; PARAMETER");
	std::cout <<"j=" << j << std::endl;
	pf("; PARAMETER");
	std::cout <<"k=" << k << std::endl;
	pf("; PARAMETER");
	std::cout <<"l=" << l << std::endl;

	pf("; PARAMETER");
	if (arg.boolOpt("-no_pruning")) {
		std::cout <<"no_prune=<true>" << std::endl;
	}
	else {
		std::cout <<"no_prune=<false>" << std::endl;
	}
	pf("; PARAMETER");
	std::cout <<"min_LS_score=" << arg.stOpt("-min_LS_score") << std::endl;
}

inline void writecol::entry_head(const sequence* const seqCurrent, const int score, const int start, const int end, const int length, char seq[], int org_pos[], int org_bp[], int ali_bp[]) {
	// Output an alignment entry

	pf("; TYPE");
	std::cout <<"RNA" << std::endl;
	pf("; COL 1");
	std::cout <<"label" << std::endl;
	pf("; COL 2");
	std::cout <<"residue" << std::endl;
	pf("; COL 3");
	std::cout <<"seqpos" << std::endl;
	pf("; COL 4");
	std::cout <<"alignpos" << std::endl;
	pf("; COL 5");
	std::cout <<"align_bp" << std::endl;
	pf("; COL 6");
	std::cout <<"seqpos_bp" << std::endl;
	pf("; ENTRY");
	std::cout << seqCurrent->getName() << std::endl;
	pf("; ALIGNMENT_ID");
	std::cout << arg.stringOpt("-ID") << std::endl;
	pf("; ALIGNMENT_LIST");
	std::cout << seq_1->getName() << " " << seq_2->getName() << std::endl;
	pf("; FOLDALIGN_SCORE");
	std::cout << score << std::endl;
	pf("; GROUP");
	std::cout << seqCurrent->getGroupNumber() << std::endl;
	pf("; FILENAME");
	std::cout << seqCurrent->getFilename() << std::endl;
	pf("; START_POSITION");
	std::cout << start << std::endl;
	pf("; END_POSITION");
	std::cout << end << std::endl;
	pf("; ALIGNMENT_SIZE");
	std::cout <<"2" << std::endl;
	pf("; ALIGNMENT_LENGTH");
	std::cout << length << std::endl;
	pf("; SEQUENCE_LENGTH");
	std::cout << seqCurrent->getLength() << std::endl;
	pf("; SEQUENCE_GC_CONTENT");
	std::cout << std::setprecision(2) << seqCurrent->getGCcontent() << std::endl;

	parameters();

	std::cout <<"; ------------------------------------------------------------------------------" << std::endl;
	for (int i=0; i<length; i++) {
		if (seq[i] != '-') {
			std::cout << "N";
			std::cout << std::setw(space) << seq[i];
			std::cout << std::setw(space) << org_pos[i];
			std::cout << std::setw(space) << (i+1);
			if (ali_bp[i] !=-1) {
				std::cout << std::setw(space) << (ali_bp[i]+1);
				std::cout << std::setw(space) << org_bp[i] << std::endl;
			}
			else {
				std::cout << std::setw(space) << ".";
				std::cout << std::setw(space) << "." << std::endl;
			}
		}
		else {
			std::cout << "G";
			std::cout << std::setw(space) << "-";
			std::cout << std::setw(space) << ".";
			std::cout << std::setw(space) << (i+1);
			std::cout << std::setw(space) << ".";
			std::cout << std::setw(space) << ".";
			std::cout << std::endl;
		}
	}
	std::cout <<"; ******************************************************************************" << std::endl;
}

inline void writecol::backtrack(const scoreType& score, const scoreType& stm_score, 
                      const positionType& i, const positionType& k,
							 const lengthType& Wi, const lengthType& Wk,
							 const stateType& state,
							 const lengthType& len1, const lengthType& len2,
							 const lengthType& len3, const lengthType& len4,
							 const std::string& explain) const {
	std::cout << "; BACKTRACK";
	helper::print_space(11);
	std::cout << int(score)  << " " << int(stm_score) << " (";
	std::cout << int(i) << " " << int(i + Wi) << ", ";
	std::cout << int(k) << " " << int(k + Wk) << ") "<< int(state) << " ";
	std::cout << int(len1) << " " << int(len2) << " ";
	std::cout << int(len3) << " " << int(len4) << " ";
	std::cout << explain << std::endl;
}
	
inline void writecol::backtrackPrintStem(const positionType& ei,
									       const positionType& ek,
											 const lengthType& eWi,
											 const lengthType& eWk,
									       const positionType& bi,
									       const positionType& bk,
											 const lengthType& bWi,
											 const lengthType& bWk,
											 const scoreType& score) const {
	
	
	helper::pc("; STEM END", ei, ek, eWi, eWk);
	helper::pc("; START", bi, bk, bWi, bWk);
	std::cout << "SCORE " << score << std::endl;

}

inline void writecol::backtrackPrintMBL(const positionType& oi,
								      const positionType& ok,
										const lengthType& oWi,
										const lengthType& oWk,
								      const positionType& li,
								      const positionType& lk,
										const lengthType& lWi,
										const lengthType& lWk,
								      const positionType& ri,
								      const positionType& rk,
										const lengthType& rWi,
										const lengthType& rWk,
										const scoreType& score) const {

	helper::pc("; MBL", oi, ok, oWi, oWk);
	helper::pc("LEFT", li, lk, lWi, lWk);
	helper::pc("RIGHT", ri, rk, rWi, rWk);
	std::cout << "SCORE " << score << std::endl;
}

inline void writecol::outputAlign(std::string name_1, std::string name_2, char sequen_1[], char sequen_2[], char struc[], int length) {
	// Output the summary alignment
	std::string head = "; ALIGN";

	const int num_lines = length/line_length; // Calculate the number of lines (-1). Integer divison intended
	for(int i=0; i < num_lines; i++) {
		printLine(sequen_1, (i*line_length), ((i+1)*line_length), head, name_1);
		printLine(struc, (i*line_length), ((i+1)*line_length), head, "Structure");
		printLine(sequen_2, (i*line_length), ((i+1)*line_length), head, name_2);
		pf(head, true);
	}
	if (num_lines*line_length != length) {
		printLine(sequen_1, (num_lines*line_length), length, head, name_1);
		printLine(struc, (num_lines*line_length), length, head, "Structure");
		printLine(sequen_2, (num_lines*line_length), length, head, name_2);
	}
}

inline void writecol::printLine(char line[], int start, int end, std::string head, std::string title) {

	// Output an alignment line
	int length = end - start +1; // The length of the line
	int num_space = length/10; // The number of spaces to be added
	
	pf(head);

	std::cout << std::left << std::setw(length_name_field) << title << std::right;

	int pos = start;
	int stop = start;
	for(int i =0; i<=num_space; i++) {
		pos = stop;
		stop+=10;
		if (stop > end) {stop = end;}
		if (pos < stop) {std::cout << ' ';}
		for(int j=pos; j<stop; j++) {
			std::cout << line[j];
		}
	}
	std::cout << std::endl;
}


#endif
