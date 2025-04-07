#ifndef WRITESTK
#define WRITESTK
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
class writestk : public outputInterface {
public:
	inline writestk(arguments& argu, 
					  const sequence* const one, 
					  const sequence* const two, scorematrix& score) 
		: arg(argu),
		  seq_1(one),
		  seq_2(two),
		  s_matrix(score), head_written(false), flip(arg.boolOpt("switch")) {};
	
	inline void head();

	inline void localscorehead() {head();};
	
	inline void parameters();

	inline void plotscoreSep() {};

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
		writeGFnN("BACKTRACK", "Bifurcation ");
		std::cout << count << std::endl;
	}

	inline void backtrackEnd() const {writeGF("BACKTRACK", "Branch end");}

	inline void errorNoGlobal() const {
		writeGF("ERROR", "No global alignment was found. This can either be due to");
		writeGF("ERROR", "the pruning, or because no structural alignment exists.");
		writeGF("ERROR", "Option: -no_prune can used to turn off pruning.");
	}
	inline void errorNoLocal() const {
		writeGF("ERROR", "No structural alignment was found.");
	}
	
	inline void saveOutput(const bool& plot_score) const {
		std::cout << "//" << std::endl;
	}

	inline void saveOutputError(const bool& plot_score,
										 const std::string& error) const {
		if (error.compare("")) {writeGF("ERROR", error);}
		else {writeGF("ERROR", "Unknown");}
		saveOutput(plot_score);
	}

	inline std::string getLShead() const {return "#=GF LS ";}

private:

	void writeLineBool(std::string opt);
	template< class type > void writeGFnN(const std::string& head, const type& feature) const;
	template< class type > void writeGF(const std::string& head, const type& feature) const;
	template< class type > void writeGS(std::string name, std::string feature, type value);
	void printStkLine(std::string head, char feat[], const int length);

	arguments& arg;
	const sequence* const seq_1; // Stores the sequences
	const sequence* const seq_2;
	scorematrix& s_matrix; // Contains the score matrixs
	bool head_written;
	const bool flip;
};

inline void writestk::head() {

	if ( ! head_written ) {
		std::cout << "# STOCKHOLM 1.0" << std::endl;
		std::cout << std::endl;
		writeGF("VERSION", arg.stringOpt("version"));
		for(int i=0; i<reference_size; i++) {
			writeGF("REFERENCE", reference[i]);
		}
		writeGF("ID", arg.stringOpt("-ID"));
		writeGF("SEQUENCE_NAMES", seq_1->getName() + " against " + seq_2->getName());

		parameters();

		writeGS(seq_1->getName(), "COMMENT", seq_1->getComment());
		writeGS(seq_2->getName(), "COMMENT", seq_2->getComment());
		writeGS(seq_1->getName(), "LENGTH", seq_1->getLength());
		writeGS(seq_2->getName(), "LENGTH", seq_2->getLength());
		writeGS(seq_1->getName(), "GC_CONTENT", seq_1->getGCcontent());
		writeGS(seq_2->getName(), "GC_CONTENT", seq_2->getGCcontent());
		writeGS(seq_1->getName(), "INPUTFILE", seq_1->getFilename());
		writeGS(seq_2->getName(), "INPUTFILE", seq_2->getFilename());
		head_written = true;
	}
}

inline void writestk::out(const float sim, const float similar, const float total,
								  char  sequence_1[], char sequence_2[],
								  char  structure[], const int total_length,
								  int align_start_1, int align_end_1,
								  int align_start_2, int align_end_2, int align_score,
								  int* org_pos_1, int* org_bp_1, int* ali_bp_1,
								  int* org_pos_2, int* org_bp_2, int* ali_bp_2) {


	std::cout << std::endl;
	writeGS(seq_1->getName(), "START_POSITION", align_start_1);
	writeGS(seq_2->getName(), "START_POSITION", align_start_2);
	writeGS(seq_1->getName(), "END_POSITION", (align_end_1));
	writeGS(seq_2->getName(), "END_POSITION", (align_end_2));
	std::cout << std::endl;

	writeGF("SCORE", align_score);

	int pres = 2;
	if (sim >= 100) {pres=3;}
	std::cout << "#=GF " << std::setw(stk_feat_size) << "ALIGNMENT_IDENTITY" << " " << std::setprecision(pres) << sim << std::setprecision(6) << " % ( " << similar << " / " << total << " )" << std::endl;
	writeGF("ALIGNMENT_LENGTH", total_length);
	std::cout << std::endl;

	printStkLine(seq_1->getName(), sequence_1, total_length);
	printStkLine(seq_2->getName(), sequence_2, total_length);
	printStkLine("#=GC SS_cons ", structure, total_length);

	std::cout << std::endl;
	std::cout << "//" << std::endl;

}

inline void writestk::parameters() {
	writeGFnN("PARAMETER", "max_length=");
	std::cout << arg.ltOpt("-max_length") << std::endl;
	writeGFnN("PARAMETER", "max_diff=");
	std::cout << arg.ltOpt("-max_diff") << std::endl;
	writeGFnN("PARAMETER", "min_loop=");
	std::cout << arg.ltOpt("-min_loop") << std::endl;
	writeGF("PARAMETER", "score_matrix=" + arg.stringOpt("-score_matrix"));
//	writeGF("PARAMETER", "seed_constraints=" + arg.stringOpt("-seed_constraints"));

	writeLineBool("nobranch");
	writeLineBool("global");
	writeLineBool("use_global_pruning");

	positionType i = arg.ptOpt("-i");
	positionType j = arg.ptOpt("-j");
	positionType k = arg.ptOpt("-k");
	positionType l = arg.ptOpt("-l");
	if (flip) {
		helper::swap(i,k);
		helper::swap(j,l);
	}

	writeGFnN("PARAMETER", "i=");
	std::cout << i << std::endl;
	writeGFnN("PARAMETER", "j=");
	std::cout << j << std::endl;
	writeGFnN("PARAMETER", "k=");
	std::cout << k << std::endl;
	writeGFnN("PARAMETER", "l=");
	std::cout << l << std::endl;
	
	writeLineBool("no_pruning");

	writeGFnN("PARAMETER", "min_LS_score=");
	std::cout << arg.stOpt("-min_LS_score") << std::endl;
}
	
inline void writestk::backtrack(const scoreType& score, const scoreType& stm_score, 
                      const positionType& i, const positionType& k,
							 const lengthType& Wi, const lengthType& Wk,
							 const stateType& state,
							 const lengthType& len1, const lengthType& len2,
							 const lengthType& len3, const lengthType& len4,
							 const std::string& explain) const {
							 
	writeGFnN("BACKTRACK", "");

	std::cout << int(score)  << " " << int(stm_score) << " (";
	std::cout << int(i) << " " << int(i + Wi) << ", ";
	std::cout << int(k) << " " << int(k + Wk) << ") "<< int(state) << " ";
	std::cout << int(len1) << " " << int(len2) << " ";
	std::cout << int(len3) << " " << int(len4) << " ";
	std::cout << explain << std::endl;
}
	
inline void writestk::backtrackPrintStem(const positionType& ei,
									       const positionType& ek,
											 const lengthType& eWi,
											 const lengthType& eWk,
									       const positionType& bi,
									       const positionType& bk,
											 const lengthType& bWi,
											 const lengthType& bWk,
											 const scoreType& score) const {
	
	writeGFnN("STEM", "");
	helper::pc("END:", ei, ek, eWi, eWk);
	helper::pc("START:", bi, bk, bWi, bWk);
	std::cout << "SCORE: " << score << std::endl;

}

inline void writestk::backtrackPrintMBL(const positionType& oi,
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
	writeGFnN("MBL", "");
	helper::pc("", oi, ok, oWi, oWk);
	helper::pc("LEFT", li, lk, lWi, lWk);
	helper::pc("RIGHT", ri, rk, rWi, rWk);
	std::cout << "SCORE " << score << std::endl;
}

inline void writestk::writeLineBool(std::string opt) {
	std::string boolTrue = "<true>";
	std::string boolFalse = "<false>";
	
	if (arg.boolOpt("-"+opt)) {
		writeGF("PARAMETER", opt + "=" + boolTrue);
	}
	else {
		writeGF("PARAMETER", opt + "=" + boolFalse);
	}
}

template< class type >
inline void writestk::writeGFnN(const std::string& feature, const type& value) const {
	std::cout << "#=GF " << std::left << std::setw(stk_feat_size) << feature;
	std::cout << " " << value;
}

template< class type >
inline void writestk::writeGF(const std::string& feature, const type& value) const {
	writeGFnN(feature, value);
	std::cout << std::endl;
}

template< class type >
inline void writestk::writeGS(std::string name, std::string feature, type value) {
	std::cout << "#=GS " << std::left << std::setw(stk_name_size) << name << " ";
	std::cout << std::setw(stk_feat_size) << feature;
	std::cout << " " << value << std::endl;
}

inline void writestk::printStkLine(std::string head, char feat[], const int length) {

	std::cout << std::left<< std::setw(stk_feat_size + 5) << head << " " << std::right;
	for(int i= 0; i < length; i++) {std::cout << feat[i];}
	std::cout << std::endl;
}

#endif /* WRITESTK */

