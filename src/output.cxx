#ifndef OUTPUT
#define OUTPUT
#include "arguments.cxx"
#include "sequence.cxx"
#include "helper.cxx"
#include "stack.cxx"
#include "foldalign.hxx"
#include "output.hxx"
#include "output.interface.cxx"
#include "output.stk.cxx"
#include "output.col.cxx"
#include "output.summary.cxx"


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
class output {
public:

	output(arguments& arg, const sequence* const seq_1,
	            const sequence* const seq_2, scorematrix& score)
				  : flip(arg.boolOpt("switch")) {
		if (!flip) {
			setOutputFormat(arg, seq_1, seq_2, score);
		}
		else {
			setOutputFormat(arg, seq_2, seq_1, score);
		}
	}



	void head() {outi->head();}

	void localscorehead() {outi->localscorehead();};

	void parameters() {outi->parameters();};

	void plotscoreSep() {outi->plotscoreSep();};

	void out(const float sim, const float similar, const float total,
				char  sequence_1[], char sequence_2[],
				char  structure[], const int total_length,
				int align_start_1, int align_end_1,
				int align_start_2, int align_end_2, int align_score,
				int* org_pos_1, int* org_bp_1, int* ali_bp_1,
				int* org_pos_2, int* org_bp_2, int* ali_bp_2) {
		if (!flip) {
			outi->out(sim, similar, total, sequence_1, sequence_2, structure,
					 	total_length, align_start_1, align_end_1,
						align_start_2, align_end_2, align_score,
					 	org_pos_1, org_bp_1, ali_bp_1,
					 	org_pos_2, org_bp_2, ali_bp_2);
		}
		else {
			outi->out(sim, similar, total, sequence_2, sequence_1, structure,
					 	total_length, align_start_2, align_end_2,
						align_start_1, align_end_1, align_score,
					 	org_pos_2, org_bp_2, ali_bp_2,
					 	org_pos_1, org_bp_1, ali_bp_1);
		}
	}


	void backtrack(const scoreType& score, const scoreType& stm_score,
	               const positionType& i, const positionType& k,
						const lengthType& Wi, const lengthType& Wk,
						const stateType& state,
						const lengthType& len1, const lengthType& len2,
						const lengthType& len3, const lengthType& len4,
						const std::string& explain) const {
		if (!flip) {
			outi->backtrack(score, stm_score, i, k, Wi, Wk, state, len1, len2, len3,
								len4, explain);
		}
		else {
			outi->backtrack(score, stm_score, k, i, Wk, Wi, state, len3, len4, len1,
								len2, explain);
		}
	}


	void backtrackPrintStem(const positionType& ei,
									const positionType& ek,
									const lengthType& eWi,
									const lengthType& eWk,
									const positionType& bi,
									const positionType& bk,
									const lengthType& bWi,
									const lengthType& bWk,
									const scoreType& score) const {
		if (!flip) {
			outi->backtrackPrintStem(ei, ek, eWi, eWk, bi, bk, bWi, bWk, score);
		}
		else {
			outi->backtrackPrintStem(ek, ei, eWk, eWi, bk, bi, bWk, bWi, score);
		}
	}

	void backtrackPrintMBL(const positionType& oi,
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
		if (!flip) {
			outi->backtrackPrintMBL(oi, ok, oWi, oWk, li, lk, lWi, lWk,
										  ri, rk, rWi, rWk, score);
		}
		else {
			outi->backtrackPrintMBL(ok, oi, oWk, oWi, lk, li, lWk, lWi,
										  rk, ri, rWk, rWi, score);
		}
	}

	void backtrackStart(const long& count) const {outi->backtrackStart(count);}

	void backtrackEnd() const {outi->backtrackEnd();};

	void errorNoGlobal() const {outi->errorNoGlobal();};
	void errorNoLocal() const {outi->errorNoLocal();};

	void saveOutput(const bool& plot_score = false) const {
		outi->saveOutput(plot_score);
	};

	void saveOutputError(const bool& plot_score = false,
								const std::string& error = "") const {
		outi->saveOutputError(plot_score, error);
	};

	std::string getLShead() const {return outi->getLShead();};

	~output() {delete outi;};
private:

	outputInterface* outi;

	const bool flip;

	void setOutputFormat(arguments& arg,
	                     const sequence* seq_1, const sequence* seq_2,
								scorematrix& score, bool useDef = false) {

		std::string format = arg.stringOpt("-output_format");
		if (useDef) {format = arg.getDefault("-output_format");}

		// Switch to uppercase
		for(unsigned int p=0; p<format.length(); p++) {format[p] = toupper(format[p]);}

		if (! format.compare("STOCKHOLM")) {
			outi = new writestk(arg, seq_1, seq_2, score);
		}
		else if (! format.compare("COLUMN")) {
			outi = new writecol(arg, seq_1, seq_2, score);
		}
		else if (! format.compare("SUMMARY")) {
			outi = new writesummary(arg, seq_1, seq_2, score);
		}
		else {

			if (useDef) {
				std::string error = program_name + ": Program error! Unknown default format:\n" + format;
				throw exception(error, true);
			}

			std::cerr << program_name << ": Unknown output format:\n" << format;
			std::cerr << std::endl << "The default format will used." << std::endl;
			setOutputFormat(arg, seq_1, seq_2, score, true);
		}


	}

	output(); // No default constructor

};

#endif
