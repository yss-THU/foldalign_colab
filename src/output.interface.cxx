#ifndef OUTPUTINTERFACE
#define OUTPUTINTERFACE
#include "arguments.cxx"
#include "sequence.cxx"
#include "helper.cxx"
#include "stack.cxx"
#include "foldalign.hxx"
#include "output.hxx"

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
class outputInterface {
public:

	virtual void head() = 0;
	
	virtual void localscorehead() = 0;
	
	virtual void parameters() = 0;

	virtual void plotscoreSep() = 0;
	
	virtual void out(const float sim, const float similar, const float total,
								  char  sequence_1[], char sequence_2[],
								  char  structure[], const int total_length,
								  int align_start_1, int align_end_1,
								  int align_start_2, int align_end_2, int align_score,
								  int* org_pos_1, int* org_bp_1, int* ali_bp_1,
								  int* org_pos_2, int* org_bp_2, int* ali_bp_2) = 0;

	virtual void backtrack(const scoreType& score, const scoreType& stm_score, 
	                      const positionType& i, const positionType& k,
								 const lengthType& Wi, const lengthType& Wk,
								 const stateType& state,
								 const lengthType& len1, const lengthType& len2,
								 const lengthType& len3, const lengthType& len4,
								 const std::string& explain) const = 0;
	
	virtual void backtrackPrintStem(const positionType& ei,
									       const positionType& ek,
											 const lengthType& eWi,
											 const lengthType& eWk,
									       const positionType& bi,
									       const positionType& bk,
											 const lengthType& bWi,
											 const lengthType& bWk,
											 const scoreType& score) const = 0;

	virtual void backtrackPrintMBL(const positionType& oi,
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
											const scoreType& score) const = 0;

	virtual void backtrackStart(const long& count) const = 0;

	virtual void backtrackEnd() const = 0;

	virtual void errorNoGlobal() const = 0;
	virtual void errorNoLocal() const = 0;

	virtual void saveOutput(const bool& plot_score = false) const = 0;

	virtual void saveOutputError(const bool& plot_score = false,
									    const std::string& error = "") const = 0;	

	virtual std::string getLShead() const = 0;

	virtual ~outputInterface() {};

};

#endif
