#ifndef SEQUENCEARG
#define SEQUENCEARG
#include "sequence.cxx"
#include "arguments.cxx"
#include "foldalign.hxx"

/***************************************************************
* This class extends the sequence class with an argument value *
* Written by Jakob Hull Havgaard, 2004 hull@bioinf.kvl.dk      *
****************************************************************/

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

class sequenceArg : public sequence {
public:
	inline sequenceArg(std::string& sek_name, std::string& sek, 
	                   scorematrix scores, arguments argu, int groupNum=1, 
							 std::string com = "", std::string filename="<commandline>")
							 : sequence(sek_name, sek, scores, groupNum, com, filename),
							   arg(argu), score(scores) {};

	inline sequenceArg& operator=(const sequenceArg& s) {
		if (this != &s) {
			arg = s.arg;
			*static_cast<sequence*>(this) = s;
		}
		return *this;
	}

	inline sequenceArg(sequenceArg& s) : sequence(s), arg(s.arg), score(s.score) {}

	inline arguments& getArguments() {return arg;}
	
	inline scorematrix& getScorematrix() {return score;}

private:
	arguments arg;
	scorematrix score;
};
#endif /* SEQUENARG */
