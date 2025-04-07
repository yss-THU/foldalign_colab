#ifndef SEQUENCE
#define SEQUENCE
#include <string>
#include <iostream>

#include "foldalign.hxx"
#include "scorematrix.cxx"


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

//*********************************************************************
// This class holds a sequence. The sequence can be supplied directly *
// or read from a fasta file                                          *
// Written by Jakob Hull Havgaard 2004, hull@bioinf.kvl.dk            *
//*********************************************************************

class sequence { // Stores a RNA sequence with gaps.
public:
	sequence(const sequence& s);

	sequence& operator=(const sequence& s);
	
	// Store the sequence without reading a file
	inline sequence(std::string& sek_name, std::string& sek, scorematrix& score, int groupNum=1, std::string com = "", std::string filename="<commandline>");// Store the sequence

	// Return the value at position pos (position 1..len allowed)
	int getPos(const int pos) const {return (seq[pos]);};

	// Return the nucleotide at this position
	char getLetter(const int pos) const {return org_seq[pos-1];}

	// Return the length of the sequence
	int getLength() const {return len;};

	// Return the sequence name
	std::string getName() const {return name;}

	// Return the sequence comment
	std::string getComment() const {return comment;}

	// Return the sequence filename
	std::string getFilename() const {return fil;}

	// Return the sequence group number
	int getGroupNumber() const {return group;}
	
	// Return the gc content
	float getGCcontent() const {return gcContent;}

	~sequence() {delete[] seq;};

private:
	int len;  		// The sequence length
	std::string name; 	// Sequence name
	int* seq; 		// The sequence
	int group; 	// The sequence's group number
	std::string comment; 	// Store the sequence comment
	std::string fil; 	// Name of the file the sequence orginated from
	std::string org_seq; // The original sequence
	float gcContent;
	sequence(); 		//No default constructor
};

inline sequence::sequence(std::string& sek_name, std::string& sek, 
                          scorematrix& score, int groupNum, std::string com, 
								  std::string filename)
								  : len(sek.length()), name(sek_name), group(groupNum),
								  comment(com), fil(filename), org_seq(sek)  {
	// Store the sequence
	seq = new int[len+1]; // Set the size of seq
	seq[0] = 0; // The 0 position is a gap
	int gc_count = 0;
	for(int n=0; n < len; n++) { // For all positions in the sequence
		seq[n+1] = score.alfa(sek[n]);
		if (toupper(sek[n]) == 'C' || toupper(sek[n]) == 'G') {
			gc_count++;
		}
	}
	gcContent = float(gc_count)/float(len);
}

inline sequence& sequence::operator=(const sequence& s) {
	if (&s != this) {
		len = s.getLength();
		name = s.getName();
		group = s.getGroupNumber();
		comment = s.getComment();
		fil = s.getFilename();
		org_seq = s.org_seq;
		delete[] seq;
		seq = new int[len+1];
		for(int n=0; n <= len; n++) {
			seq[n] = s.getPos(n);
		}
		gcContent = s.gcContent;
	}
	return *this;
}

inline sequence::sequence(const sequence& s) : group(1) {
	len = s.getLength();
	name = s.getName();
	group = s.getGroupNumber();
	comment = s.getComment();
	fil = s.getFilename();
	org_seq = s.org_seq;
	seq = new int[len+1];
	for(int n=0; n <= len; n++) {
		seq[n] = s.getPos(n);
	}
	gcContent = s.gcContent;
}	
#endif /* SEQUENCE */
