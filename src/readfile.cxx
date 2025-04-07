#ifndef READFILE
#define READFILE

#include "exception.cxx"

#include <string>
#include <fstream>
#include <iostream>

class readfile {
	

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
	/**************************************
	*                                     *
	* This class is used for reading.     *
	* It can either read a named file,    *
	* or read from stdin.                 *
	*                                     *
	**************************************/


public:

	readfile( std::string fname ="" ) : 
	   filename("<STDIN>"),
		file(&std::cin),
		must_del(false) {

		// By default files are read from stdin
		// stdin values are set above.

		if (fname.compare("")) {
			// Read a name file

			// Store the name
			filename = fname;

			// Open the new file
			try {
				file = new std::ifstream(filename.c_str());
			}
			catch (...) {
				std::string error = "Trouble opening the file: " + filename;
				throw exception(error, true);
			}
			
			// Set the delete variable
			must_del = true;
		}

		// Check that the file was openend correctly
		if (file->fail()) {
			delete file;
			std::string error = "Could not open the file: " + filename;
			throw exception(error, true);
		}
	};
	
	// The named file case the file pointer must be deleted
	~readfile() {
		if (must_del) {
			delete file;
		}
	};

	void skip_to_empty() {

		// Ignore all line until an empty line is found

		std::string line;

		getline(*file, line);

		while (line.compare("") && (!file->fail())) {
			getline(*file,line);
		}
	
	}


	bool get_line(std::string& line){
		// Read a line from the file
		// Empty lines and comment lines starting with # are ignored.
		// The line is returned through the line argument.
		// The function returns true if a line could be read
		// and false if no line was read.
		
		getline(*(file), line);

		while (((line[0] == '#') || (line.length() == 0)) && (!file->fail())){
			getline(*(file), line);
		}

		if (file->fail()) {
			line = "";
			return false;
		}

		return true;
	};

	bool get_line_failEmpty(std::string& line){
		// Read a line from the file
		// Comment lines starting with # are ignored.
		// The line is returned through the line argument.
		// The function returns true if a line could be read
		// and false if no line was read or the line was empty.
		
		getline(*(file), line);

		while (line[0] == '#' && !file->fail()){
			getline(*(file), line);
		}

		if (file->fail() || !line.compare("")) {
			line = "";
			return false;
		}

		return true;
	};

	std::string name() {return filename;}

private:
	// The name of the file
	std::string filename;
	
	// Pointer to the input file stream. Must be deleted in the named file case
	std::istream* file;

	// True if file must be deleted, false otherwise (i.e. it is cin)
	bool must_del;	

	// A getline which also remove \r at the end of the line
	void getline(std::istream& stream, std::string& line) {
		std::getline(stream, line);
		if (*line.rbegin() == '\r') {
                    line.erase(line.end()-1);
                }
	}
};

#endif /* readfile */
