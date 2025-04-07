#ifndef HELPER
#define HELPER
#include <iostream>

/******************************************************************************
*                                                                             *
*   Copyright 2004 Jakob Hull Havgaard, hull@bioinf.kvl.dk                    *
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

//************************************************************//
// This class implements a number of helper functions.        //
//                                                            //
// Written by Jakob Hull Havgaard, 2004, hull@bioinf.kvl.dk   //
//************************************************************//

/*! \brief A library of little helper functions. Meant to make life easier
*/

class helper {
public:
	// Assigns the values i,j,k,l to the array to
	static inline void assign(int to[], int i, int j, int k, int l);

	// Swaps one and two
	template<class type> static inline void swap(type& one, type& two);

	// Prints <number> spaces to stdout
	static inline void print_space(const int number);
	
	static inline void print_line(const std::string& line, const int& start, const int& length);

	// Prints an array to std error.
	template<class type> static inline void printArray(type array[], int size) {for(int i=0; i< size; i++) {std::cerr << i << " " << array[i] << std::endl;}}

	// Make a new array of the new size and copies the values from the old array
   template<class type> static inline void expandArray(type*& oldAr, int& oldSize, int newSize);

	template<class type> static inline void copyArray(type& to, const type& from, int size)
	{ for(int i = 0; i < size; i++) {to[i] = from[i];} }

	template<class type, class num> static inline void init_array(type array[], num size, type val)
	{ for(num i = 0; i < size; i++) {array[i] = val;} }
};

inline void helper::print_space(const int number) {
	// Print number spaces
	for(int c = 0; c < number; c++) {std::cout << " ";}
}

inline void helper::print_line(const std::string& line, const int& start, const int& length) {
	// Prints line. The string is printed with length characters per line, and
	// if possible line breaks are placed between words. With the exception of
	// the first line start spaces are printed before each line.
	
	// Start position of substring
	int pos = 0;
	// The length of the string
	const int tot_len = line.length();

	// The substr of line which is to be printed.
//	std::string substr;

	while (tot_len - pos > length) {
		// Print all subsequences except the last one

		// The current substr of line which is to be printed.
		std::string substr = line.substr(pos, length+1);
		// Look for a nice place cut
		int len = substr.find_last_of(" ", length);

		if (len > 0) {
			// A nice place to cut was found
			// change the substring to the end at the new cut position
			substr = line.substr(pos, len);
			// update pos to the position of the next substring
			pos += len+1;
		}
		else {
			// No nice place to cut was found
			// make room for a "-" at the of substring
			substr = line.substr(pos, length-1);
			substr += '-';
			// update pos to the position of the next substring
			pos += length-1;
		}

		// Print the substring
		std::cout << substr << std::endl;
//std::cerr << pos << " " << length << " " << tot_len << " " << len << std::endl;
		// Print the spaces which start the next line
		print_space(start);
	}

	// Print the last substring
	std::cout << line.substr(pos) << std::endl; 

}

template<class type> inline void helper::swap(type& one, type& two) {
	type tmp = one;
	one = two;
	two = tmp;
}

inline void helper::assign(int to[], int i, int j, int k, int l) {
	to[0] = i;
	to[1] = j;
	to[2] = k;
	to[3] = l;
}

template<class type> inline void helper::expandArray(type*& oldAr, int& oldSize, int newSize) {
	// Expands an array. In a crude but simple way.
	// elem must be of the type which the array stores.
	// the elem element is not touched.
	
	// The new max size of the array must be handled by the caller.


	type* newAr = new type[newSize];
	if (newSize < oldSize) {std::cerr << "Program error. New array size is smaller than the old array size." << std::endl; throw -1;}
	for(int i=0; i < oldSize; i++) {
		newAr[i] = oldAr[i];
	}
	delete[] oldAr;
	oldAr = newAr;
}

#endif /*HELPER*/
