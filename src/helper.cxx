#ifndef HELPER
#define HELPER

#include "foldalign.hxx"
#include "exception.cxx"

#include <iostream>
#include <iomanip>
#include <string>
#include <stdlib.h>

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

//************************************************************//
// This class implements a number of helper functions.        //
//                                                            //
// Written by Jakob Hull Havgaard, 2004, hull@bioinf.kvl.dk   //
//************************************************************//
class helper {
public:
  // Assigns the values i,j,k,l to the array to
  static inline void assign(int to[], int i, int j, int k, int l)
  {to[0] = i; to[1] = j; to[2] = k; to[3] = l;}

  // Swaps one and two
  template<class type> static inline void swap(type& one, type& two)
  {type tmp = one; one = two; two = tmp;}

  // Prints <number> spaces to stdout
  static inline void print_space(const int number)
  {for(int c = 0; c < number; c++) {std::cout << " ";}}
	 
  static inline void print_line(const std::string& line, const int& start,
				const int& length);
  static inline void print_array(const std::string array[], const int size,
				 const std::string& head, const int start);

  // Prints an array to std error.
  template<class type> static inline void printArray(type array[], int size) {
    for(int i=0; i< size; i++) {std::cerr << i << " " << array[i] << std::endl;}
  }

  // Prints a matrix to std error.
  template<class type> static inline void printMatrix(type**& matrix,
						      int size) {
    for(int i=0; i< size; i++) {
      for(int j=0; j< size; j++) {
	std::cerr << i << " " << j << " "<< matrix[i][j] << std::endl;
    } }
  }

  // Make a new array of the new size and copies the values from the old array
  template<class type, class index> static inline void expandArray(type*& oldAr,
						      index& oldSize,
						      index newSize);

  template<class type> static inline void copyArray(type& to, const type& from,
						    int size)
  { for(int i = 0; i < size; i++) {to[i] = from[i];} }

  template<class type, class num> static inline void init_array(type array[],
								num size,
								type val)
  { for(num i = 0; i < size; i++) {array[i] = val;} }


  static inline void pc(const std::string& start, const positionType i,
			const positionType k, const lengthType Wi,
			const lengthType Wk, const bool end = false) {
    if (start.compare("")) {
      std::cout << start << " ";
    }
    std::cout << i << " " << i + Wi << " " << k << " " << k + Wk;
    if (end) {
      std::cout << std::endl;
    }
    else {
      std::cout << " " << std::flush;
    }
  }

  static inline void pcerr(const std::string& start, const positionType i,
			   const positionType k, const lengthType Wi,
			   const lengthType Wk, const bool end = false) {
    if (start.compare("")) {std::cerr << start << " ";}
    std::cerr << i << " " << i + Wi << " " << k << " " << k + Wk;
    if (end) {std::cerr << std::endl;}
    else {std::cerr << " " << std::flush;}
  }

  static inline void findWhiteSpace(const std::string& line, int& first,
				    int& last, int offset = 0);

  static inline int nextSpace(const std::string& line,
			      const int& offset = 0);

  static inline int nextNotSpace(const std::string& line,
				 const int& offset = 0);

  static std::string findNextString(const std::string& line,
				    const int& start, int& last);

  template<class type> static inline type min(const type a, const type b) {
    if (a < b) {return a;}
    return b;
  }

  template<class type> static inline type max(const type a, const type b) {
    if (a > b) {return a;}
    return b;
  }
	
  static inline int getValue(int& prev, std::string& line) {
    std::string substr = findName(prev, line);
    if (substr.compare("none")) {return atoi(substr.c_str());}
    else {return big_neg;}
  }
	
  static inline std::string findName(int& prev, std::string& line) {

		int pos;
		int end_pos = -1;
		int len = line.length();
		for(pos=prev; pos < len; pos++) {
			if ((line[pos] != ' ') && (line[pos] != '\t')) {
				end_pos = pos+1;
				while ((line[end_pos] != ' ') &&
				       (line[end_pos] != '\t') && 
					   (end_pos < len)) {
					end_pos++;
				}
				break;
			}
		}
		if (end_pos == -1) {
			std::string error = "Could not parse line. Somethings wrong with this line\n";
			error +=	line;
			throw exception(error, false);
		}
		prev = end_pos+1;
		return line.substr(pos, (end_pos - pos));
	}

};

inline void helper::print_line(const std::string& line, const int& start,
			       const int& length) {
  // Prints line. The string is printed with length characters per line, and
  // if possible line breaks are placed between words. With the exception of
  // the first line start spaces are printed before each line.
	
  // Start position of substring
  int pos = 0;
  // The length of the string
  const int tot_len = line.length();

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

    // Print the spaces which start the next line
    print_space(start);
  }

  // Print the last substring
  std::cout << line.substr(pos) << std::endl; 
}

inline void helper::print_array(const std::string array[], const int size,
				const std::string& head, const int start) {
  for(int i=0; i < size; i++) {
    std::cout << std::left << std::setw(start) << head << std::right;
    std::cout << array[i] << std::endl;
  }
}

template<class type, class index> inline void helper::expandArray(type*& oldAr, 
				index& oldSize,
				index newSize) {
  // Expands an array. In a crude but simple way.
  // elem must be of the type which the array stores.
  // the elem element is not touched.
	
  // The new max size of the array must be handled by the caller.

  // Listen to me, Lucy Moran, you just listen. 
  // When the Tacoma Sperm Bank was looking for donors, naturally I applied.
  // It's my civic duty and I like whales.
  // A routine physical examination revealed that I'm sterile.
  // Sure I thought it meant that I didn't have to take a bath, but the doctors
  // told me the truth. They told me I can't have babies.
  // So what I wanna know now is why are you having one and how? 

#ifdef DEBUG
	// It is necessary to keep this check in DEBUG otherwise the compiler
	// will complain in the cases where newSize is proveable < oldSie
	if (newSize < oldSize) {
		std::string error = "Program error. New array size is smaller than the old array size.";
		throw exception(error, false);
	}
#endif

  type* newAr = new type[newSize];
  for(int i=0; i < oldSize; i++) {
    newAr[i] = oldAr[i];
  }
  delete[] oldAr;
  oldAr = newAr;
}

inline void helper::findWhiteSpace(const std::string& line, int& first,
				   int& last,  int offset) {
  // Find the first substring of space or tabs after character number offset
  // first becomes the start of the space sub string and last becomes the first
  // non tab or space position
  first = last = -1;
  first = nextSpace(line, offset);
  if (first == -1) {return;}
  last = first+1;
  last = nextNotSpace(line, last);
}

inline int helper::nextSpace(const std::string& line, const int& offset) {

  // Find the next space or tab position after offset
  // returns -1 if no position was found

  const int length = line.length();
  
  for(int i = offset; i < length; i++) {
    if (!line.substr(i,1).compare(" ") ||
	!line.substr(i,1).compare("\t")) {return i;}
  }
  return -1;
}

inline int helper::nextNotSpace(const std::string& line,const int& offset) {

  // Find the next non space or tab position after offset
  // returns -1 if no position was found

  const int length = line.length();

  for(int i = offset; i < length; i++) {
    if (line.substr(i,1).compare(" ") &&
	line.substr(i,1).compare("\t")) {return i;}
  }
  return -1;
}

inline std::string helper::findNextString(const std::string& line,
					  const int& start, int& last) {

  // If the character at position start is a ' (or ") then the function returns
  // the substring from start+1 until the position before the next ' (or ") (or 
  // the end of the string). Likewise for ". If the character at position start 
  // is not a ' or " then it is the substring from start untill the next space
  // or tab which is returned.
  std::string value;
  int pos = -1;
  last = -1;
  int s_pos = 0;

  const char let = char(line[start]);

  if (let == '\'') {
    pos = line.find("'", start+1);
    s_pos++;
  }
  else if (let == '"') {
    pos = line.find('"', start+1);
    s_pos++;
  }
  else {
    pos = nextSpace(line, start);
  }

  if (pos < 0) {value = line.substr(start+s_pos);}
  else {
    value = line.substr(start+1, pos - start -s_pos);
    last = nextNotSpace(line, pos+1);
  }
  return value;
}
	
#endif /*HELPER*/
