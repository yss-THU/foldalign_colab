#ifndef EXCEPTION
#define EXCEPTION

#include <string>

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

class exception {

	// This is an attempt at better error handling.
	// This object holds an error message and a fatal or nonfatal indication
	
public:

	exception() : message("Unknown error encountered. Foldalign will try to continue."), die(false) {};

	exception(const std::string mess, const bool fatal = false) : message(mess), die(fatal) {};

	std::string getMessage() const {return message;}
	
	bool getFatal() const {return die;}

private:

	const std::string message;
	
	const bool die;
	
};

#endif /* EXCEPTION */
