#ifndef NOHITHXX
#define NOHITHXX

#include <string>
// This file defines the subrevison number. The file is made by the make file
// which tries to get the subversion revision number.
#include "../../src/revnumber.hxx"

/*! \brief Used for types which are of the same size as Foldalign scores
*/
typedef int scoreType;

/*! \brief Used for types which are of the same range as positions along the
	sequences/genomes
*/
typedef long posType;

/*! \brief Needed by arguments.cxx
*/
typedef posType positionType; // Needed by arguments.cxx

/*! \brief Used for window sizes along the sequences
*/
typedef int lengthType;


// Information about the program
const std::string program_name = "locateHits";
const std::string version = "2.5.2"; //rev is defined in revnumber.hxx
const std::string copyright = "Copyright by Jakob Hull Havgaard, 2011-2019";

/*! \brief This constant is mainly used as an indicator of a task being
	finished or impossible
*/
const scoreType bad = -1000;

#endif
