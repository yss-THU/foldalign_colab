#ifndef ARGUMENTS
#define ARGUMENTS

#include<string>
//#include<iostream>
//#include<iomanip>
//#include<stdlib.h>
#include<vector>

#include "options.cxx"

#include "foldalign.hxx"
#include "helper.cxx"
#include "exception.cxx"
#include "convertString.cxx"

/******************************************************************************
*                                                                             *
*   Copyright 2004 - 2012 Jakob Hull Havgaard, hull@bioinf.ku.dk              *
*                                                                             *
*   This file is part of FOLDALIGN                                            *
*                                                                             *
*   FOLDALIGN is free software; you can redistribute it and/or modify         *
*   it under the terms of the GNU General Public License as published by      *
*   the Free Software Foundation; either version 2 of the License, or         *
*   (at your option) any later version.                                       *
*                                                                             *
*   FOLDALIGN is distributed in the hope that it will be useful,              *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with FOLDALIGN; if not, write to the Free Software                  *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA *
*                                                                             *
******************************************************************************/

/* This class is intended to handle the parameters given to the main function
   Written by Jakob Hull Havgaard 2004 hull@bioinf.kvl.dk

   The option parameters are set in the top of arguments function
*/

class arguments {
public:
	//Set up and pass the parameters
	inline arguments(int argc, char* argv[]);

	// Copy constructor
	inline arguments(const arguments& arg) {copyassign(arg);};
	
	// Assignment constructor
	inline arguments& operator=(const arguments& arg);
	
	// These functions return the value of argument.
	// One function pr type of argument
	inline bool boolOpt(std::string opt) {return boolOpts->get(opt);}
	inline positionType ptOpt(std::string opt) {return positionOpt->get(opt);}
	inline lengthType ltOpt(std::string opt) {return lengthOpt->get(opt);}
	inline scoreType stOpt(std::string opt) {return scoreOpt->get(opt);}
	inline std::string stringOpt(std::string opt) {return stringOpts->get(opt);};

	// Return the default value of an option
	inline std::string getDefault(std::string name);

	// Functions used for returning information about the options
	inline int numberOfOptions() const {return num_arg;};
	inline std::string optionNumber(const int num) const;
	inline std::string descriptionNumber(const int num) const;

	// Returns the argument indicated by num
	inline std::string argument(const int num) const {return arg.at(num);};

	// Returns the number of arguments
	inline int numberOfArguments() const {return arg.size();};

	// Set the value of an option
	inline void setSt(const std::string& opt, const scoreType& val) {scoreOpt->set(opt, val);};
	inline void setPt(const std::string& opt, const positionType& val) {positionOpt->set(opt, val);};
	inline void setLt(const std::string& opt, const lengthType& val) {lengthOpt->set(opt, val);};
	inline void setString(const std::string& opt, const std::string& val) {stringOpts->set(opt, val);};
	inline void setBool(const std::string& opt, const bool& val) {boolOpts->set(opt, val);};

	// Add a new parameter
	inline void addSt(const std::string& opt, const scoreType& val) {scoreOpt->add(opt, val);};
	inline void addPt(const std::string& opt, const positionType& val) {positionOpt->add(opt, val);};
	inline void addLt(const std::string& opt, const lengthType& val) {lengthOpt->add(opt, val);};
	inline void addBool(const std::string& opt, const bool& val) {boolOpts->add(opt, val);};
	inline void addString(const std::string& opt, const std::string& val) {stringOpts->add(opt, val);};

	// Check for the existance of an argument of the given name
	inline bool existStOpt(const std::string& opt) const {return scoreOpt->exists(opt);};
	inline bool existPtOpt(const std::string& opt) const {return positionOpt->exists(opt);};
	inline bool existLtOpt(const std::string& opt) const {return lengthOpt->exists(opt);};
	inline bool existStringOpt(const std::string& opt) const {return stringOpts->exists(opt);};
	inline bool existBoolOpt(const std::string& opt) const {return boolOpts->exists(opt);};

	// Print all options and descriptions
   inline void printOptions() const;

	// Cleans up
	inline ~arguments() {clean();};

private:

	// number of characters pr line. Used when printing the help text.
	static const int line_length = 80;

	int num_arg; // The number of possible options. The size of options, defaults and types
//	int n_arg;				// The number of arguments
	std::string* optNames; // The name of the option
	std::string* value;   // The default value of the option
	std::string* types;   // The type of the option. int, bool and string allowed
	std::string* description; // The description of the option

	options<scoreType>* scoreOpt;
	options<positionType>* positionOpt;
	options<lengthType>* lengthOpt;
	options<bool>* boolOpts;
	options<std::string>* stringOpts;

	std::vector<std::string> arg;

	inline bool errorIfLastElement(const int j, const int argc,
									std::string option, std::string& error);

	// This function defines the options
	inline void setOptions();
	inline void copyassign(const  arguments& a); // The copy constructor in a form the assignment operator can use
	inline void clean(); // The destructor in a form that can be used by the assignment operator
};

inline void arguments::setOptions() {
	//************************************************************
	// This function defines the options                         *
	// When an option is added five things must be added/changed *
	// The int num must be set to the number of possible options *
	// When adding a new option add one to this number.          *
	// An opt line must be made. This line defines the name      *
	// of the option.                                            *
	// A val line must be added. This is the default value of    *
	// the option                                                *
	// A typ line indicating the type of the option must be      *
	// added. Possible types are bool, int and string.           *
	// Finally a description text should be written.             *
	// Remember of to use the correct index when writting these  *
	// lines                                                     *
	//************************************************************
	// This number must be updated when new options are added
	num_arg = 29;

	// This array contains the name of the option
	
	optNames = new std::string[num_arg];
	optNames[0] = "-max_diff";
	optNames[1] = "-max_length";
	optNames[2] = "-min_loop";
	optNames[3] = "-score_matrix";
//	optNames[4] = "-seed_constraints";
//	optNames[5] = "-min_seed_expand_length";
//	optNames[6] = "-backtrack_seed_constraints";
	optNames[4] = "-chunk_size";
	optNames[5] = "-no_pruning";
	optNames[6] = "-ID";
	optNames[7] = "-nobranch";
	optNames[8] = "-no_backtrack";
	optNames[9] = "-global";
	optNames[10] = "-use_global_pruning";
	optNames[11] = "-i";
	optNames[12] = "-j";
	optNames[13] = "-k";
	optNames[14] = "-l";
	optNames[15] = "-format";
	optNames[16] = "-output_format";
	optNames[17] = "-plot_score";
//	optNames[21] = "-print_seed_constraints";
	optNames[18] = "-min_LS_score";
	optNames[19] = "-print_all_LS_scores";
	optNames[20] = "-number_of_processors";
	optNames[21] = "-align_self";
	optNames[22] = "-backtrack_info";
	optNames[23] = "-all_scores";
	optNames[24] = "-memory_roof";
	optNames[25] = "-memory_info";
	optNames[26] = "-version";
	optNames[27] = "-help";
	optNames[28] = "-h";

	// This array contains the default value

	value = new std::string[num_arg];
	value[0] = "25";
	value[1] = "-1";
	value[2] = "3";
	value[3] = "<default>";
//	value[4] = "<none>";
//	value[5] = "-1";
//	value[6] = "<none>";
	value[4] = "-1";
	value[5] = "false";
	value[6] = "n.a.";
	value[7] = "false";
	value[8] = "false";
	value[9] = "false";
	value[10] = "false";
	value[11] = "1";
	value[12] = "-1";
	value[13] = "1";
	value[14] = "-1";
	value[15] = "fasta";
	value[16] = "column";
	value[17] = "false";
//	value[21] = "false";
	value[18] = "0";
	value[19] = "false";
	value[20] = "1";
	value[21] = "false";
	value[22] = "false";
	value[23] = "false";
	value[24] = "-1";
	value[25] = "false";
	value[26] = "false";
	value[27] = "false";
	value[28] = "false";

	// This array defines the type of the option
	// There three possible types: bool, int and string
	// There are five types: string, bool, scoreType, positionType, lengthType

	types = new std::string[num_arg];
	types[0] ="lengthType";
	types[1] ="lengthType";
	types[2] ="lengthType";
	types[3] ="string";
//	types[4] ="string";
//	types[5] ="lengthType";
//	types[6] ="string";
	types[4] ="lengthType";
	types[5] ="bool";
	types[6] ="string";
	types[7] ="bool";
	types[8] ="bool";
	types[9] ="bool";
	types[10] ="bool";
	types[11] ="positionType";
	types[12] ="positionType";
	types[13] ="positionType";
	types[14] ="positionType";
	types[15] ="string";
	types[16] ="string";
	types[17] ="bool";
//	types[21] ="bool";
	types[18] ="scoreType";
	types[19] ="bool";
	types[20] ="lengthType";
	types[21] ="bool";
	types[22] ="bool";
	types[23] ="bool";
	types[24] ="string";
	types[25] ="bool";
	types[26] ="bool";
	types[27] ="bool";
	types[28] ="bool";

	// This array holds the description of the option
	// used by the help option. Note that the default value is sometimes
	// added at the end of the description string. You must update the index
	// numbers of these when you insert af new option.

	description = new std::string[num_arg];
	description[0] = "The maximum length difference between two sequences being aligned. Also known as delta. High values will make FOLDALIGN much slower and use much more memory. Default value: 25.";
	description[1] = "The maximum length of the alignment, not counting gaps. Also known as lambda. Default value: Sequence length.";
	description[2] = "The minimum hairpin loop length. Default value 3.";
	description[3] = "Path and filename of a score matrix file. There are two buildt in matrices one for local alignment (the default) and one for global alignment. The global alignment matrix is used when option -global is used.";
//	description[4] = "Path and filename of a seed constraints file. If this option is used only hairpin loop stems which coordinates are given in the file can aligned. An alignment is not forced to go through all the given coordinates, but it must pass through at least one set. Note: This option can only be used safely when only one set of sequences are aligned, or when the pair format is used and the constraints file is specified as an option for each of sequence pairs. Default: " + value[4];
//	description[5] = "Seeds with a length shorter than this are not expanded into longer alignments. They can however be part of bifurcated loops. Use -1 to indicate no minimum length. Default value: " + value[5];
//	description[6] = "Path and filename of a backtrack seed constraints file. If this option is used then the seed constraints (see the -seed_constraints option) in this file is used during the backtrack of the stem segments. This is of use when local alignments are used iteratively to build global alignments. If this option is not used then the -seed_constraints are used (if there are any).";
	description[4] = "A memory control option. This option controls the length of the sections which one of the sequences is split into. Low values make FOLDALIGN use less memory but run slower. High values make FOLDALIGN use more memory but run faster. Chunk_size must be at least -max_length. Default value: -max_length.";
	description[5] = "Turns off pruning";
	description[6] = "Alignment id. Default value: n.a.";
	description[7] = "Turns FOLDALIGN into a stemloop only algorithm. This speeds up the program. Default: off.";
	description[8] = "Turns off backtracking. The best local alignment scores will be print as if option -plot_score had been used.";
	description[9] = "Do global alignments. When this option is used the default score matrix is switched to the default global score matrix.";
	description[10] = "Use global alignment pruning during local alignment. During global alignment the pruning cut off is compensated for the minimum number of gaps needed to align two sequences. Use this option makes the same compensations during local alignment. Default: Not used";
	description[11] = "Start aligning from this position in the first sequence. Default: Position one.";
	description[12] = "Stop aligning at this position in the first sequence. Default: End of the sequence.";
	description[13] = "Start aligning from this position in the second sequence. Default: Position one.";
	description[14] = "Stop aligning at this position in the second sequence. Default: End of the sequence.";
	description[15] = "Specifies the input format. Allowed values: fasta, tab, pair, stockholm, and commandline. Default: fasta.";
	description[16] = "Specifies the output format. Allowed values: column, summary and stockholm. Default: "+value[17];
	description[17] = "Print the best local structural alignment score and its coordinates at all pair of positions between the two sequences for which a structural alignment was found.";
//	description[21] = "Print positions, score and state for subalignments which can be used as seed_constraints.";
	description[18] = "Local structural alignments with a score below this cut off are not printed, unless option -print_all_LS_scores is used. Default value: " + value[22];
	description[19] = "Also print a dummy alignment score for those coordinates where no structural alignment was found. If option -plot_score is not used then this option has no effect.";
	description[20] = "Sets the number of processors (threads) which the program uses. Default: "+value[24];
	description[21] = "When set, a sequence is also aligned against it self (this is only relevant if only one fasta or tab file is used as input). Default value: off.";
	description[22] = "Adds extra backtrack information to the output. General format score, positions, state. Default value: off.";
	description[23] = "Dumps all the scores and states before storage in the matrixes. Warning: Using this option will generate an extreme amount of output.";
	description[24] = "Currently not working properly. When this option is used the program will try to estimate its memory consumption. If the estimated memory consumption is above the option's value the alignment of the current pair of sequences is abandoned. One of the letter M, G or T can used to indicate mega, giga, or tetra bytes. The default value is no limit.";
	description[25] = "Print the estimated memory consumption";
	description[26] = "Prints the version number and exits.";
	description[27] = "Prints this text,";
	description[28] = "also prints this text.";
	
	//***********************************************************
	// End of option setting section.                            *
	// Nothing below this point needs changing when options are  *
	// added or changed                                          *
	//************************************************************
}

// Constructor which parses the main arguments
inline arguments::arguments(int argc, char* argv[]) {

	// Define all allowed options and setup the arrays
	setOptions();
 
	scoreOpt = new options<scoreType>("scoreType", num_arg, optNames,
						value, types, description,
						&convertString::toNumeric< scoreType > );

	positionOpt = new options<positionType>("positionType", num_arg, optNames,
						value, types, description,
						&convertString::toNumeric<positionType>);

	lengthOpt = new options<lengthType>("lengthType", num_arg, optNames,
						value, types, description,
						&convertString::toNumeric<lengthType>);

	boolOpts = new options<bool>("bool", num_arg, optNames,
						value, types, description,
						&convertString::toBool);

	stringOpts = new options<std::string>("string", num_arg, optNames,
						value, types, description,
						&convertString::toString);

	std::string error = "";

	for(int j = 1; j < argc; j++) {

		std::string argOption = std::string(argv[j]);
		std::string argValue = "";
		if (j < argc -1) {
			argValue = std::string(argv[j+1]);
		}
		
		if (argv[j][0] != '-') {
			// This must be a file name
			arg.push_back(argOption);
		}
		else if (boolOpts->exists(argOption)) {
			bool current = boolOpts->get(argOption);
			boolOpts->set(argOption, !current);
		}
		else if (scoreOpt->exists(argOption)) {
			if (errorIfLastElement(j, argc, argOption, error)) { continue; }
			scoreOpt->set(argOption, convertString::toNumeric<scoreType>(argValue));
			j++; // j+1 is the value to option j
		}
		else if (positionOpt->exists(argOption)) {
			if (errorIfLastElement(j, argc, argOption, error)) { continue; }
			positionOpt->set(argOption, convertString::toNumeric<positionType>(argValue));
			j++; // j+1 is the value to option j
		}		
		else if (lengthOpt->exists(argOption)) {
			if (errorIfLastElement(j, argc, argOption, error)) { continue; }
			lengthOpt->set(argOption, convertString::toNumeric<lengthType>(argValue));
			j++; // j+1 is the value to option j
		}
		else if (stringOpts->exists(argOption)) {
			if (errorIfLastElement(j, argc, argOption, error)) { continue; }
			stringOpts->set(argOption, argValue);
			j++;
		}
		else {
			error += " Unknown option: " + argOption;
		}
	}
	
	if (error.compare("") != 0) {
		throw exception("Could not parse the commandline:"+error, true);
	}
}

inline bool arguments::errorIfLastElement(const int j, const int argc,
							std::string option, std::string& error) {

	if (j == argc -1) {
		error += " option: " + option + " is missing its value";
		return true;
	}
	return false;
}

// Assignment constructor

inline arguments& arguments::operator=(const arguments& a) {
	if (this != &a) { // If not the same object
		// Clean up
		clean();
		// Assign the new values
		copyassign(a);
	}
	return *this;
}


inline void arguments::copyassign(const  arguments& a) {

	scoreOpt = new options<scoreType>(*a.scoreOpt);
	positionOpt = new options<positionType>(*a.positionOpt);
	lengthOpt = new options<lengthType>(*a.lengthOpt);
	boolOpts = new options<bool>(*a.boolOpts);
	stringOpts = new options<std::string>(*a.stringOpts);

	arg = a.arg;
	num_arg = a.num_arg;
	optNames = new std::string[num_arg];
	value = new std::string[num_arg];
	types = new std::string[num_arg];
	description = new std::string[num_arg];
	
	helper::copyArray(optNames, a.optNames, num_arg);
	helper::copyArray(value, a.value, num_arg);
	helper::copyArray(types, a.types, num_arg);
	helper::copyArray(description, a.description, num_arg);
}


inline std::string arguments::optionNumber(const int num) const {
	if (num >= num_arg) {
		std::string error = "Option number too high. Program error!";
		throw exception(error, true);
	}
	return optNames[num];
}


inline std::string arguments::descriptionNumber(const int num) const {
	if (num >= num_arg) {
		std::string error = "Description number too high. Program error!";
		throw exception(error, true);
	}
	return description[num];
}


// Return the default argument of an option

inline std::string arguments::getDefault(std::string name) {
	for(int i = 0; i < num_arg; i++) {
		if (optNames[i].compare(name) == 0) {return value[i];}
	}
	std::string error = "Program error! No default value for option: " + name;
	throw exception(error, true);
}
	

/**********************************************************************
* This function prints all the options and their descriptions to cout *
**********************************************************************/


inline void arguments::printOptions() const {
	unsigned int length_name = 0;
	for (int i=0; i < num_arg; i++) {
		if (optNames[i].length() > length_name) {length_name = optNames[i].length();}
	}
	length_name++;
	const int length_desc = line_length - length_name;
	for (int i=0; i < num_arg; i++) {
		std::cout << std::left << std::setw(length_name) << optNames[i] << std::right;
		helper::print_line(description[i], length_name, length_desc);
	}
}


inline void arguments::clean() {

	delete scoreOpt;
	delete positionOpt;
	delete lengthOpt;
	delete boolOpts;
	delete stringOpts;
	delete[] optNames;
	delete[] value;
	delete[] types;
	delete[] description;

}

#endif /* ARGUMENTS */
