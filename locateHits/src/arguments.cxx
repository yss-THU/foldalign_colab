#ifndef ARGUMENTS
#define ARGUMENTS
#include<string>
#include<iostream>
#include<stdlib.h>

#include "nohit.hxx"
#include "helper.cxx"

/******************************************************************************
*                                                                             *
*   Copyright 2004 Jakob Hull Havgaard, hull@bioinf.kvl.dk                    *
*                                                                             *
*   This file is part of NOHIT                                            *
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

/*! \brief This class is intended to handle the parameters given to the main
	function

   The option parameters are set in the top of arguments function
*/

class arguments {
public:
//Set up and pass the parameters
	inline arguments(int argc, char* argv[]);

// Copy constructor
	inline arguments(const arguments& arg);
	
// Assignment constructor
	inline arguments& operator=(const arguments& arg);
	
// These functions return the value of argument.
// One function pr type of argument
	inline bool boolOpt(std::string opt);
	inline positionType ptOpt(std::string opt);
	inline lengthType ltOpt(std::string opt);
	inline scoreType stOpt(std::string opt);
	inline std::string stringOpt(std::string opt);

// Return the default value of an option
	inline std::string getDefault(std::string name);

// Functions used for returning information about the options
	inline int numberOfOptions() const {return num_arg;};
	inline std::string optionNumber(const int num) const;
	inline std::string descriptionNumber(const int num) const;

// Returns the argument indicated by num
	inline std::string argument(const int num) const;

// Returns the number of arguments
	inline int numberOfArguments() const;

// Set the value of an option
	inline void setSt(const std::string opt, const scoreType val);
	inline void setPt(const std::string opt, const positionType val);
	inline void setLt(const std::string opt, const lengthType val);
	inline void setString(const std::string opt, const std::string val);
	inline void setBool(const std::string opt, const bool val);

// Add a new parameter
	inline void addSt(const std::string opt, const scoreType val);
	inline void addPt(const std::string opt, const positionType val);
	inline void addLt(const std::string opt, const lengthType val);
	inline void addBool(const std::string opt, const bool val);
	inline void addString(const std::string opt, const std::string val);

// Check for the existance of an argument of the given name
	inline bool existStOpt(const std::string& opt) const;
	inline bool existPtOpt(const std::string& opt) const;
	inline bool existLtOpt(const std::string& opt) const;
	inline bool existStringOpt(const std::string& opt) const;
	inline bool existBoolOpt(const std::string& opt) const;

// Print all options and describtions
   inline void printOptions(const int length_name = 20, const int length_desc = 59) const;

// Cleans up
	inline ~arguments();

private:
	int num_arg; // The number of possible arguments. The size of options, defaults and types
	int n_st; // Number of possible scoreType arguments
	int n_pt; // Number of possible positionType arguments
	int n_lt; // Number of possible lengthType arguments
	int n_bool; // Number of possible bool arguments
	int n_string; // Number of possible string arguments
	std::string* st_opt; // Array of scoreType option names
	scoreType* st_val;         // Array with the scoreType option values
	std::string* pt_opt;   // etc
	positionType* pt_val;  // etc
	std::string* lt_opt;   // etc
	lengthType* lt_val;    // etc
	std::string* bool_opt;// Etc bool
	bool* bool_val;       // Etc bool
	std::string* string_opt; // Etc string
	std::string* string_val; // Etc string
	std::string* arg;		// The arguments
	int n_arg;				// The number of arguments
	std::string* options; // The name of the option
	std::string* value;   // The default value of the option
	std::string* types;   // The type of the option. int, bool and string allowed
	std::string* description; // The description of the option
	// This function defines the options
	inline void setOptions(int& n, std::string*& opt, std::string*& val, std::string*& typ, std::string*& des);
	inline void copyassign(const  arguments& a); // The copy constructor in a form the assignment operator can use
	inline void clean(); // The destructor in a form that can be used by the assignment operator
	template<class type> void assignNums(std::string name, std::string*& opt, type*& val, int argc, char* argv[], bool* notFound, int& c_type, int n_type);
	template<class type> type getOpt(std::string opt, std::string*& type_opt, type*& type_val, int type_size);
	template<class type> inline void set(const std::string opt, const type val, std::string*& type_opt, type*& type_val, int type_size);
	template<class type> inline void add(const std::string opt, const type val, std::string*& type_opt, type*& type_val, int& type_size);
	inline bool existOpt(const std::string& opt, std::string* type_opt, const int& type_size) const;

};

inline void arguments::setOptions(int& n, std::string*& opt, std::string*& val, std::string*& typ, std::string*& des) {
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
	n = 8;
	opt = new std::string[n];
// This array contains the name of the option
	opt[0] = "-score_cut_off";
	opt[1] = "-lambda_border";
	opt[2] = "-max_number_of_hits";
	opt[3] = "-estimate_parameters";
	opt[4] = "-version";
	opt[5] = "-help";
	opt[6] = "-h";
	val = new std::string[n];
// This array contains the default value
	val[0] = "0";
	val[1] = "false";
	val[2] = "1000";
	val[3] = "false";
	val[4] = "false";
	val[5] = "false";
	val[6] = "false";
	typ = new std::string[n];
// This array defines the type of the option
// There three possible types: bool, int and string
// There are five types: string, bool, scoreType, positionType, lengthType
	typ[0] ="scoreType";
	typ[1] ="bool";
	typ[2] ="lengthType";
	typ[3] ="bool";
	typ[4] ="bool";
	typ[5] ="bool";
	typ[6] ="bool";
	des = new std::string[n];
// This array holds the description of the option
// used by the help option.
	des[0] = "Alignments with scores lower than this are ignored. Default value: " + val[0];
	des[1] = "If set then alignments cannot start less than lambda nucleotides from the end of the sequences. Furthermore an alignment cannot overlap a better alignment or a nucleotide in a distance of lambda nucleotides from the better alignment. Default value: " + val[1];
	des[2] = "Output this number or less non-overlapping alignments for each FOLDALIGN alignment. Default: 1,000. Set to -1 to get all alignments. This can take significantly longer time";
	des[3] = "If used then the parameters of the extreme value distribution (lambda and K) are estimated from the non-overlapping alignments. One set of parameters are estimated for each alignment. The default is to use the build in values";
	des[4] = "Prints the version number and exits.";
	des[5] = "Prints this text,";
	des[6] = "also prints this text.";
//************************************************************
// End of option setting section.                            *
// Nothing below this point needs changing when options are  *
// added or changed                                          *
//************************************************************
}

// Constructor which parses the main arguments
inline arguments::arguments(int argc, char* argv[]) {

// Define all allowed options and setup the arrays
	setOptions(num_arg, options, value, types, description);

// Array used to check if an option has been handled
// needed because bool option only take one value from argv
// and int and string takes two.
	bool *notFound = new bool[argc];
	for(int j =0; j < argc; j++) {notFound[j] = true;}

// Count the number of different types of options
	n_st = n_pt = n_lt = n_bool = n_string = 0;
	for(int i=0; i < num_arg; i++) { // At the begining nothing has been found
		if (types[i].compare("scoreType")==0)         {n_st++;}
		else if (types[i].compare("positionType")==0) {n_pt++;}
		else if (types[i].compare("lengthType")==0)   {n_lt++;}
		else if (types[i].compare("bool")==0)         {n_bool++;}
		else if (types[i].compare("string")==0)       {n_string++;}
	}

// Allocate memory for the different types of options
	st_opt = new std::string[n_st];
	st_val = new scoreType[n_st];
	pt_opt = new std::string[n_pt];
	pt_val = new positionType[n_pt];
	lt_opt = new std::string[n_lt];
	lt_val = new lengthType[n_lt];
	bool_opt = new std::string[n_bool];
	bool_val = new bool[n_bool];
	string_opt = new std::string[n_string];
	string_val = new std::string[n_string];

//************************
// Handle the bool options

// Counter for how many bool arg were found.
// Later use to calculate the number of input files
	int c_bool = 0;
	int next_j=0; // count how far into the option array the count has reached
	for(int i=0; i < n_bool; i++) {
// Setup the default bool options from the arrays
		for(int j=next_j; j < num_arg; j++) {
			if (types[j].compare("bool")==0) {
				bool_opt[i] = options[j];
				if (value[j] .compare( "true")==0) {
					bool_val[i] = true;
				}
				else {
					bool_val[i] = false;
				}
				next_j = j +1;
				break;
			}
		}
// Check the commandline inputs and change the value if necessary
		for(int j =0; j < argc; j++) { // and check if it has been set
			if (bool_opt[i].compare(argv[j])==0) {
				bool_val[i] = !bool_val[i]; // If set change value to the opposite of the default value
				notFound[j] = false;
				c_bool++;
				break;
			}
		}
	}

// Handle the types
//	std::string type = "scoreType";
	int c_st = 0;
	int c_pt = 0;
	int c_lt = 0;
	assignNums("scoreType", st_opt, st_val, argc, argv, notFound, c_st, n_st);
	assignNums("positionType", pt_opt, pt_val, argc, argv, notFound, c_pt, n_pt);
	assignNums("lengthType", lt_opt, lt_val, argc, argv, notFound, c_lt, n_lt);
	
// Handle the string options
	int c_string = 0;
	next_j=0;
	for(int i=0; i < n_string; i++) {
// Setup the default values
		for(int j=next_j; j < num_arg; j++) {
			if (types[j].compare("string")==0) {
				string_opt[i] = options[j];
				string_val[i] = value[j];
				next_j = j +1;
				break;
			}
		}
// Handle the commandline options
		for(int j =0; j < argc; j++) { // and check if it has been set
			if (string_opt[i].compare(argv[j])==0) {
				string_val[i] = argv[j+1];
				notFound[j] = false;
				notFound[j+1] = false;
				c_string++;
				break;
			}
		}
	}
// Calculate the number of files
	n_arg = argc - c_bool - 2*c_st - 2*c_pt - 2*c_lt - 2*c_string-1;
// If any filenames is given - store them
	if (n_arg > 0) {
		arg = new std::string[n_arg];
		for(int i = (argc-n_arg); i < argc; i++) {
			if (argv[i][0] != '-') { // if its starts with a '-' then its not a file but an option
				arg[i - argc + n_arg] = argv[i];
				notFound[i] = false;
			}
		}
	}
// Check for unknown options			
	bool exitFlag = false;
	for(int j =1; j < argc; j++) {
		if(notFound[j]) {
			std::cerr << "Unknown option: " << argv[j] << std::endl;
			exitFlag=true;
		}
	}
	if (exitFlag) {throw -1; exit(-1);}
	delete[] notFound;
}


template<class type>
inline void arguments::assignNums(std::string name, std::string*& opt, type*& val, int argc, char* argv[], bool* notFound, int& c_type, int n_type) {
// Handle type options
	int next_j=0; // How far into the option array the count has come
	for(int i=0; i < n_type; i++) {
// Setup the default int options
		for(int j=next_j; j < num_arg; j++) {
			if (types[j].compare(name)==0) {
				opt[i] = options[j];
				const char* c_string = value[j].c_str();// new char[len+1];
				val[i] = atoi(c_string);
				next_j = j +1;
				break;
			}
		}
// Check the commandline options
		for(int j =1; j < argc; j++) { // and check if it has been set
			if (opt[i].compare(argv[j])==0) {
				val[i] = atoi(argv[j+1]);
				notFound[j] = false;
				notFound[j+1] = false;
				c_type++;
				break;
			}
		}
	}

}

// Copy constructor

inline arguments::arguments(const arguments& a) {
	copyassign(a);
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
	// Assign the new values
	setOptions(num_arg, options, value, types, description);
	n_st = a.n_st;
	n_pt = a.n_pt;
	n_lt = a.n_lt;
	n_bool = a.n_bool;
	n_string = a.n_string;
	n_arg = a.n_arg;
	if (n_arg > 0) {arg = new std::string[n_arg];}
	helper::copyArray(arg, a.arg, n_arg);
//	for(int i=0; i<n_arg; i++) {arg[i] = a.arg[i];}
	st_opt = new std::string[n_st];
	st_val = new scoreType[n_st];
	helper::copyArray(st_opt, a.st_opt, n_st);
	helper::copyArray(st_val, a.st_val, n_st);
	pt_opt = new std::string[n_pt];
	pt_val = new positionType[n_pt];
	helper::copyArray(pt_opt, a.pt_opt, n_pt);
	helper::copyArray(pt_val, a.pt_val, n_pt);
	lt_opt = new std::string[n_lt];
	lt_val = new lengthType[n_lt];
	helper::copyArray(lt_opt, a.lt_opt, n_lt);
	helper::copyArray(lt_val, a.lt_val, n_lt);
	string_opt = new std::string[n_string];
	string_val = new std::string[n_string];
	helper::copyArray(string_opt, a.string_opt, n_string);
	helper::copyArray(string_val, a.string_val, n_string);
	bool_opt = new std::string[n_bool];
	bool_val = new bool[n_bool];
	helper::copyArray(bool_opt, a.bool_opt, n_bool);
	helper::copyArray(bool_val, a.bool_val, n_bool);
}


inline std::string arguments::optionNumber(const int num) const {
	if (num >= num_arg) {
		std::string error = "Option number too high. Program error!";
		std::cerr << error << std::endl;
		throw error;
	}
	return options[num];
}


inline std::string arguments::descriptionNumber(const int num) const {
	if (num >= num_arg) {
		std::string error = "Description number too high. Program error!";
		std::cerr << error << std::endl;
		throw error;
	}
	return description[num];
}


inline bool arguments::boolOpt(std::string opt) {
	return getOpt(opt, bool_opt, bool_val, n_bool);
}
	

inline scoreType arguments::stOpt(std::string opt) {
	return getOpt(opt, st_opt, st_val, n_st);
}
	

inline positionType arguments::ptOpt(std::string opt) {
	return getOpt(opt, pt_opt, pt_val, n_pt);
}
	

inline lengthType arguments::ltOpt(std::string opt) {
	return getOpt(opt, lt_opt, lt_val, n_lt);
}
	

inline std::string arguments::stringOpt(std::string opt) {
	return getOpt(opt, string_opt, string_val, n_string);
}
	

template<class type>
inline type arguments::getOpt(std::string opt, std::string*& type_opt, type*& type_val, int type_size) {
	for(int i =0; i<type_size; i++) {
		if (opt.compare(type_opt[i])==0) {
			return type_val[i];
		}
	}
	std::string error = "Program error! Option: " + opt + " does not exists.";
	std::cerr << error << std::endl;
	throw error;
}


// Return the default argument of an option

inline std::string arguments::getDefault(std::string name) {
	for(int i = 0; i < num_arg; i++) {
		if (options[i].compare(name) == 0) {return value[i];}
	}
	std::cerr << "Program error: No default value for option: " << name << std::endl;
	throw -1;
}


inline std::string arguments::argument(const int num) const {
	return arg[num];
}


inline int arguments::numberOfArguments() const {
	return n_arg;
}
	
/**************************************************************
* These functions are used to change the value of an argument *
***************************************************************/


template<class type>
inline void arguments::set(const std::string opt, const type val, std::string*& type_opt, type*& type_val, int type_size) {
	bool found_flag = false; // False - option not found. True - option found
	for(int i=0; i<type_size; i++) { // Loop through all int options
		if (opt.compare(type_opt[i]) == 0) {
			type_val[i] = val;
			found_flag = true;
			break;
		}
	}
	if (!found_flag) {
		std::string error = "Program error! Option not found: ";
		std::cerr << error << opt << std::endl;
		throw error;
	}
}


inline void arguments::setSt(const std::string opt, const scoreType val) {
	set(opt, val, st_opt, st_val, n_st);
}


inline void arguments::setPt(const std::string opt, const positionType val) {
	set(opt, val, pt_opt, pt_val, n_pt);
}


inline void arguments::setLt(const std::string opt, const lengthType val) {
	set(opt, val, lt_opt, lt_val, n_lt);
}


inline void arguments::setString(const std::string opt, const std::string val) {
	set(opt, val, string_opt, string_val, n_string);
}


inline void arguments::setBool(const std::string opt, const bool val) {
	set(opt, val, bool_opt, bool_val, n_bool);
}


/************************************************************
* The functions adds a new parameter to the list of options *
*************************************************************/


template<class type>
inline void arguments::add(const std::string opt, const type val, std::string*& type_opt, type*& type_val, int& type_size) {
	// Make temporary arrays and store the old values
	std::string* tmp_opt = new std::string[type_size];
	type* tmp_val = new type[type_size];
	helper::copyArray(tmp_opt, type_opt, type_size);
	helper::copyArray(tmp_val, type_val, type_size);

	// Delete the old arrays and make new bigger ones
	delete[] type_opt;
	delete[] type_val;
	type_opt = new std::string[type_size+1];
	type_val = new type[type_size+1];
	
	// Copy the temporary arrays back into the permanet arrays
	helper::copyArray(type_opt, tmp_opt, type_size);
	helper::copyArray(type_val, tmp_val, type_size);

	// Delete the tmp arrays
	delete[] tmp_opt;
	delete[] tmp_val;
	
	// Add the new value and raise the size counter
	type_opt[type_size] = opt;
	type_val[type_size] = val;
	type_size++;
}
	

inline void arguments::addSt(const std::string opt, const scoreType val) {
	add(opt, val, st_opt, st_val, n_st);
}
	

inline void arguments::addPt(const std::string opt, const positionType val) {
	add(opt, val, pt_opt, pt_val, n_pt);
}
	

inline void arguments::addLt(const std::string opt, const lengthType val) {
	add(opt, val, lt_opt, lt_val, n_lt);
}
	

inline void arguments::addString(const std::string opt, const std::string val) {
	add(opt, val, string_opt, string_val, n_string);
}
	

inline void arguments::addBool(const std::string opt, const bool val) {
	add(opt, val, bool_opt, bool_val, n_bool);
}



/****************************************************************************
* The functions check the existance of an option of the given type and name *
*****************************************************************************/


inline bool arguments::existOpt(const std::string& opt, std::string* type_opt, const int& type_size) const {
	for(int i = 0; i < type_size; i++) {
		if (!type_opt[i].compare(opt)) {return true;}
	}
	return false;
}


inline bool arguments::existStOpt(const std::string& opt) const {
	return existOpt(opt, st_opt, n_st);
}


inline bool arguments::existPtOpt(const std::string& opt) const {
	return existOpt(opt, pt_opt, n_pt);
}


inline bool arguments::existLtOpt(const std::string& opt) const {
	return existOpt(opt, lt_opt, n_lt);
}


inline bool arguments::existStringOpt(const std::string& opt) const {
	return existOpt(opt, string_opt, n_string);
}


inline bool arguments::existBoolOpt(const std::string& opt) const {
	return existOpt(opt, bool_opt, n_bool);
}

/**********************************************************************
* This function prints all the options and their descriptions to cout *
**********************************************************************/


inline void arguments::printOptions(const int length_name, const int length_desc) const {
	for (int i=0; i < num_arg; i++) {
		const int len_name = options[i].length();
		std::cout << options[i];
		helper::print_space((length_name - len_name));
		helper::print_line(description[i], length_name, length_desc);
	}
}


inline void arguments::clean() {
	if (n_arg > 0) {delete[] arg;}
	delete[] options;
	delete[] value;
	delete[] types;
	delete[] description;
	delete[] st_opt;
	delete[] st_val;
	delete[] pt_opt;
	delete[] pt_val;
	delete[] lt_opt;
	delete[] lt_val;
	delete[] bool_opt;
	delete[] bool_val;
	delete[] string_opt;
	delete[] string_val;
}

/*******************************
* Destructor - nothing special *
********************************/


inline arguments::~arguments() {
	clean();
}

#endif /* ARGUMENTS */
