#include "foldalign.hxx"
#include "local.cxx" // The local alignment algorithm
#include "backtrack.cxx" // The global alignment and backtrack algorithm
#include "sequence.cxx" // Stores the sequences
#include "scorematrix.cxx"
#include "seqs.cxx"
#include "arguments.cxx"
#include "helper.cxx"
#include "exception.cxx"
#include "output.stk.cxx"
#include "output.col.cxx"
#include "output.summary.cxx"
#include "cell.cxx"
#include "longCell.cxx"
#include "constraints.cxx"

#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <new>

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

constraints* setConstraints(std::string fileName, arguments& arg) {

	
	if (fileName.compare("<none>")) {
		std::string error = "Could not read the constraints file: ";
		error += fileName;

		arg.setBool("hpstart", false);
			
		constraints* cons = new constraints(fileName, arg.boolOpt("switch"));

		if (cons == 0) {
			throw exception(error, false);
		}
		
		return cons;
	}

	return 0;
}

template<class type> inline void deletePtr(type ptr) {
	if (ptr != 0) {
		delete ptr;
	}
}

void cleanUp(output*& out, constraints*& cons, constraints*& bt_cons) {
	deletePtr(out);
	deletePtr(cons);
	deletePtr(bt_cons);
}

int align(sequence* seq_1, sequence* seq_2, arguments arg, scorematrix& score) {

	//*****************************************************
	//
	// Check the coordinates and the sequence lengths
	//
	
	output* out = 0;
	constraints* cons = 0;
	constraints* bt_cons = 0;
	
	lengthType len_1 = seq_1 ->getLength(); // Short names for the length of the sequences
	lengthType len_2 = seq_2 ->getLength();

	positionType i = arg.ptOpt("-i");
	positionType j = arg.ptOpt("-j");
	positionType k = arg.ptOpt("-k");
	positionType l = arg.ptOpt("-l");

	if (len_1 < len_2) {
		helper::swap<sequence*>(seq_1, seq_2);
		arg.setBool("switch", true);
		helper::swap(len_1, len_2);
		helper::swap(i,k);
		helper::swap(j,l);
		arg.setPt("-i", i);
		arg.setPt("-j", j);
		arg.setPt("-k", k);
		arg.setPt("-l", l);
	}
		

	if ((i != 1) || (j != -1)) {
		if (i > seq_1->getLength()) {
			std::string error = "The start position given by option -i is higher than the length of the sequence";
			throw exception(error, false);
		}
		if (i < 1) {
			std::cerr << "The start position -i " << i << " must be 1 or larger. The start position has been set to 1" << std::endl;
			i = 1;
		}
		if ((j < i) && (j != -1)) {
			std::string error = "The end position given by option -j is less than the start position given by option -i";
			throw exception(error, false);
		}
		if (j > seq_1->getLength()) {
			std::cerr << "The end position -j " << j << " is higher than the length of sequence: " << seq_1->getName() << "'s length. The end position has been set to the sequence length: " << seq_1->getLength() << std::endl;
			j = -1;
		}
		if (j == -1) {j = seq_1->getLength(); arg.setPt("-j", j);}
		
		len_1 = j - i +1;
	}
	else {
		len_1 = seq_1->getLength(); // Short names for the length of the sequences
		arg.setPt("-j", len_1);
		j = len_1;
	}

	if ((k != 1) || (l != -1)) {
		if (k > seq_2->getLength()) {
			std::string error = "The start position given by option -k is higher than the length of the sequence.";
			throw exception(error, false);
		}
		if (k < 1) {
			std::cerr << "The start position -k " << k << " must be 1 or larger. The start position has been set to 1" << std::endl;
			k = 1;
		}
		if ((l < k) && (l != -1)) {
			std::string error = "The end position given by option -l is less than the start position given by option -k";
			throw exception(error, false);
		}
		if (l > seq_2->getLength()) {
			std::cerr << "The end position -l " << l << " is higher than the length of sequence: " << seq_2->getName() << "'s length. The end position has been set to the sequence length: " << seq_2->getLength() << std::endl;
			l = -1;
		}
		if (l == -1) {l = seq_2->getLength(); arg.setPt("-l", l);}
		len_2 = l - k +1;
	}
	else {
		len_2 = seq_2->getLength(); // Short names for the length of the sequences
		arg.setPt("-l", len_2);
		l = len_2;
	}

	if (len_1 == 0) {
		std::cerr << "Sequence " << seq_1->getName() << " is empty. Skipping" << std::endl;
		return -1;
	}
	if (len_2 == 0) {
		std::cerr << "Sequence " << seq_2->getName() << " is empty. Skipping" << std::endl;
		return -1;
	}

	lengthType real_max_length = j - i +1;
	lengthType real_min_length = l - k +1;
	if (real_min_length > real_max_length) {
		helper::swap(real_max_length, real_min_length);
	}

	//**********************************************************************
	//
	// Handle the paramters, check lambda, delta etc bounds due lengths etc

	// Handle lambda and delta < 1
	if (arg.boolOpt("-global") && arg.ltOpt("-max_length") > 0 && 
	    arg.ltOpt("-max_length") < real_max_length) {
		 	std::cerr << "Using option -global with a -max_length shorter than the ";
			std::cerr << "longest of the two sequence is not possible. ";
			std::cerr << "Please remove either option -global or -max_length. ";
			std::cerr << "Skipping sequences" << std::endl;
			return -1;
	}
	if (arg.ltOpt("-max_diff") < 1) {arg.setLt("-max_diff", real_max_length);}
	if (arg.ltOpt("-max_diff") > real_max_length) {
		arg.setLt("-max_diff", real_max_length);
	}
	if (arg.ltOpt("-max_length") < 1) {
		arg.setLt("-max_length", len_2+arg.ltOpt("-max_diff"));
	}
	if (arg.ltOpt("-max_length") > (len_2 + arg.ltOpt("-max_diff")) && 
	    !arg.boolOpt("-global") ) {
		 	arg.setLt("-max_length", len_2+arg.ltOpt("-max_diff"));
	}
	if (arg.ltOpt("-max_length") > real_max_length ) {
		arg.setLt("-max_length", real_max_length);
	}
	if ( arg.ltOpt("-min_loop") > arg.ltOpt("-max_length") ) {
		arg.setLt("-min_loop", arg.ltOpt("-max_length"));
	}
	if (arg.ltOpt("-chunk_size") < 1) {
		arg.setLt("-chunk_size",2*arg.ltOpt("-max_length"));
	}
	if (arg.ltOpt("-chunk_size") < arg.ltOpt("-max_length")) {
		arg.setLt("-chunk_size",arg.ltOpt("-max_length"));
	}
	if (arg.ltOpt("-chunk_size") < 2*arg.ltOpt("-max_length") && 
	    arg.ltOpt("-chunk_size") < len_2) {
		std::cerr << "Using a chunk_size of less than two times -max_length ";
		std::cerr << "make the program very slow. Chunk_size is therefore ";
		std::cerr << "adjusted to: ";
		if (2*arg.ltOpt("-max_length") < len_2) {
			std::cerr << int(2*arg.ltOpt("-max_length")) << std::endl;
			arg.setLt("-chunk_size",2*arg.ltOpt("-max_length"));
		}
		else {
			std::cerr << int(len_2) << std::endl;
			arg.setLt("-chunk_size",len_2);
		}

	}
	if (arg.ltOpt("-chunk_size") > len_2) { arg.setLt("-chunk_size",len_2);}
	if (arg.boolOpt("-global")) {arg.setLt("-max_length", real_max_length);}
	if (arg.boolOpt("-global") && 
	    (arg.ltOpt("-max_diff") + real_min_length < real_max_length)) {

			lengthType new_delta = lengthType(1.1*(real_max_length - real_min_length));
			if (new_delta > real_max_length) {new_delta = real_max_length;}
			std::cerr << "The length difference between the two sequences ";
			std::cerr << seq_1->getName() << " and " << seq_2->getName() << " is ";
			std::cerr << (len_1 - len_2);
			std::cerr << " nucleotides which is more than the maximum length";
			std::cerr << "difference -max_diff " << arg.ltOpt("-max_diff");
			std::cerr << ". It is therefore not possible to use option -global.";
			std::cerr << " The -max_diff is therefore increased to: " << new_delta;
			std::cerr << std::endl;

			arg.setLt("-max_diff", new_delta);
	}


	// Since it makes no sense to run foldalign with option -no_backtrack without
	// also using option -plot_score option -plot_score is set when -no_backtrack
	// is set
	if (arg.boolOpt("-no_backtrack")) {
		if (!arg.boolOpt("-plot_score")) {
			arg.setBool("-plot_score", true);
		}
	}

	// This makes sure that the loop table length are big enough.
	score.checkSize(arg.ltOpt("-max_length")+1);

	arg.setPt("lenSeq1", len_1);
	arg.setPt("lenSeq2", len_2);
			
        bool global = arg.boolOpt("-global");
//        if (defaultScorematrix) {
            score.setGCparameters(seq_1->getGCcontent(), seq_2->getGCcontent(),
				global, seq_1->getLength(), seq_2->getLength());
//        }
	if (global) {
            // This is necessary to get the right length difference in the global pruning
            score.updatePruning(seq_1->getLength(), seq_2->getLength(), arg.ltOpt("-max_diff"));
        }

	out = new output(arg, seq_1, seq_2, score);

	std::string error;
	
	std::string bt_cons_file = "<none>"; //arg.stringOpt("-backtrack_seed_constraints");
	std::string cons_file = "<none>"; //arg.stringOpt("-seed_constraints");
	
	try {
		bt_cons = setConstraints(bt_cons_file, arg);
		cons = setConstraints(cons_file, arg);
	}
	catch ( exception ) {
		cleanUp(out, cons, bt_cons);
		throw;
	}
	catch ( ... ) {
		cleanUp(out, cons, bt_cons);
		error = "Could not read constraints files";
		throw exception(error, false);
	}
	
	if ( !arg.boolOpt("-global") ) {

		try {

			local<cell, longCell, longCell, false, false, false> 
			     (i, j, k, l, seq_1, seq_2, arg, score, out, cons, bt_cons);

		}
		catch ( exception& exc ) {
			cleanUp(out, cons, bt_cons);
			throw;
		}
		catch ( ... ) {
			cleanUp(out, cons, bt_cons);
			error = "Could not align sequences.";
			throw exception(error, false);
		}
	}
	else {
		try {
			backtrack(big_neg, big_neg, i, j, k, l, seq_1, seq_2, arg, score, out, cons, bt_cons);
		}
		catch ( exception ) {
			cleanUp(out, cons, bt_cons);
			throw;
		}
		catch ( ... ) {
			error = "Could not globally align the sequences.";
			out->saveOutput(arg.boolOpt("-plot_score"));
			cleanUp(out, cons, bt_cons);
			throw exception(error, false);
		}
	}

	cleanUp(out, cons, bt_cons);

	return 0;
}

arguments* setupArguments(int argc, char* argv[]) {
  
  arguments* arg = new arguments(argc, argv);

	// Store information in arguments class
	// Version number defined above
	arg->addString("version", version);
	// true if the two sequences has switch places. False otherwise
	arg->addBool("switch", false);
	// The length of the sequences
	arg->addPt("lenSeq1", 0);
	arg->addPt("lenSeq2", 0);
	// true during the realignment runs
	arg->addBool("realigning", false);
	// true during the mblrealignment run
	arg->addBool("mblrealign", false);
	// True if the folding can start from a zero length hp stem
	// false if the start coordinate is given and no new seed alignment is needed
	arg->addBool("hpstart", true);

	if (arg->stringOpt("-memory_roof").compare("-1")) {

		// There is a roof. Find the size.
		std::string line = arg->stringOpt("-memory_roof");
		const char end = line[line.length()-1];

		long factor = 1;
		switch (end) {
			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
				break;
			case 'm':
			case 'M':
				factor = 1;
				line = line.substr(0, line.length() -1);
				break;
			case 'g':
			case 'G':
				factor = 1000;
				line = line.substr(0, line.length() -1);
				break;
			case 't':
			case 'T':
				factor = 1000000;
				line = line.substr(0, line.length() -1);
				break;
			default:
				std::cerr << "Could not process the -memory_roof option." << std::endl;
				exit(1);
		}

		positionType size = atoi(line.c_str()); // * factor;
		size *= factor;
		positionType min_size = 0; //20000000;
		if (size < min_size) {
			std::cerr << "The minimum -memory_roof is 20M. The roof has been raised to the minimum" << std::endl;
			size = min_size;
		}
		size -= min_size/2;
		arg->addPt("memory_roof", size);
	}
	else {
		arg->addPt("memory_roof", -1);
	}

  return arg;
}

void handleHelpAndVersionOptions(arguments* arg) {


	// Handle -help and -h options
	if (arg->boolOpt("-help") || arg->boolOpt("-h") || arg->numberOfArguments() > 2 ) {
		std::cout << program_name << " " << version << std::endl;
		std::cout << description << std::endl;
		std::cout << "Usage:\n" << program_name << " [<options>] <file_1> [<file_2>]\n";
		std::cout << program_name << " makes pairwise alignments between two sequence.\n";
		std::cout << "The options are:\n";
		arg->printOptions();
		delete arg;
		exit(0);
	}

	// Handle -version option
	if ( arg->boolOpt("-version") ) {
		std::cout << program_name << " " << version << std::endl;
		std::cout << description << std::endl;
		delete arg;
		exit(0);
	}
  
}

scorematrix* setupScorematrix(bool& defaultScorematrix, arguments* arg) {

  // Handle the -score_matrix option
	std::string score_name = arg->stringOpt("-score_matrix");
	const bool global = arg->boolOpt("-global");

	scorematrix* score = 0;

	if (!score_name.compare("<default>")) {

	  defaultScorematrix = true;
	  score = new scorematrix(arg->ltOpt("-max_diff"), global, arg->boolOpt("-no_pruning"));
	}
	else {

		defaultScorematrix = false;
		// Read the score matrix from a file
		score = new scorematrix(score_name, arg->ltOpt("-max_diff"), global, arg->boolOpt("-no_pruning"));
	}

	return score;
}

int main(int argc, char* argv[]) {

	// Don't search for all the answers at once.
	// A path is formed by laying one stone at a time.
	// One person saw the third man that night.
	// Three have seen him, yes, but not his body.
	// Only one, known to you, ready now to talk.
	// One more thing, you forgot something. 

	arguments* arg = 0;
	scorematrix* score = 0;
	seqs* input = 0;
	bool defaultScorematrix;
	try {
	  // Setup the posible options
	  arg = setupArguments(argc, argv);

	  handleHelpAndVersionOptions(arg);
	
	  score = setupScorematrix(defaultScorematrix, arg);

	  // Read the sequences
	  input = new seqs(*arg, *score);
	}
	catch ( exception& exc) {
	
	  std::cerr << exc.getMessage() << std::endl;
	  if (exc.getFatal()) {
	    deletePtr(arg);
	    deletePtr(score);
	    deletePtr(input);
	    return 1;
	  }
	}
	catch ( ... ) {
	
	  std::cerr << program_name << ": Unknown initialization error. Aborting." << std::endl;
	  deletePtr(arg);
	  deletePtr(score);
	  deletePtr(input);
	  return 1;
	  
	}
	
	bool more_seq = true;
	while (more_seq) {
		sequence* seq1 = 0;
		sequence* seq2 = 0;
		try {
			more_seq = input->get_next_pair(seq1, seq2, *arg, *score);

			if (seq1 == 0) {
				std::cerr << "There is no first sequence." << std::endl;
				continue;
			}
			else if (seq2 == 0) {
				std::cerr << "There is no second sequence. Check that the input file format and the -format option are the same." << std::endl;
				continue;
			}

			align(seq1, seq2, *arg, *score);
		}
		catch ( std::bad_alloc& bad ) {
			std::cerr << "Memory problem:" << std::endl;
			std::cerr << bad.what() << std::endl;
		}
		catch ( exception& exc ) {
			if (exc.getFatal()) {
				std::cerr << "The fatal error:" << std::endl;
				std::cerr << exc.getMessage() << std::endl;
				std::cerr << "was incountered while aligning the sequences:" << std::endl;
				std::cerr << seq1->getName() << std::endl;
				std::cerr << seq2->getName() << std::endl;
				delete arg; delete score; delete input;
				return 1;
			}
			else {
				std::cerr << "The error:" << std::endl;
				std::cerr << exc.getMessage() << std::endl;
				std::cerr << "was encountered while aligning the sequences:" << std::endl;
				std::cerr << seq1->getName() << std::endl;
				std::cerr << seq2->getName() << std::endl;
				if (more_seq) {
					std::cerr << program_name << " will attempt to continue with any remaining alignments." << std::endl;
				}
			}
		}
		catch ( ... ) {
			std::cerr << "An error was encountered while aligning sequences: " << std::flush;
			std::cerr << seq1->getName() << " and " << seq2->getName() << std::endl;
			std::cerr << program_name << " will attempt to continue with any remaining alignments." << std::endl;
		}
	}

	delete arg;
	delete score;
	delete input;
	
	return 0;

	// One day my log will have something to say about this.
	// My log saw something that night.

}
