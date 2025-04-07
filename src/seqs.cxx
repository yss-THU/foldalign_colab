#ifndef SEQS
#define SEQS

#include "foldalign.hxx"
#include "sequence.cxx"
#include "sequenceArg.cxx"
#include "helper.cxx"
#include "scorematrix.cxx"
#include "arguments.cxx"
#include "readfile.cxx"
#include "exception.cxx"

#include <string>
#include <fstream>
#include <iostream>

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

class seqs {
public:

  inline seqs(arguments& argum, scorematrix& sm);
	
  inline bool get_next_pair(sequence*& seq1, sequence*& seq2, arguments& argum, scorematrix& sm);
	
  inline ~seqs();

private:
  arguments& argu;
  scorematrix& scores;              // Passed on to the sequence objects
  seqs();
  sequence** seq1;
  sequence** seq2;
  sequenceArg** seqA1;
  sequenceArg** seqA2;
  int num_seq1;
  int num_seq2;
  int curr_seq1;
  int curr_seq2;
  std::string fileType;
  int als;
  template<class type> void delArray(type**& seq, int num);

  void readfasta(sequence**& seq, int& num_seq, std::string file,
		 scorematrix& score, int group, int max_seq);
  void header(const std::string& line, std::string& name, std::string& comment);

  void readtab(sequence**& seq, int& num_seq, std::string file,
	       scorematrix& score, int group, int max_seq);

  void readpair(sequenceArg**& seq1, sequenceArg**& seq2, int& num_seq,
		std::string file, scorematrix& score, arguments arg,
		int group, int max_seq);
  void getNamesAndSeqs(std::string& name, std::string& sequen,
		       std::string& line, int& start);
  void parseArg(std::string line, int pos, arguments local_Arg,
		std::string& name1, std::string& name2, std::string& sequen1,
		std::string& sequen2, std::string file, int num_seq, int group);

  readfile* openfile(std::string name);
};

inline seqs::seqs(arguments& argum, scorematrix& sm) : argu(argum), scores(sm), fileType(argu.stringOpt("-format")) {
  
  // Align self? If so index 2 starts 0 positions from index 1
  // Otherwise index 2 starts 1 positions from index 1
  if (argu.boolOpt("-align_self")) {als = 0;}
  else {als = 1;}
	
  // An old number. Not really necessary anymore but needed by some functions.
  int max_number_seq = 1000;

  if (argu.numberOfArguments() == 2) {
    if (!fileType.compare("commandline")) {
      std::string sek_1 = argu.argument(0);
      std::string sek_2 = argu.argument(1);
      std::string name_1 ="Sequence_1";
      std::string name_2 ="Sequence_2";

      num_seq1 = 1;
      num_seq2 = 1;
      seq1 = new sequence*[num_seq1];
      seq2 = new sequence*[num_seq2];
      seq1[0] = new sequence(name_1, sek_1, scores, 1);
      seq2[0] = new sequence(name_2, sek_2, scores, 2);
      curr_seq1 = 0;
      curr_seq2 = 0;
    }
    else if (!fileType.compare("fasta")) {
      readfasta(seq1, num_seq1, argu.argument(0), scores, 1, max_number_seq);
      readfasta(seq2, num_seq2, argu.argument(1), scores, 2, max_number_seq);
      curr_seq1 = 0;
      curr_seq2 = 0;
    }
    else if (!fileType.compare("tab")) {
      readtab(seq1, num_seq1, argu.argument(0), scores, 1, max_number_seq);
      readtab(seq2, num_seq2, argu.argument(1), scores, 2, max_number_seq);
      curr_seq1 = 0;
      curr_seq2 = 0;
    }
    else if (!fileType.compare("pair")) {
      std::string error = "The pair format can only be used with one file at the time";
      throw exception(error, true);
    }
    else {
      std::string error = "Unknown input format: " + fileType;
      throw exception(error, true);
    }
  }
  else if (argu.numberOfArguments() == 1) {
    if (!fileType.compare("commandline")) {
      std::string error = "Only one sequence specified. " + program_name + " needs two sequences";
      throw exception(error, true);
    }
    else if (!fileType.compare("fasta")) {
      readfasta(seq1, num_seq1, argu.argument(0), scores, 1, max_number_seq);
      num_seq2 = 0;
      curr_seq1 = 0;
      curr_seq2 = 0 + als;
    }
    else if (!fileType.compare("tab")) {
      readtab(seq1, num_seq1, argu.argument(0), scores, 1, max_number_seq);
      num_seq2 = 0;
      curr_seq1 = 0;
      curr_seq2 = 0 + als;
    }
    else if (!fileType.compare("pair")) {
      readpair(seqA1, seqA2, num_seq1, argu.argument(0), scores, argu, 1, max_number_seq);
      num_seq2 = 0;
      curr_seq1 = 0;
      curr_seq2 = 0;
    }
    else {
      std::string error = "Unknown input format: " + fileType;
      throw exception(error, true);
    }
  }
  else if (argu.numberOfArguments() == 0) {
    if (!fileType.compare("commandline")) {
      std::string error = "No sequences specified. " + program_name + " needs two sequences";
      throw exception(error, true);
    }
    else if (!fileType.compare("fasta")) {
      readfasta(seq1, num_seq1, "", scores, 1, max_number_seq);
      num_seq2 = 0;
      curr_seq1 = 0;
      curr_seq2 = 0 + als;
    }
    else if (!fileType.compare("tab")) {
      readtab(seq1, num_seq1, "", scores, 1, max_number_seq);
      num_seq2 = 0;
      curr_seq1 = 0;
      curr_seq2 = 0 + als;
    }
    else if (!fileType.compare("pair")) {
      readpair(seqA1, seqA2, num_seq1, "", scores, argu, 1, max_number_seq);
      num_seq2 = 0;
      curr_seq1 = 0;
      curr_seq2 = 0;
    }
    else {
      std::string error = "Unknown input format: " + fileType;
      throw exception(error, true);
    }
  }
  else {
    // This should never be reached
    std::string error = "To many arguments. Only two files allowed.";
    throw exception(error, true);
  }
}


inline bool seqs::get_next_pair(sequence*& sequen1, sequence*& sequen2,
				arguments& argum, scorematrix& sm) {
	
  if (!fileType.compare("pair")) {
    sequen1 = seqA1[curr_seq1];
    sequen2 = seqA2[curr_seq1];  // curr_seq1 on purpose
    argum   = seqA1[curr_seq1]->getArguments();
    sm      = seqA1[curr_seq1]->getScorematrix();
    curr_seq1++;
    if (curr_seq1 >= num_seq1) {return false;}
    return true;
  }
  else if (num_seq2 == 0) {
    // This is the one file case
    sequen1 = seq1[curr_seq1];
    sequen2 = seq1[curr_seq2];
    argum   = argu;
    sm      = scores;
    curr_seq2++;
    if (curr_seq2 >= num_seq1) {
      curr_seq1++;
      if (curr_seq1 < num_seq1 -1) {
	curr_seq2 = curr_seq1 + als;
	return true;
      }
      else if ((curr_seq1 < num_seq1) && (als == 0)) {
	curr_seq2 = curr_seq1;
	return true;
      }
      else { return false; }
    }
    return true;
  }
  else {
    // This is the two files case
    sequen1 = seq1[curr_seq1];
    sequen2 = seq2[curr_seq2];
    argum   = argu;
    sm      = scores;
    curr_seq2++;
    if (curr_seq2 >= num_seq2) {
      curr_seq2 = 0;
      curr_seq1++;
      if (curr_seq1 >= num_seq1) {
	return false;
      }
    }
    return true;
  }
}

inline seqs::~seqs() {

  if (!fileType.compare("pair")) {
    delArray(seqA1, num_seq1);
    delArray(seqA2, num_seq1); // num_seq1 on purpose
  }
  else {
    delArray(seq1, num_seq1);

    if (num_seq2 != 0) {delArray(seq2, num_seq2);}
  }
}

template<class type>
inline void seqs::delArray(type**& seq, int num) {
  for(int i=0; i<num; i++) {delete seq[i];}
  delete[] seq;
}


//==========================================================================
// These functions reads FASTA files

inline readfile* seqs::openfile(std::string name) {

  try {
    return new readfile(name);
  }
  catch ( exception ) {throw;}
  catch ( ... ) {
    std::string error_name = name;
    if (!error_name.compare("")) {error_name = "<stdin>";}
    std::string error = "Could not open the sequence file: " + error_name;
    throw exception(error, true);
  }
}

inline void seqs::readfasta(sequence**& seq, int& num_seq, std::string file,
			    scorematrix& score, int group, int max_seq) {
  // Initilizing
  num_seq = 0;
  seq = new sequence*[max_seq];
  std::string line;
  std::string name;
  std::string comment;
  std::string sequen = "";

  // Opening and checking the file
  readfile* fil = openfile(file);

  fil->get_line(line);
  std::string filename = fil->name();

  if (line.substr(0,1).compare(">")) {
    std::string error = "File " + filename + " is not in fasta format";
    delete fil;
    throw exception(error, true);
  }

  // Process the first header
  header(line, name, comment);

  while (fil->get_line(line)) {
    if (!line.substr(0,1).compare(">")) {
      // check sequence
      if (sequen.length() == 0) {
				std::cerr << "Warning sequence " << name << " is empty" << std::endl;
      }
      // Check for overflow
      if (num_seq >= max_seq) {
				helper::expandArray(seq, max_seq, 2*max_seq);
				max_seq = 2*max_seq;
      }
      // Store the sequence
      seq[num_seq] = new sequence(name, sequen, score, group, comment, filename);
      // and get ready for the next sequence
      sequen = "";
      header(line, name, comment);
      num_seq++;
    }
    else {
      // Store the line
      sequen+=line;
    }
  }

  // Check for overflow
  if (num_seq >= max_seq) {
    helper::expandArray(seq, max_seq, max_seq+1);
    max_seq = max_seq+1;
  }

  // Store the last sequence
  seq[num_seq] = new sequence(name, sequen, score, group, comment, filename);
  num_seq++;

  // Close the file
  delete fil;
}

inline void seqs::header(const std::string& line, std::string& name,
			 std::string& comment) {

  // Get the name
  int pos = helper::nextSpace(line, 0);
  if (pos > 0) {name = line.substr(1, pos-1);}
  else {name = line.substr(1, pos);}

  // If there is a comment store it
  comment = "";
  if (pos > 0) {
    pos = helper::nextNotSpace(line, pos);
    if (pos >= 0) {
      comment = line.substr((pos));
    }
  }
}

// ==========================================================================
// These functions reads tab formated files

inline void seqs::readtab(sequence**& seq, int& num_seq, std::string file,
			  scorematrix& score, int group, int max_seq) {

  // Initilizing
  num_seq = 0;
  seq = new sequence*[max_seq];
  std::string line;
  std::string name;
  std::string comment = "";
  std::string sequen = "";

  // Opening and checking the file
  readfile* fil = openfile(file);
  std::string filename = fil->name();

  int pos=0;
  int last;
  int old_last;
  while (fil->get_line(line)) {
    comment = "";
    helper::findWhiteSpace(line, pos, last);
    if (pos == -1) {
      std::string error = "No whitespace detected. Sequence is not in tab format aborting";
      delete fil;
      throw exception(error, true);
    }
    name = line.substr(0,pos);
    old_last = last;
    helper::findWhiteSpace(line, pos, last, old_last);
    if (pos >= 0) {
      sequen = line.substr(old_last,(pos-old_last));
      if (sequen.length() == 0) {std::cerr << "Warning sequence " << name << " is empty" << std::endl;}
      helper::findWhiteSpace(line, pos, last, last);
      if (pos >= 0) {
	comment = line.substr(last);
      }
    }
    else {
      sequen = line.substr(old_last);
    }
    // Check for overflow
    if (num_seq >= max_seq) {
      helper::expandArray(seq, max_seq, 2*max_seq);
      max_seq = 2*max_seq;
    }
    seq[num_seq] = new sequence(name, sequen, score, group, comment, filename);
    num_seq++;
  }
  // Close the file
  delete fil;
}

// ==========================================================================
// These functions reads pair formated files



inline void seqs::readpair(sequenceArg**& seq1, sequenceArg**& seq2,
			   int& num_seq, std::string file, scorematrix& score,
			   arguments arg, int group, int max_seq) {

  // Initilizing
  num_seq = 0;
  seq1 = new sequenceArg*[max_seq];
  seq2 = new sequenceArg*[max_seq];
  std::string line;
  std::string name1;
  std::string name2;
  std::string sequen1 = "";
  std::string sequen2 = "";
  std::string error = "Unknown error found during reading of pair file";

  // Opening and checking the file
  readfile* fil = openfile(file);
  std::string filename = fil->name();

  // Keep reading while the file is ok
  while (fil->get_line(line)) {
    int start = 0;

    getNamesAndSeqs(name1, sequen1, line, start);

    if (start < 0) {
      error = "Unabel to read pair file " + filename;
      delete fil;
      throw exception(error, true);
    }

    getNamesAndSeqs(name2, sequen2, line, start);

    // Check for overflow and expand the array if necessary.
    if (num_seq >= max_seq) {
      helper::expandArray(seq1, max_seq, 2*max_seq);
      helper::expandArray(seq2, max_seq, 2*max_seq);
      max_seq = 2*max_seq;
    }

    if (start >= 0) {
      parseArg(line, start, arg, name1, name2, sequen1, sequen2, filename,
	       num_seq, group);
    }
    else {
      seq1[num_seq] = new sequenceArg(name1, sequen1, score, arg, group, "",
				      filename);
      seq2[num_seq] = new sequenceArg(name2, sequen2, score, arg, group, "",
				      filename);
    }
    num_seq++;
  }
  // Close the file
  delete fil;
}

inline void seqs::getNamesAndSeqs(std::string& name, std::string& sequen,
				  std::string& line, int& start) {

  int pos;
  int last;
  std::string error = "Unknown error?";
  helper::findWhiteSpace(line, pos, last, start);
  if (pos == -1) {
    error = "No whitespace detected. Sequence is not in pair format";
    throw exception(error, true);
  }

  name = line.substr(start,(pos-start));

  start = last;
  helper::findWhiteSpace(line, pos, last, start);
  if (pos >= 0) {
    sequen = line.substr(start,(pos-start));
    if (sequen.length() == 0) {
      std::cerr << "Warning sequence " << name << " is empty" << std::endl;
    }
  }
  else {
    sequen = line.substr(start);
  }

  start = last;
}


inline void seqs::parseArg(std::string line, int start, arguments local_Arg,
			   std::string& name1, std::string& name2,
			   std::string& sequen1, std::string& sequen2,
			   std::string file, int num_seq, int group) {
  int last;
  int pos;
  std::string opt;
  std::string value;
  while (start >= 0) {
    // Get the option
    helper::findWhiteSpace(line, pos, last, start);
    if (pos >= 0) {opt = line.substr(start, (pos-start));}
    else {opt = line.substr(start);}

    // Get the value
    start = last;
    helper::findWhiteSpace(line, pos, last, start);
    if (pos < 0) {value = line.substr(start);}
    else {value = line.substr(start, (pos-start));}

    if (local_Arg.existStOpt(opt)) {
      scoreType val = atoi(value.c_str());
      local_Arg.setSt(opt, val);
    }
    else if (local_Arg.existPtOpt(opt)) {
      positionType val = atoi(value.c_str());
      local_Arg.setPt(opt, val);
    }
    else if (local_Arg.existLtOpt(opt)) {
      lengthType val = atoi(value.c_str());
      local_Arg.setLt(opt, val);
    }
    else if (local_Arg.existStringOpt(opt)) {
      if(!opt.compare("-ID")) {
				value = helper::findNextString(line, start, last);
      }

      local_Arg.setString(opt, value);
    }
    else if (local_Arg.existBoolOpt(opt)) {
      if (!value.compare("true")) {local_Arg.setBool(opt, true);}
      else {local_Arg.setBool(opt, false);}
    }
    else {
      std::string error = "Unknown option: " + opt + " in pair file " + file;
      throw exception(error, true);
    }
    //    std::cout << opt << " : " << value <<  " "  << start << " " << pos << " " << last << std::endl;
    start = last;
  }

  if (argu.stringOpt("-score_matrix").compare(local_Arg.stringOpt("-score_matrix")) ||
      argu.boolOpt("-global") != local_Arg.boolOpt("-global")) {

    // There is a new score matrix
    std::string score_name = local_Arg.stringOpt("-score_matrix");
    const bool global = local_Arg.boolOpt("-global");
    try {
      if (score_name.compare("<default>")) {
				scorematrix local_score(score_name, global);
				seqA1[num_seq] = new sequenceArg(name1, sequen1, local_score, local_Arg, group, "",
																					file);
				seqA2[num_seq] = new sequenceArg(name2, sequen2, local_score, local_Arg,
																				 group, "", file);
      }
      else {
				scorematrix local_score(global);
				seqA1[num_seq] = new sequenceArg(name1, sequen1, local_score, local_Arg,
																					 group, "", file);
				seqA2[num_seq] = new sequenceArg(name2, sequen2, local_score, local_Arg,
																					 group, "", file);
      }	
      return;
    }
    catch ( exception& exc ) {
	
      if ( exc.getFatal() ) {
				std::string error = "Fatal error while reading the score_matrix: ";
				error += score_name;
				error += "\nError:\n";
				error += exc.getMessage();
		
				throw exception(error, true);
      }
		
      std::string error = "Trouble reading the score_matrix for " + name1;
      error += " and " + name2 + ". The score matrix: " + scores.getName();
      error += " will be used instead. The error was: " + exc.getMessage();
			
      std::cerr << error << std::endl;
	
    }
    catch ( ... ) {
		
      std::cerr << "Unknown error found while reading score matrix: " + score_name;
      std::cerr << ". Sequences " + name1 + " and " + name2 + " will be aligned ";
      std::cerr << "using the score matrix: " + scores.getName() << std::endl;
			
    }
  }
  seqA1[num_seq] = new sequenceArg(name1, sequen1, scores, local_Arg, group, "",
				   file);
  seqA2[num_seq] = new sequenceArg(name2, sequen2, scores, local_Arg, group, "",
				   file);
}

#endif /* SEQS */
