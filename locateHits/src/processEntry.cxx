#ifndef PROCESSENTRY
#define PROCESSENTRY

#include "nohit.hxx"
#include "datalist.cxx"
#include "alignmentCells.cxx"
#include "nonoverlap.cxx"
#include "pValue.cxx"

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdlib.h>


class processEntry {

public:

	// Reads a foldalign alignment. Parses the header and the LS lines.
	// The LS lines are stored in a datalist object which is sendt to the
	// nonoverlap.cxx algorithm. nonoverlap returns another datalist containing
	// the non-overlapping alignments. These are then printed.
	processEntry(std::istream& file, const scoreType& score_cut = bad,
	             const bool border = false, const lengthType ma = -1, const bool estimateParameters = false);
	
private:

	// The informative part of the header lines starts at this position
	static const unsigned int parameter_start_pos = 22;

	// Parse the header of the foldalign output
	bool parseHead(std::istream& input_file, bool& read_ok);


	// Read the LS lines and store those with a score above cutoff.
	// Then call nonoverlap to get the list of non-overlapping alignments.
	// Finally call printAlign to output the non-overlapping alignments.
	void parseLS(std::istream& input_file, datalist<alignmentcells>& alignments, const scoreType score_cut,
					 const bool border);


	// Setup the p-value calculation
        pValue setupPvalues(datalist<alignmentcells>& hitlist, const scoreType scoreCutoff, const bool estimateParameters);

	// Print the first line with the run info
	void printHeaderLine(pValue& pval, const scoreType score_cut);

	// Print the non-overlapping alignments.
	void printAlign(datalist<alignmentcells>& hitlist, pValue& pval);

	// Skip the parts of the foldalign output which is not needed.
	void skipTheRest(std::istream& input_file);

	// Read a line from the file
	bool get_line(std::istream& file, std::string& line, const bool warn = true);

	// In case of unexpected end of file print an error and throw one too.
	void readError(const int error_number) {
		std::cerr << program_name << ": Unexpected end of file" << std::endl;
		throw error_number;
	}

	// Find the first "word" after position prev in line and translate it to int
	scoreType getValue(int& prev, std::string& line);

	// Find the first word after position prev in line
	std::string findName(int& prev, std::string& line);

	bool parseNumber(posType& value, const std::string& line, const std::string& target,
                    const unsigned int& position, const unsigned int& len) const;

	// No defaults constructors
	processEntry();

	processEntry(const processEntry&);
	
	processEntry& operator=(const processEntry&);

	// The width of the position numbering
	static const int width = 6;
	static const posType preDef = 3;

	// Precission of lambda and k
	static const unsigned int precission = 9;

	// The max number of alignments to be printed. Values < 0 means all alignments
	const lengthType max_alignments;

	std::string name_I;
	std::string name_K;
	
	posType lambda;
	
	posType i_start;
	posType j_end;
	posType k_start;
	posType l_end;
	float gc1;
	float gc2;

	std::string id;

};

inline processEntry::processEntry(std::istream& file, 
                                  const scoreType& score_cut,
				  const bool border,
				  const lengthType ma, const bool estimateParameters) : max_alignments(ma) {

	
	bool ls_head = false; // Becomes true for LS header false otherwise
	bool file_ok = true; // True if the file is alright
	while (!ls_head && file_ok) {
		ls_head = parseHead(file, file_ok);
	}

	if (file_ok) {
		datalist<alignmentcells> alignments;
	
		parseLS(file, alignments, score_cut, border);
	
		datalist<alignmentcells> hitlist;
		nonoverlap<alignmentcells> no(alignments, hitlist, border, lambda, max_alignments);
	
		pValue pval = setupPvalues(hitlist, score_cut, estimateParameters);

		printHeaderLine(pval, score_cut);
		printAlign(hitlist, pval);
	}
}

inline bool processEntry::parseHead(std::istream& input_file, bool& read_ok) {

//std::cout << "........................................" << std::endl;
	std::string line = "";

	lambda = 0;
	i_start = 1;
	k_start = 1;
	j_end = 0;
	l_end = 0;
	gc1 = 0;
	gc2 = 0;
	id = "";
	
	bool ls_head = false;
	
	// It is ok for this read to fail if this is the end the file
	read_ok = get_line(input_file, line);
	if (!read_ok) {return false;}

	while (line.substr(0,5).compare("; ---") && line.substr(0,5).compare("; ===")) {

		unsigned int len = line.length();
		unsigned int lo = 11;
		unsigned int co = 2;
//std::cout << line << " " << read_ok << " " << ls_head <<std::endl;

		if (len <= parameter_start_pos) {}
		else if ((len > parameter_start_pos) && 
		         (!line.substr(0, parameter_start_pos).compare("; ALIGNING            "))) {

			for(unsigned int pos = parameter_start_pos; pos < len; pos++) {
				
				if (!line.substr(pos, 9).compare(" against ")) {
					name_I = line.substr(parameter_start_pos, pos - parameter_start_pos);
					name_K = line.substr(pos + 9);
					break;
				}
			}
		}
		else if ((len > parameter_start_pos) && 
		         (!line.substr(0, parameter_start_pos).compare("; ALIGNMENT_ID        "))) {
			id = line.substr(parameter_start_pos).c_str();
		}
		else if ((len > parameter_start_pos) && 
				(!line.substr(0, parameter_start_pos).compare("; GC_CONTENT_SEQ_1    "))) {
			gc1 = atof(line.substr(parameter_start_pos).c_str());
		}
		else if ((len > parameter_start_pos) &&
				(!line.substr(0, parameter_start_pos).compare("; GC_CONTENT_SEQ_2    "))) {
			gc2 = atof(line.substr(parameter_start_pos).c_str());
		}
		
		parseNumber(lambda, line, std::string("; PARAMETER           max_length="), parameter_start_pos+lo, len);

		parseNumber(i_start, line,  std::string("; PARAMETER           i="), parameter_start_pos+co, len);
									
		parseNumber(k_start, line,  std::string("; PARAMETER           k="), parameter_start_pos+co, len);
									
		parseNumber(j_end, line,  std::string("; PARAMETER           j="), parameter_start_pos+co, len);
									
		parseNumber(l_end, line,  std::string("; PARAMETER           l="), parameter_start_pos+co, len);
									
		if (!line.compare("; TYPE                Foldalign_local_scores")) {
			ls_head = true;
		}
		else if (!line.compare("; TYPE                RNA")) {
			skipTheRest(input_file);
			return false;
		}
		
		read_ok = get_line(input_file, line);
		if (!read_ok) {readError(-1);}
	}

	return ls_head;
}
	
inline bool processEntry::parseNumber(posType& value, const std::string& line,
                        const std::string& target,
                        const unsigned int& position,
								const unsigned int& len) const {

	if (len > position && !line.substr(0, position).compare(target)) {
		value = atoi(line.substr(position).c_str());
		return true;
	}
	return false;
}


inline void processEntry::parseLS(std::istream& input_file, datalist<alignmentcells>& alignments,
                                  const scoreType score_cut, 
				  const bool border) {
	

	std::string line;
//std::cout << "Start reading" << std::endl;
	bool read_ok = get_line(input_file, line);
	if (!read_ok) {readError(-1);}

	const posType i_border = lambda;
	const posType k_border = lambda;
	const posType j_border = j_end - lambda;
	const posType l_border = l_end - lambda;

	if (border && (j_border <= i_border || l_border <= k_border)) {
		std::cerr << "Warning! The sequences " << name_I << " and " << name_K << " are too short to produce an alignment using the -lambda_border option." << std::endl;
		return;
	}

	try {
		while (line.substr(0,5).compare("; ***")) {

			if (line.substr(0, 2).compare("LS")) {

				if (line.substr(0,4).compare("; AS")) {
					std::cerr << program_name <<": Input file is not in the right format" << std::endl;
				}
			
				read_ok = get_line(input_file, line);
				if (!read_ok) {readError(-1);}
			
				continue;
			}
		
			int pre = preDef;
			posType i = getValue(pre, line);
			posType j = getValue(pre, line);
			posType k = getValue(pre, line);
			posType l = getValue(pre, line);
			alignmentcells align(i, j, k, l);
			scoreType score = getValue(pre, line);

			if (score >= score_cut) {
				// It is not necessary to check for i < j_border since i < j < j_border
				// and likewise for the other three "missing" cases
				if (!border || (i > i_border && k > k_border &&
				                j < j_border && l < l_border)) {
					alignments.add(score, align);
				}
			}
	
			read_ok = get_line(input_file, line);
			if (!read_ok) {readError(-1);}
		
		}
	}
	catch ( ... ) {
		std::cerr << program_name <<  ": An error happened while reading the LS lines for the alignment of " << name_I << " against " << name_K << std::endl;
	}
//std::cout << "Finished reading" << std::endl;
}

inline pValue processEntry::setupPvalues(datalist<alignmentcells>& hitlist, const scoreType scoreCutoff, const bool estimateParameters) {
	
	pValue pval((j_end-i_start+1), (l_end - k_start+1), gc1, gc2);

	if (estimateParameters) {
		pval.estimateParameters(hitlist, scoreCutoff);
	}

	return pval;
}

inline void processEntry::printHeaderLine(pValue& pval, const scoreType score_cut) {

	std::cout << "# New alignment: ";
	std::cout << name_I << " " << i_start << " " << j_end << " GC: " << gc1 << " ";
	std::cout << name_K << " " << k_start << " " << l_end << " GC: " << gc2 << " ";
	std::cout << " max_length: " << lambda;
	std::cout << " Lambda: " << std::setprecision(precission) << pval.getLambda();
	std::cout << " K: " << std::setprecision(precission) << pval.getK();
	std::cout << " Score cut: " << score_cut;
	std::cout << " ID: " << id  << std::endl;
}

inline void processEntry::printAlign(datalist<alignmentcells>& hitlist, 
		pValue& pval) {

	scoreType sa;
	alignmentcells a;

	posType i;
	posType j;
	posType k;
	posType l;
	
	hitlist.initCheckNext();
	
	hitlist.checkNext(sa, a);
	
	//pValue pval = pValue(len1, len2, gc1, gc2);
	
	long count = 1;
	while (sa > bad) {
		
		a.get(i, j, k, l);

		std::cout << name_I << " " << std::setw(width) << i;
		std::cout << " " << std::setw(width) << j << " ";
		std::cout << name_K << " " << std::setw(width) << k;
		std::cout << " " << std::setw(width) << l << " ";
		std::cout << std::setw(width) << sa;
		std::cout.precision(3);
		std::cout << "   " << std::fixed << pval.calcPvalue(sa);
		std::cout << std::setw(width) << count << std::endl;
		std::cout.precision(precission);
//		std::cout << i << ": " << m << " " << n << " " << k << " " << l <<  " Score:" << sa << std::endl;
		
		hitlist.checkNext(sa, a);

		if (max_alignments > -1 && count >= max_alignments) {break;}
		count++;
	}
	
//	std::cout << sa << std::endl;
//	std::cout << "============================================="<< std::endl;
}	



inline void processEntry::skipTheRest(std::istream& input_file) {

	std::string line;
	
	
	int line_count = 0;
	
	while (line_count < 2) {

		bool read_ok = get_line(input_file, line);
		if (!read_ok) {readError(-1);}
		
		if (!line.substr(0,5).compare("; ***")) {
			line_count++;
		}
	}
}

//==========================================================================
// General file line reader and processors

inline bool processEntry::get_line(std::istream& file, std::string& line,
                                   const bool warn) {
	std::getline(file,line);
	while (((line[0] == '#') || (line.length() == 0)) && (!file.fail())){
		std::getline(file,line);
	}

	if (file.fail()) {return false;}

	return true;
}

inline scoreType processEntry::getValue(int& prev, std::string& line) {
	return atoi(findName(prev, line).c_str());
}
	
inline std::string processEntry::findName(int& prev, std::string& line) {
	// Find the first word after position prev in line

	int pos;
	int end_pos = -1;
	int len = line.length();
	for(pos=prev; pos < len; pos++) {
		// Find the first non space or tab after the prev position.
		if ((line[pos] != ' ') && (line[pos] != '\t')) {
			// Now find the first space or tab after the previously found position.
			end_pos = pos+1;
			while ((line[end_pos] != ' ') && (line[end_pos] != '\t') && 
			       (end_pos < len)) {end_pos++;}

			// The next non-space -tab character can start no earlier than
			//end_pos+1
			prev = end_pos +1;
			
			// Return the word
			return line.substr(pos, (end_pos - pos));
		}
	}
//	if (end_pos == -1) {
		std::cerr << "Could not process the LS line.\n" << line << std::endl;
		throw -1;
//	}
//	prev = end_pos+1;
//	return line.substr(pos, (end_pos - pos));
	return "";
}


#endif /* PROCESSENTRY */
