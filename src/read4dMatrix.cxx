#ifndef READ4DMATRIX
#define READ4DMATIRX

#include "matrix4d.cxx"
#include "readfile.cxx"
#include "scorematrix.cxx"
#include "helper.cxx"


template<class data>
class read4dMatrix {
public:
	read4dMatrix() {};
	
	matrix4d<data>* read(const int& size, readfile*& file, int (scorematrix::*alfa)(char) const );
	
private:

	inline void parseLine(std::string line, int*& i_letters, int& len, int (scorematrix::*alfa)(char) const);
	int getLetterIndex(std::string line, int (scorematrix::*alfa)(char) const, int& prev);
	

};

template<class data>
matrix4d<data>* read4dMatrix<data>::read(const int& size, readfile*& file, int (scorematrix::*alfa)(char) const) {

	matrix4d<data>* matrix = new matrix4d<data>(size);

	std::string line;
	int len = size*size;

	int prev = 0;
	int pos;
	int* i_letters = new int[len];
	int* j_letters = new int[len];
	int k_letter;
	int l_letter;
	char letter;

	file->get_line(line);
	parseLine(line, i_letters, len);

	file->get_line(line);
	parseLine(line, j_letters, len);

	for(int i=0; i<len; i++) {
		file->get_line(line);
		prev=line.length();
		
		l_letter = getLetterIndex(line, alfa, prev);
		k_letter = getLetterIndex(line, alfa, prev);
		
		prev = 0;
		for(int j=0; j<len; j++) {
			matrix.set(i_letters[j], j_letters[j], k_letter, l_letter, helper::getValue(prev, line));
		}
	}
	delete[] i_letters;
	delete[] j_letters;	
}

template<class data>
inline void read4dMatrix<data>::parseLine(std::string line, int*& i_letters, int& len, int (scorematrix::*alfa)(char) const) {
	char letter;
	int pos=0;
	int length = line.length();
	for(int i=0; i<len; i++) {
		while ((line[pos] == ' ') || (line[pos] == '\t')) {pos++;}
		if (pos > length) {
			std::string error = "Could not 4d matrix (matrix4d.cxx). Somethings wrong with this line\n";
			error += line;
			throw exception(error, false);
		}
		letter = line[pos];
		i_letters[i] = alfa(letter);
		pos++;
	}
}

template<class data>
int read4dMatrix<data>::getLetterIndex(std::string line, int (scorematrix::*alfa)(char) const, int& prev) {
		int pos = line.find_last_not_of(" ", prev);
		char letter = line[pos];
		prev = pos -1;
		return alfa(letter);
}

#endif
