#ifndef READCONSTRAINTS
#define READCONSTRAINTS

#include<map>
#include <vector>
#include <algorithm>

#include "foldalign.hxx"
#include "readfile.cxx"
#include "helper.cxx"
#include "combineSimilarityEnergy.cxx"

template< class jl >
class readconstraints {

public:

	typedef std::multimap< positionType, jl > mmap_k;
	typedef std::map< positionType, mmap_k > map_i;


	readconstraints(std::string& filename, map_i& cons, bool swap, 
		scoreType& minScore, 
		std::vector<lengthType>& maxLengths) {
	
		std::string line;

		// Opening the new file
		readfile* file;
		try {
			file = new readfile(filename);
		}
		catch ( ... ) {
			std::string error = "It is not possible to open the constraint file: ";
			error += filename;
			throw exception(error, false);
		}

		minScore = -1*big_neg;

		while (file->get_line(line)) {

			if (!line.substr(0,11).compare("MAX_LENGTH:")) {
				int start = 12;
				maxLengths.push_back(helper::getValue(start, line));
				continue;
			}
			
			jl con;

			int prev = 0;
			
			positionType i = helper::getValue(prev, line);
			con.j =          helper::getValue(prev, line);
			positionType k = helper::getValue(prev, line);
			con.l =          helper::getValue(prev, line);
			con.similarityScore =      helper::getValue(prev, line);
			con.energyScore =      helper::getValue(prev, line);
			con.state =      helper::getValue(prev, line);

			if (i > con.j || k > con.l) {
				lengthType Wi = con.j - i;
				lengthType Wk = con.l - k;
				helper::pcerr("Found seed on the wrong strand?", i, k, Wi, Wk);
				std::cerr << "Ignoring this seed" << std::endl;
				continue;
			}

			if (swap) {
				std::swap(i,k);
				std::swap(con.j, con.l);
			}

			mmap_k& consK = cons[i];
			consK.insert(std::pair< positionType, jl >(k, con));
			
			scoreType score = combineSimilarityEnergy(con.similarityScore,
													  con.energyScore);

			if (score < minScore) {
				minScore = score;
			}
		}
	
		delete file;
		
		std::sort(maxLengths.begin(), maxLengths.end());
	}

};

#endif /* READCONSTRAINTS */
