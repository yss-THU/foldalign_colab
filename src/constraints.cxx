#ifndef CONSTRAINTS
#define CONSTRAINTS

#include <map>
#include <vector>

#include "foldalign.hxx"
#include "readConstraints.cxx"
#include "jl.cxx"
#include "helper.cxx"

class constraints {

public:
	constraints(std::string fname, bool swap) : minScore(-1*big_neg) {

		readconstraints< jl >(fname, cons, swap, minScore, maxLengths);
		reset();
		rvals = 0;
	}

	inline scoreType getMinScore() const {
		return minScore;
	}

	inline void reset() {
		it_i = cons.end();
		--it_i;
	}
	
	//! \brief Alignments with a length shorter than the minExpandLength
	//! are not expanded into new alignments. They can however be part of
	//! multiloops and there by become part of longer alignments
	inline lengthType getSeedMinExpandLength(lengthType alignLength);
		
	
	inline void get(const positionType i, const positionType k, jl*& retv, 
				positionType& size);

	inline bool existsEnclosedConstraints(const positionType i,
					const positionType k,
					const lengthType Wi,
					const lengthType Wk);

	~constraints() {
		if (rvals != 0) {
			delete[] rvals;
		}
	}
	
	inline void printConstraints();
	
private:

	typedef std::multimap< positionType, jl > mmap;
	typedef std::map< positionType, mmap > map;

	map cons;
	
	map::iterator it_i;
	
	jl* rvals;
	
	scoreType minScore;

	std::vector<lengthType> maxLengths;

	constraints(const constraints&);
	constraints& operator=(const constraints&);
};

inline lengthType constraints::getSeedMinExpandLength(lengthType alignLength) {
	
	if (maxLengths.size() < 2 || maxLengths[1] >= alignLength) {
		return 0;
	}
	
	if (maxLengths.size() == 2) {
		return maxLengths[0];
	}
		
	for(unsigned int i=2; i < maxLengths.size(); i++) {
		if (maxLengths[i] >= alignLength) {
			return maxLengths[i-2];
		}
	}
		
	return maxLengths[maxLengths.size()-2];
}


inline void constraints::get(const positionType i, const positionType k, 
							 jl*& retv, positionType& size) {
			
	size = 0;
	
	if (retv == rvals && retv != 0) {
		delete[] retv;
		retv = 0;
		rvals = 0;
	}
	else {
		if (retv != 0)  {
			delete[] retv;
			retv = 0;
		}
		if (rvals != 0) {
			delete[] rvals;
			rvals = 0;
		}
	}
	
		
	positionType current_i = it_i->first;

	// i is larger than the current i, this should only happen when a new
	// k-section is started. In this case use the reset function.
	if (current_i < i) {
		return;
	}
		
	// If current_i == i then everything is ok.
	// It is already known the it_i->fist is smaller than i so
	if (current_i > i) {
		if (it_i != cons.begin()) {
			// There exists a smaller it_i->first
			--it_i;
				
			current_i = it_i->first;
			if (current_i != i) {
				// Taking one step back did not locate the right it_i.
				// Use the find() function to locate the right it_i
				map::iterator tmp = it_i;
					
				it_i = cons.find(i);
				if (it_i == cons.end()) {
					// There is no i constraint available
					// Keep the old it_i to make sure it_i is valid next time
					it_i = tmp;
//					if (it_i != cons.begin() && it_i->first > i) {it_i--;}
					return;
				}
				current_i = it_i->first;
			}
		}
		else {
			// it_i is already as small as possible. ie there is no constraint
			// available
			return;
		}
	}

	size = it_i->second.count(k);
		
	// size == 0 if there is no k for this i
	if (size <= 0) {return;}

	rvals = new jl[size];
	
	const std::pair<mmap::iterator, mmap::iterator> range = 
												it_i->second.equal_range(k);
	const mmap::iterator start = range.first;
	const mmap::iterator end = range.second;
	
	positionType pos = 0;
	for(mmap::iterator it = start; it != end; ++it) {
	
		rvals[pos].j = it->second.j;
		rvals[pos].l = it->second.l;
		rvals[pos].similarityScore = it->second.similarityScore;
		rvals[pos].energyScore = it->second.energyScore;
		rvals[pos].state = it->second.state;
	
		pos++;
	}

	retv = rvals;
}

inline bool constraints::existsEnclosedConstraints(const positionType i,
					const positionType k,
					const lengthType Wi,
					const lengthType Wk) {
	
	for (map::iterator itI = cons.lower_bound(i);
			itI != cons.upper_bound(i+Wi); ++itI) {

		for (mmap::iterator itK = itI->second.lower_bound(k);
				itK != itI->second.upper_bound(k+Wk); ++itK) {

			if(itK->second.j <= i+Wi && itK->second.l <= k+Wk) {
				return true;
			}
		}
	}

	return false;
}

inline void constraints::printConstraints() {
	
	for(map::iterator itI = cons.begin(); itI != cons.end(); ++itI) {
		
		for(mmap::iterator itK = itI->second.begin();
				itK != itI->second.end(); ++itK) {
			
			const std::pair<mmap::iterator, mmap::iterator> range = 
									itI->second.equal_range(itK->first);
				
			for(mmap::iterator itjl = range.first; itjl != range.second;
						++itjl) {
				std::cout << itI->first << " " << itjl->second.j << " ";
				std::cout << itK->first << " " << itjl->second.l;
				std::cout << std::endl;
			}
			std::cout << "============" << std::endl;
		}
	}
}

#endif /* CONSTRAINTS */
