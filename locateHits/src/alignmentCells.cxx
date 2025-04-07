#ifndef ALIGNMENTCELLS
#define ALIGNMENTCELLS

#include "nohit.hxx"

/*! \brief Stores the coordinates of an alignment

*/
class alignmentcells {

public:
	
	alignmentcells(posType bi=0, posType bj=0, posType bk=0, posType bl=0) 
		: i(bi), j(bj), k(bk), l(bl) {};
	
	void get(posType& bi, posType& bj, posType& bk, posType& bl) {
		bi = i; bj = j; bk = k; bl = l;
	}
	
	void set(posType bi, posType bj, posType bk, posType bl) {
		i = bi; j = bj; k = bk; l = bl;
	}

private:

	posType i;
	posType j;
	posType k;
	posType l;
};

#endif /* ALIGNMENTCELLS */

