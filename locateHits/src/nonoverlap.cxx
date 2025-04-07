#ifndef NONOVERLAP
#define NONOVERLAP

#include "nohit.hxx"
#include "alignmentCells.cxx"
#include "datalist.cxx"

#include<iostream>

template< class alignment >
class nonoverlap {

public:
	nonoverlap(datalist<alignment>& alignments, 
				  datalist<alignment>& hitlist, 
				  const bool border, const posType lambda, const lengthType ma);

private:

	nonoverlap();
	
	void checkCoords(scoreType& score, 
			const posType i,
			const posType j,
			const posType k,
			const posType l,
			alignment align,
			datalist<alignment>& hitlist,
			const bool border,
			const posType lambda,
                        lengthType& alignmentCount
		);

};

template< class alignment >
inline nonoverlap< alignment >::nonoverlap(datalist<alignment>& alignments,
                                           datalist<alignment>& hitlist,
					   const bool border,
					   const posType lambda,
					   const lengthType ma
					  ) {
	scoreType score;
	alignment align;
	posType i;
	posType j;
	posType k;
	posType l;
	
	alignments.getNext(score, align);

        lengthType alignmentCount = 0;
	while(score > bad) {
	
		align.get(i,j,k,l);

		checkCoords(score, i, j, k, l, align, hitlist, border, lambda, alignmentCount);
		
                if (ma > -1 && alignmentCount > ma) {
                    break;
                }
                alignments.getNext(score, align);

	}
}

template< class alignment >
inline void nonoverlap< alignment >::checkCoords(scoreType& score, 
					 const posType i,
					 const posType j,
					 const posType k,
					 const posType l,
					 alignment align,
					 datalist<alignment>& hitlist,
					 const bool border,
					 const posType lambda,
					 lengthType& alignmentCount
				 ) {

	scoreType s;
	alignment a;
	posType ci;
	posType cj;
	posType ck;
	posType cl;

	hitlist.initCheckNext();
	hitlist.checkNext(s, a);

	while(s > bad) {

		if (s < score) {
			std::cerr << "Program error. Warning: score s is larger than score" << std::endl;
		}
		
		a.get(ci,cj,ck,cl);
		
		if (border) {
			if (i - lambda < ci && ci < j + lambda &&
			    k - lambda < ck && ck < l + lambda) {

				return;
			}
		}
		else {
			if (i <= cj && j >= ci && k <= cl && l >= ck) {
				return;
			}
		}
			
		hitlist.checkNext(s,a);
	}
	
	hitlist.add(score, align);
        alignmentCount++;
	
}
#endif /* NONOVERLAP */
