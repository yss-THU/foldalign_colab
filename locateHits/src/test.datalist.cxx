#include "nohit.hxx"
#include "datalist.cxx"
#include "alignmentCells.cxx"
#include "nonoverlap.cxx"

#include <iostream>

/* The correct output should be:
11 31 32 33 34
11 41 42 43 44
11 51 52 53 54
7 61 62 63 64
7 71 72 73 74
5 81 82 83 84
3 1 2 3 4
3 11 12 13 14
3 21 22 23 24
-10
*/

void printCheck(datalist<alignmentcells>& testList) {

	scoreType sa;
	alignmentcells a;

	posType m;
	posType n;
	posType k;
	posType l;
	
	testList.initCheckNext();
	
	testList.checkNext(sa, a);
	
	int i = 1;
	while (sa > bad) {
		
		a.get(m, n, k, l);

		std::cout << i << ": " << m << " " << n << " " << k << " " << l <<  " Score:" << sa << std::endl;
		
		testList.checkNext(sa, a);

		i++;
	}
	
	std::cout << sa << std::endl;
	std::cout << "============================================="<< std::endl;
}	

int main() {
//	scoreType bad = -10;

	datalist<alignmentcells> testList;

	alignmentcells a(1,2,3,4);
	alignmentcells b(11,12,13,14);
	alignmentcells c(21,52,23,84);
	alignmentcells d(31,32,33,34);
	alignmentcells e(11,42,43,44);
	alignmentcells f(51,52,53,54);
	alignmentcells g(61,62,63,64);
	alignmentcells h(71,72,73,74);
	alignmentcells i(81,82,83,84);
	
	
	scoreType sa = 3;
	scoreType sb = 3;
	scoreType sc = 3;
	scoreType sd = 11;
	scoreType se = 11;
	scoreType sf = 11;
	scoreType sg = 7;
	scoreType sh = 7;
	scoreType si = 5;
	
	testList.add(sa, a);
	testList.add(sb, b);
	testList.add(sc, c);
	testList.add(sd, d);
	testList.add(se, e);
	testList.add(sf, f);
	testList.add(sg, g);
	testList.add(sh, h);
	testList.add(si, i);

	std::cout << "PrintCheck 1 testList" << std::endl;
	printCheck(testList);
	std::cout << "PrintCheck 2 testList" << std::endl;
	printCheck(testList);
	
	std::cout << "PrintCheck testList done" << std::endl;
	
	datalist<alignmentcells> testList2;
	std::cout << "Nonoverlap" << std::endl;
	nonoverlap<alignmentcells> no(testList, testList2);

	std::cout << "PrintCheck testList2" << std::endl;
	printCheck(testList2);

}
