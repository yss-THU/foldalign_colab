#ifndef TESTPRUNETABLE
#define TESTPRUNETABLE

#include "test.cxx"
#include "../../src/pruneTable.cxx"
#include "../../src/foldalign.hxx"
#include "../../src/readfile.cxx"
#include <iostream>
#include <string>

class testPruneTable {
public:
	testPruneTable(int& passed, int& ran, int& expected, std::string& messages) : tester(test(0, "pruneTable.cxx test:")) {
	
		lengthType size = 30;
		scoreType linearPrune = 3;
		scoreType start = -10;
		
		pruneTable prune(start, linearPrune, size, false);
		testFullTable(prune, start, linearPrune, size);		

		linearPrune = 1;
		prune.setPruneCoefficient(linearPrune);
		testFullTable(prune, start, linearPrune, size);		

		size = 42;
		prune.reSizeTable(size);
		testFullTable(prune, start, linearPrune, size);		

		// Test assign
		lengthType newSize = 100;
		scoreType newLinearPrune = 5;
		scoreType newStart = 7;
		pruneTable newPrune(newStart, newLinearPrune, newSize, false);
		testFullTable(newPrune, newStart, newLinearPrune, newSize);		

		prune = newPrune;
		testFullTable(prune, newStart, newLinearPrune, newSize);		
		testFullTable(newPrune, newStart, newLinearPrune, newSize);		

		linearPrune = 3;
		prune.setPruneCoefficient(linearPrune);
		testFullTable(newPrune, newStart, newLinearPrune, newSize);		

		prune = prune;
		testFullTable(prune, newStart, linearPrune, newSize);		
		
		// Test copy
		pruneTable newestPrune = prune;
		testFullTable(newestPrune, newStart, linearPrune, newSize);		

		newLinearPrune = 7;
		newestPrune.setPruneCoefficient(newLinearPrune);
		testFullTable(newestPrune, newStart, newLinearPrune, newSize);	
		testFullTable(prune, newStart, linearPrune, newSize);		


		// Test reading from file
		std::string filename = "test/auxfiles/pruneTable.cxx.aux";
		readfile file(filename);
		readfile* p2file = &file;
		prune.buildTableFromFile(p2file);
		lengthType tableSize = prune.getLength();
		tester.equal(tableSize, lengthType(18), "Readfile table size");

		for(int i = 0; i < tableSize; i++){
		
			if (i < 2) {
				tester.equal(prune.get(i), scoreType(400), "readfile < 2");
			}
			else if (i < 10) {
				tester.equal(prune.get(i), scoreType(-300), "readfile < 10");
			}
			else if (i < 17) {
				tester.equal(prune.get(i), scoreType(-5), "readfile < 17");
			}
			else {
				tester.equal(prune.get(i), scoreType(20), "readfile > 17");
			
			}
		}

		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	};
private:

	test tester;

	void testFullTable(pruneTable& prune, const scoreType start, const scoreType linear,
			const lengthType size);
};

void testPruneTable::testFullTable(pruneTable& prune, const scoreType start, const scoreType linear,
			const lengthType size) {

	tester.equal(prune.get(0), start);
	for(int i=1; i < size; i++) {
		scoreType expect = start + i*linear;
		tester.equal(prune.get(i), expect );
	}
	tester.equal(prune.getLength(), lengthType(size+1));
}

#endif
