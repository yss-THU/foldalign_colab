#ifndef TESTCELL
#define TESTCELL

#include "test.cxx"
#include "../../src/cell.cxx"
#include "../../src/foldalign.hxx"
#include "../../src/mbllist.cxx"

class testCell {
public:
	testCell(int& passed, int& ran, int& expected, std::string& messages) {

		// Test setup
		int nTest = 17;
		test master(nTest, "cell.cxx test");
	
		// Cell setup
		scoreType simSc = 0;
		scoreType enSc = 1;
		stateType st = 2;
		lengthType len1 = 3;
		lengthType len2 = 4;
		lengthType len3 = 5;
		lengthType len4 = 6;
	
		cell tCell(simSc, enSc, st, len1, len2, len3, len4);
	
		// Test scores getters
		master.equal(tCell.getSimilarityScore(), simSc);
		master.equal(tCell.getEnergyScore(), enSc);

		// Test state getters
		master.equal(tCell.getState(), st);
		tCell.setState(10);
		master.equal(tCell.getState(), stateType(10));
	
		// Test one length getters
		master.equal(tCell.getLength1(), len1);
		master.equal(tCell.getLength2(), len2);
		master.equal(tCell.getLength3(), len3);
		master.equal(tCell.getLength4(), len4);

		// Test the multiple length getters and setters
		scoreType sS = 11;
		scoreType eS = 13;
		stateType state = 15;
		lengthType l1 = 17;
		lengthType l2 = 19;
		lengthType l3 = 21;
		lengthType l4 = 23;
		tCell.set(sS, eS, state, l1, l2, l3, l4);

		master.equal(tCell.getSimilarityScore(), sS);
		master.equal(tCell.getEnergyScore(), eS);
		master.equal(tCell.getState(), state);

		tCell.getLengths(len1, len2, len3, len4);
		master.equal(l1, len1);
		master.equal(l2, len2);
		master.equal(l3, len3);
		master.equal(l4, len4);

		// Test the mblList pointer getters and setters
		mbllist* null = 0;
		master.equal(tCell.getPointer(), null);
		mbllist* mblPointer = new mbllist;
		tCell.setPointer(mblPointer);
		master.equal(tCell.getPointer(), null);
		delete mblPointer;

		// output results
		master.printResult();

		passed += master.getTestsPassed();
		ran += master.getTestsRan();
		expected += master.getTestsExpected();
	}
private:
};




#endif
