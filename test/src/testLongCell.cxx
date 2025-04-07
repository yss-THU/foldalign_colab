#ifndef TESTLONGCELL
#define TESTLONGCELL

#include "test.cxx"
#include "../../src/longCell.cxx"
#include "../../src/foldalign.hxx"
#include "../../src/mbllist.cxx"

class testLongCell {
public:
	testLongCell(int& passed, int& ran, int& expected, std::string& messages) {

		// Test setup
		int nTest = 12;
		test tester(nTest, "longCell.cxx test");
	
		// Cell setup
		scoreType simSc = 0;
		scoreType enSc = 1;
		stateType st = 2;
		mbllist* mblNull = 0;
		
		longCell defaultCell;
		tester.equal(defaultCell.getSimilarityScore(), big_neg);
		tester.equal(defaultCell.getEnergyScore(), big_neg);
		tester.equal(defaultCell.getState(), noState);
		tester.equal(defaultCell.getPointer(), mblNull);
		
		longCell tCell(simSc, enSc, st);
		tester.equal(tCell.getSimilarityScore(), simSc);
		tester.equal(tCell.getEnergyScore(), enSc);
		tester.equal(tCell.getState(), noState);
		tester.equal(tCell.getPointer(), mblNull);

		simSc = -10;
		enSc = 200;
		st = 5;
		mbllist* mbl = new mbllist();
		
		defaultCell.set(simSc, enSc, st, mbl);
		tester.equal(defaultCell.getSimilarityScore(), simSc);
		tester.equal(defaultCell.getEnergyScore(), enSc);
		tester.equal(defaultCell.getState(), noState);
		tester.equal(defaultCell.getPointer(), mblNull);

		delete mbl;
		
		// output results
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	}
private:
};




#endif
