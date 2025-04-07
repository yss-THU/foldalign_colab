#ifndef TESTLONGTERMMEMORY
#define TESTLONGTERMMEMORY

#include "../../src/longTermMemory.cxx"
#include "test.cxx"
#include "../../src/longCell.cxx"
#include "../../src/foldalign.hxx"

class testLongTermMemory {
public:
	testLongTermMemory(int& passed, int& ran, int& expected, std::string& messages)
		: tester(test(302, "longTermMemory.cxx test")){
	
		// Test local alignment matrix
		testMatrix(false);
		
		// test global alignment matrix
		testMatrix(true, 1, 1000);
		
		// And now for the result
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	}


private:
	test tester;

	void testMatrix(const bool global,
							const positionType k_off = 0,
							const positionType i_fix = 0);
	void store(longTermMemory<longCell>& ltm,
			positionType i, positionType k,
			lengthType Wi, lengthType Wk,
			scoreType similarityScore, scoreType energyScore);
	void testGetNext(longTermMemory<longCell>& ltm,
			positionType i, positionType k,
			lengthType WiExpected, lengthType WkExpected,
			scoreType similarityScore, scoreType energyScore);
	void testGetNextNoReset(longTermMemory<longCell>& ltm,
			positionType i, positionType k,
			lengthType WiExpected, lengthType WkExpected,
			scoreType similarityScore, scoreType energyScore);
	void testGetNextReturnNullPtr(longTermMemory<longCell>& ltm,
			positionType i, positionType k);
};


void testLongTermMemory::testMatrix(const bool global,
							const positionType k_off,
							const positionType i_fix) {

	const lengthType i_dimension = 3;
	const positionType k_dimension = 4;
	positionType i_position = 3;
	positionType k_position = 4;
	const int number_of_threads = 0;
	
	longTermMemory<longCell> ltm(k_dimension, i_dimension,
			i_position, k_position, global, number_of_threads, k_off, i_fix);

	if (global) {
		ltm.setPositions(1,1);
	}

	store(ltm, 3, 4, 1, 2, 1, 2);	
	store(ltm, 3, 3, 0, 1, -1, 2);	
	store(ltm, 3, 2, 1, 0, 1, -2);	

	testGetNext(ltm, 3, 4, 1, 2, 1, 2);	
	testGetNextReturnNullPtr(ltm, 3, 4);
	testGetNext(ltm, 3, 3, 0, 1, -1, 2);	
	testGetNext(ltm, 3, 2, 1, 0, 1, -2);	

	// Test that the values are still there
	testGetNext(ltm, 3, 4, 1, 2, 1, 2);
	
	// Check that transfer works and that new data can added with deleting the
	// old
	ltm.transfer();
	
	// Add more data
	store(ltm, 2, 3, 45, 54, 23, -54);
	store(ltm, 2, 3, 51, 54, 21, 4);
	store(ltm, 2, 3, 51, 154, 1, 4);
	store(ltm, 2, 3, 52, 171, 21, 14);
	store(ltm, 2, 1, 5, 4, 51, -44);
	store(ltm, 2, 2, 6, 256, -41, -544);
	
	ltm.resetCurrent(2,3);
	testGetNextNoReset(ltm, 2, 3, 45, 54, 23, -54);
	testGetNextNoReset(ltm, 2, 3, 51, 54, 21, 4);
	testGetNextNoReset(ltm, 2, 3, 51, 154, 1, 4);
	testGetNextNoReset(ltm, 2, 3, 52, 171, 21, 14);
	testGetNextReturnNullPtr(ltm, 2, 3);
	ltm.unlock(2,3);

	testGetNext(ltm, 2, 1, 5, 4, 51, -44);
	testGetNext(ltm, 2, 2, 6, 256, -41, -544);

	// The old data is still there
	testGetNext(ltm, 3, 4, 1, 2, 1, 2);	
	testGetNext(ltm, 3, 3, 0, 1, -1, 2);	
	testGetNext(ltm, 3, 2, 1, 0, 1, -2);	

	// Check getPos
	lengthType Wi = 45;
	lengthType Wk = 54;
	longCell* get = ltm.getPos(2, 3, Wi, Wk);
	tester.equal(get->getSimilarityScore(), scoreType(23));
	tester.equal(get->getEnergyScore(), scoreType(-54));
	
	// Check getPos with invalid position
	Wi = 1; Wk = 1;
	get = 0;
	tester.equal(ltm.getPos(2, 3, Wi, Wk), get, "getPos did not return null ptr");
	
	// Check switchWi
	ltm.resetCurrent(2,3);
	testGetNextNoReset(ltm, 2, 3, 45, 54, 23, -54);
	testGetNextNoReset(ltm, 2, 3, 51, 54, 21, 4);
	ltm.switchWi(2,3);
	testGetNextNoReset(ltm, 2, 3, 52, 171, 21, 14);
	testGetNextReturnNullPtr(ltm, 2, 3);
	ltm.unlock(2,3);
	
	// Check deletePos
	ltm.resetCurrent(2,3);
	testGetNextNoReset(ltm, 2, 3, 45, 54, 23, -54);
	testGetNextNoReset(ltm, 2, 3, 51, 54, 21, 4);
	ltm.deletePos(2,3);
	testGetNextNoReset(ltm, 2, 3, 45, 54, 23, -54);
	testGetNextNoReset(ltm, 2, 3, 51, 54, 21, 4);
	testGetNextNoReset(ltm, 2, 3, 51, 154, 1, 4);
	testGetNextReturnNullPtr(ltm, 2, 3);
	ltm.unlock(2,3);

	// Check getLockAndResetSubMatrix
	// Get an empty structure
	ltm.transfer();
	two_link<longCell>* newTwoLink = ltm.getLockAndResetSubMatrix(1, 3);
	newTwoLink->resetCurrent();
	get = 0;
	tester.equal(newTwoLink->getNext(Wi, Wk), get, "Data structure not empty"); 
	ltm.unlock(1, 3);
	
	// Get a structure with data
	newTwoLink = ltm.getLockAndResetSubMatrix(2,3);
	testGetNextNoReset(ltm, 2, 3, 45, 54, 23, -54);
	testGetNextNoReset(ltm, 2, 3, 51, 54, 21, 4);
	testGetNextNoReset(ltm, 2, 3, 51, 154, 1, 4);
	testGetNextReturnNullPtr(ltm, 2, 3);
	ltm.unlock(2,3);
	
	if (global) {
		tester.equal(ltm.getIStart(), positionType(1000), "istart");
	}
	else {
		tester.equal(ltm.getIStart(), positionType(0), "istart");
	}

	// Check that all the data is still there and indirectly setPositions()
	testGetNext(ltm, 3, 4, 1, 2, 1, 2);	
	testGetNext(ltm, 3, 3, 0, 1, -1, 2);	
	testGetNext(ltm, 3, 2, 1, 0, 1, -2);	
	testGetNext(ltm, 2, 1, 5, 4, 51, -44);
	testGetNext(ltm, 2, 2, 6, 256, -41, -544);
	ltm.resetCurrent(2,3);
	testGetNextNoReset(ltm, 2, 3, 45, 54, 23, -54);
	testGetNextNoReset(ltm, 2, 3, 51, 54, 21, 4);
	testGetNextNoReset(ltm, 2, 3, 51, 154, 1, 4);
	testGetNextReturnNullPtr(ltm, 2, 3);
	ltm.unlock(2,3);	
}

void testLongTermMemory::store(longTermMemory<longCell>& ltm,
			positionType i, positionType k,
			lengthType Wi, lengthType Wk,
			scoreType similarityScore, scoreType energyScore) {
			
	longCell* store;
	store = ltm.putPos(i, k, Wi, Wk);
	store->set(similarityScore, energyScore);
}

void testLongTermMemory::testGetNext(longTermMemory<longCell>& ltm,
			positionType i, positionType k,
			lengthType WiExpected, lengthType WkExpected,
			scoreType similarityScore, scoreType energyScore) {
	
	ltm.resetCurrent(i,k);

	testGetNextNoReset(ltm, i, k, WiExpected, WkExpected, 
					   similarityScore, energyScore);
	ltm.unlock(i,k);
}

void testLongTermMemory::testGetNextNoReset(longTermMemory<longCell>& ltm,
			positionType i, positionType k,
			lengthType WiExpected, lengthType WkExpected,
			scoreType similarityScore, scoreType energyScore) {
	
	lengthType Wi;
	lengthType Wk;
	longCell* read = ltm.getNextPos(i, k, Wi, Wk);
	tester.equal(Wi, WiExpected, "Wi");
	tester.equal(Wk, WkExpected, "Wk");
	tester.equal(read->getSimilarityScore(), similarityScore, "similarityScore");
	tester.equal(read->getEnergyScore(), energyScore, "energyScore");
}

void testLongTermMemory::testGetNextReturnNullPtr(longTermMemory<longCell>& ltm,
			positionType i, positionType k) {

	lengthType Wi;
	lengthType Wk;
	longCell* nullLongCellPtr = 0;
	tester.equal(ltm.getNextPos(i, k, Wi, Wk), nullLongCellPtr, "getNextPos did not return null pointer");
	tester.equal(Wi, lengthType(0));
	tester.equal(Wk, lengthType(0));
}

#endif
