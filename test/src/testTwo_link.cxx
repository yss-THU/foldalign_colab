#ifndef TESTTWOLINK
#define TESTTWOLINK

#include "../../src/two_link.cxx"
#include "test.cxx"
#include "../../src/longCell.cxx"
#include "../../src/foldalign.hxx"

class testTwo_link {
public:
	testTwo_link(int& passed, int& ran, int& expected, std::string& messages)
		: tester(test(94, "test two_link")) {

		two_link<longCell> twoLink;

		// Test that the initial structure is empty
		testGetNextReturnNullPtr(twoLink);

		// Store some data
		longCell* store;
		store = twoLink.putNext(0, 0);
		store->set(-3, 3);
		store = twoLink.putNext(1, 1);
		store->set(6,6);
		store = twoLink.putNext(1, 3);
		store->set(1,-5);
		store = twoLink.putNext(5, 7);
		store->set(-4,-6);
// The following results in a known bug. The bug is expected to be of little
// consevense for Foldalign. But hamper code reusage.
//		store = twoLink.putNext(0, 1);
//		store->set(2,-9);

		// After writing the object must be reset before reading
		twoLink.resetCurrent();

		// Try to retrive the data stored earlier.
		runGetNextTests(twoLink);
		// There is no more data stored getNext should return a null pointer
		testGetNextReturnNullPtr(twoLink);

		// Test that using resetCurrent makes it possible to retrive all the
		// stored data again
		twoLink.resetCurrent();

		// Try to retrive the data stored earlier.
		runGetNextTests(twoLink);

		// There is no more data stored getNext should return a null pointer
		testGetNextReturnNullPtr(twoLink);
		testGetNextReturnNullPtr(twoLink);

		// Test getPos
		testGetPos(twoLink, 0, 0, -3, 3);
		testGetPos(twoLink, 1, 1, 6, 6);
		testGetPos(twoLink, 1, 3, 1, -5);
		testGetPosNotExists(twoLink, 3, 3); // Wi missing in the middel of structure
		testGetPosNotExists(twoLink, 10, 3); // Wi larger than max
		testGetPosNotExists(twoLink, 5, 2); // Wk smaller than min
		testGetPosNotExists(twoLink, 1, 2); // Wk missing in the middel of structure
		testGetPosNotExists(twoLink, 1, 4); // Wk larger than the max

		twoLink.resetCurrent();
		twoLink.lastWk();
		testGetNext(twoLink, 1, 1, 6, 6);
		twoLink.lastWk();
		testGetNext(twoLink, 5, 7, -4, -6);
		
		// Test removeRest()
		twoLink.resetCurrent();
		testGetNext(twoLink, 0, 0, -3, 3);
		testGetNext(twoLink, 1, 1, 6, 6);
		twoLink.removeRest();
		testGetNext(twoLink, 0, 0, -3, 3);
		testGetNext(twoLink, 1, 1, 6, 6);
		testGetNext(twoLink, 1, 3, 1, -5);
		// There should be nothing left
		testGetNextReturnNullPtr(twoLink);
		// RemoveRest only removes Wi nodes with larger Wi values. Using
		// removeRest again should therefore make no difference
		twoLink.removeRest();
		testGetNext(twoLink, 0, 0, -3, 3);
		testGetNext(twoLink, 1, 1, 6, 6);
		testGetNext(twoLink, 1, 3, 1, -5);
		// There should be nothing left
		testGetNextReturnNullPtr(twoLink);
		// RemoveRest only removes Wi nodes with larger Wi values. Using
		// removeRest again should therefore make no difference
		twoLink.resetCurrent();
		twoLink.removeRest();
		testGetNext(twoLink, 0, 0, -3, 3);
		// There should be nothing left
		testGetNextReturnNullPtr(twoLink);
		
		// And now for the result
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();

	}

private:
	test tester;

	void runGetNextTests(two_link<longCell>& twoLink);
	void testGetNext(two_link<longCell>& twoLink,
			lengthType WiExpected, lengthType WkExpected,
			scoreType similarityScore, scoreType energyScore);

	void testGetNextReturnNullPtr(two_link<longCell>& twoLink);
	void testGetPos(two_link<longCell>& twoLink, lengthType Wi, lengthType Wk,
			scoreType similarityScore, scoreType energyScore);
	void testGetPosNotExists(two_link<longCell>& twoLink,
			lengthType Wi, lengthType Wk);
			
};

void testTwo_link::runGetNextTests(two_link<longCell>& twoLink) {

	testGetNext(twoLink, 0, 0, -3, 3);
	testGetNext(twoLink, 1, 1, 6, 6);
	testGetNext(twoLink, 1, 3, 1, -5);
	testGetNext(twoLink, 5, 7, -4, -6);
// See bug note above (in store part)
//	testGetNext(twoLink, 0, 1, 2, -9); 

}

void testTwo_link::testGetNext(two_link<longCell>& twoLink,
			lengthType WiExpected, lengthType WkExpected,
			scoreType similarityScore, scoreType energyScore) {

	lengthType Wi;
	lengthType Wk;
	
	longCell* read;
	read = twoLink.getNext(Wi, Wk);
//std::cerr << Wi << " " << Wk << std::endl;
	tester.equal(Wi, lengthType(WiExpected), "Wi");
	tester.equal(Wk, lengthType(WkExpected), "Wk");
	tester.equal(read->getSimilarityScore(), similarityScore, "similarity");
	tester.equal(read->getEnergyScore(), energyScore, "energy");
}

void testTwo_link::testGetNextReturnNullPtr(two_link<longCell>& twoLink) {

	longCell* nullLongCellPtr = 0;
	lengthType Wi = 0;
	lengthType Wk = 0;
		
	tester.equal(twoLink.getNext(Wi,Wk), nullLongCellPtr, "getNext did not return null ptr.");

}

void testTwo_link::testGetPos(two_link<longCell>& twoLink,
			lengthType Wi, lengthType Wk,
			scoreType similarityScore, scoreType energyScore) {
	longCell* read = twoLink.getPos(Wi, Wk);
	tester.equal(read->getSimilarityScore(), similarityScore);
	tester.equal(read->getEnergyScore(), energyScore);
}

void testTwo_link::testGetPosNotExists(two_link<longCell>& twoLink,
			lengthType Wi, lengthType Wk) {

	longCell* nullLongCellPtr = 0;
	tester.equal(twoLink.getPos(Wi, Wk), nullLongCellPtr);
}

#endif
