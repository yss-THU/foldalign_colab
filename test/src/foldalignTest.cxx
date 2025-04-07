#ifndef FOLDALIGNTEST
#define FOLDALIGNTEST

#include "test.cxx"
#include "testTest.cxx"
#include "testArguments.cxx"
#include "testCell.cxx"
#include "testCombineSimilarityEnergy.cxx"
#include "testConstraints.cxx"
#include "testException.cxx"
#include "testGlobalPrune.cxx"
#include "testHelper.cxx"
#include "testLockList.cxx"
#include "testLongCell.cxx"
#include "testLongTermMemory.cxx"
#include "testMatrix2d.cxx"
#include "testMatrix4d.cxx"
#include "testOptions.cxx"
#include "testPrune.cxx"
#include "testPruneTable.cxx"
#include "testReadConstraints.cxx"
#include "testReadfile.cxx"
#include "testScorematrix.cxx"
#include "testTwo_link.cxx"

int main() {

	int passed = 0;
	int ran = 0;
	int expected = 0;
	std::string messages = "";
	try {
		testTest(passed, ran, expected, messages);
		testException(passed, ran, expected, messages);

		testOptions(passed, ran, expected, messages);
		testArguments(passed, ran, expected, messages);
		testCell(passed, ran, expected, messages);
		testCombineSimilarityEnergy(passed, ran, expected, messages);
		testConstraints(passed, ran, expected, messages);
		testGlobalPrune(passed, ran, expected, messages);
		testHelper(passed, ran, expected, messages);
//		testLockList lockList(passed, ran, expected, messages);
		testLongCell(passed, ran, expected, messages);
		testLongTermMemory(passed, ran, expected, messages);
		testMatrix2d(passed, ran, expected, messages);
		testMatrix4d(passed, ran, expected, messages);
		testPrune(passed, ran, expected, messages);
		testPruneTable(passed, ran, expected, messages);
		
		testReadConstraints(passed, ran, expected, messages);
		testReadfile(passed, ran, expected, messages);
	
		testScorematrix(passed, ran, expected, messages);
/*		testTwo_link two_link(passed, ran, expected, messages);*/
	}
	catch (exception excep) {
		std::cerr << std::endl;
		std::cerr << "Caught exception with message:" << std::endl;
	 	std::cerr << excep.getMessage() << std::endl;
		if (excep.getFatal()) {
			std::cerr << "The exception is set to be fatal" << std::endl;
		}
		else {
			std::cerr << "The exception is set to be non-fatal" << std::endl;
		}
		throw;
	}
	catch ( ... ) {
		std::cerr << "Unknown exception thrown during testing" << std::endl;
	}

	std::cout << std::endl;
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << "Summary. Passed: " << passed << " out of " << ran << " tests. Expected to run: " << expected << " tests" << std::endl;
	if (messages.compare("")) {
		std::cout << std::endl << "The following messages were detected:" << std::endl;
		std::cout << "-------------------------------------" << std::endl;
		std::cout << messages << std::endl;
	}
}


#endif
