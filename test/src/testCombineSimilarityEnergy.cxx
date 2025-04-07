#ifndef TESTCOMBINESIMILARITYENERGY
#define TESTCOMBINESIMILARITYENERGY

#include "../../src/combineSimilarityEnergy.cxx"
#include "../../src/foldalign.hxx"
#include "test.cxx"

class testCombineSimilarityEnergy {
public:
	testCombineSimilarityEnergy(int& passed, int& ran, int& expected, std::string& messages) {
		test tester(4, "combineSimilarityEnergy test");
		
		tester.equal(combineSimilarityEnergy(1, 2), scoreType(3));
		tester.equal(combineSimilarityEnergy(big_neg, 2), big_neg);
		tester.equal(combineSimilarityEnergy(2, big_neg), big_neg);
		tester.equal(combineSimilarityEnergy(big_neg, big_neg), big_neg);
		
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	};

};


#endif
