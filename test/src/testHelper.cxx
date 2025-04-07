#ifndef TESTHELPER
#define TESTHELPER

#include "test.cxx"
#include "../../src/helper.cxx"

class testHelper {
public:
	testHelper(int& passed, int& ran, int& expected, std::string& messages) {

		messages = "helper.cxx is not completely tested yet\n";
	
		test tester(10, "helper.cxx test");
		
		// Test assign
		const int size = 4;
		int t[size] = {11, 12, 13, 14};
		helper::assign(t, 2, 5, 7, 12);
		tester.equal(2, t[0]);
		tester.equal(5, t[1]);
		tester.equal(7, t[2]);
		tester.equal(12, t[3]);
		
		// Test swap
		helper::swap(t[0], t[1]);
		tester.equal(5, t[0]);
		tester.equal(2, t[1]);
		
		// Test copy array
		int s[size];
		helper::copyArray(s, t, size);
		tester.equal(5, t[0]);
		tester.equal(2, t[1]);
		tester.equal(7, t[2]);
		tester.equal(12, t[3]);
		
		// Print the results;
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	}

};
#endif
