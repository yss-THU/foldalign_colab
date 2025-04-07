#ifndef TESTEXCEPTION
#define TESTEXCEPTION

#include <string>
#include <iostream>
#include "../../src/exception.cxx"
#include "test.cxx"

class testException {
public:
	testException(int& passed, int& ran, int& expected, std::string& messages) {
	
		test tester(6, "exception.cxx");
		
		exception excep;
		tester.equal(excep.getMessage(), std::string("Unknown error encountered. Foldalign will try to continue."));
		tester.isFalse(excep.getFatal());
		
		std::string message = "TEST TEST";
		exception excep2(message);
		tester.equal(excep2.getMessage(), message);
		tester.isFalse(excep2.getFatal());

		
		message = "TEST TEST 2";
		exception excep3(message, true);
		tester.equal(excep3.getMessage(), message);
		tester.isTrue(excep3.getFatal());
		
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	}
};

#endif
