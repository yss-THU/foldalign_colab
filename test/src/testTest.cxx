#ifndef TESTTEST
#define TESTTEST

#include "test.cxx"
#include "../../src/exception.cxx"

class testTest {
public:
	testTest(int& passed, int& ran, int& expected, std::string& messages) {
		test check1(0, "Dummy");

		std::cout << "------------------------------------------------------------------------------" << std::endl;	
		std::cout << "Testing the test.cxx code. This should results in some warnings" << std::endl;
		std::cout << std::endl;
		std::cout << "All results from the dummy test should be ignored" << std::endl;
		int nTest = 16;
		test master(nTest, "test.cxx test");
	
		master.equal(1, 1);
		master.equal(check1.equal(1, 2, "Wrong on purpose"), false);

		master.notEqual(1, 2);
		master.notEqual(check1.notEqual(2,2, "Wrong on purpose"), true);

		master.isTrue(check1.isTrue(true));
		master.isFalse(check1.isTrue(false, "Wrong on purpose"));

		master.isFalse(check1.isFalse(true, "Wrong on purpose"));
		master.isTrue(check1.isFalse(false));

		master.equal("ABC", "ABC");
		master.notEqual("ABC", "CBA");

		master.equal(master.getTestsPassed(), 10, "10 passed");
		master.equal(master.getTestsRan(), 11, "11 ran");

		check1.failedTest("failedTest: Wrong on purpose");
		master.equal(check1.getTestsPassed(), 2);
		master.equal(check1.getTestsRan(), 7, "Check1 ran");
		
		master.equal(master.getTestsExpected(), nTest);
		std::cout << "Note: The exception handling tests are not well tested" << std::endl;
//		master.constructorThrowsException<test, exception, int>(5, "Supposed to fail");
		master.constructorDoesNotThrowException<test, int>(5);

		std::cout << "------------------------------------------------------------------------------" << std::endl;	
		master.printResult();

		passed += master.getTestsPassed();
		ran += master.getTestsRan();
		expected += master.getTestsExpected();
	};
};
#endif
