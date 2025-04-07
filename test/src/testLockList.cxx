#ifndef TESTLOCKLIST
#define TESTLOCKLIST

#include "../../src/lockList.cxx"
#include "../../src/exception.cxx"
#include "test.cxx"

class testLockList {
public:
	testLockList(int& passed, int& ran, int& expected, std::string& messages) 
		: tester(2, "lockList.cxx test") {

		messages += "lockList.cxx is being build and does not currently cover much of the code\n";
	
		int numberLocks = 1;
		{
			lockList lock(numberLocks);
			try {
				// Unlock with no lock
				lock.unlock(1,1);
			}
			catch (exception& exp) {
				std::string expected = "Program error! lockList unlock. i=1 k=1\nlock 0 ";
				std::string errorMessage = "LockList exception message:\n" + exp.getMessage() + " Expected:\n" + expected;
				tester.equal(exp.getMessage().substr(0,47), expected, errorMessage);
			}

			lock.lock(1,1);
			lock.unlock(1,1);
			try {	
				// Unlock the same lock more times than it was locked
				lock.unlock(1,1);
			}
			catch (exception& exp) {
				std::string expected = "Program error! lockList unlock. i=1 k=1\nlock 0 1 1 0\n";
				std::string errorMessage = "LockList exception message:\n" + exp.getMessage() + " Expected:\n" + expected;
				tester.equal(exp.getMessage(), expected, errorMessage);
			}

			try {
				// Lock more locks than available
				lock.lock(1,1);
				lock.lock(2,1);
			}
			catch (exception& exp) {
				std::cerr << "Got unexpected exception of normal kind" << std::endl;
			}
			catch (...) {
				std::cerr << "Got unexpected exception of unknow kind" << std::endl;
			}

		}

		numberLocks = 2;
		lockList lock(numberLocks);
		lock.lock(1,1);
		lock.lock(1,2);
		lock.unlock(1,1);
		lock.unlock(1,2);
		

		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();

	}

private:
	test tester;

	
};

#endif
