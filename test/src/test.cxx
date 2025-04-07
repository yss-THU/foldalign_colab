#ifndef TEST
#define TEST

#include<string>
#include<iostream>

//! \file test.cxx \brief Holds the test class

//! \class test
//! \brief Keeps track of the number of tests done, the number tests passed, and
//!	the expected number of tests

class test {
public:
	/*!	\brief Constructor
	
		\param tests sets the expected number of tests. This is used to warm if more
		   or less tests than expected are performed
		\param testName Names the test. This is to make it easy to grep out test
			results from other output
	*/
	test(int tests = 0, std::string testName = "unknown test") 
		: passed(0), numberOfTests(tests), actualTestCount(0),
		  numberOfTestsDefined(numberOfTests == 0 ? false : true),
		  name(testName)
		{};
	
	//! \brief Number of tests passed with correct result
	inline int getTestsPassed() const {return passed;}
	
	//! \brief Number of tests which have actually been run
	inline int getTestsRan() const {return actualTestCount;}

	//! \brief The expected number of tests to run
	inline int getTestsExpected() const {return numberOfTests;}
	
	//! \brief Compare value and expected. The test passes if they are equal.
	//! returns true is they are equal
	template< class testType > bool equal(
		const testType& value, const testType& expected,
		std::string failMessage = "");
	
	//! \brief Compare value and expected. The test passes if they are not equal.
	//! returns true is they are not equal
	template< class testType > bool notEqual(
		const testType& value, const testType& expected,
		std::string failMessage = "");

	//! \brief The test passes if the value is true. Then the return value is
	//! also true
	template< class testType > bool isTrue(
		const testType& value,
		std::string failMessage = "") 
		{return equal(value, true, failMessage);}

	//! \brief The test passes if the value is false. Then the return value is
	//! true
	template< class testType > bool isFalse(
		const testType& value,
		std::string failMessage = "") 
		{return equal(value, false, failMessage);}
	
	//! \brief Increase the actualTestCount by one without chaning the passed
	//! number of tests
	inline bool failedTest(std::string failMessage = "");
	
	//! \brief The test passes if the constructor throws an exception of type 
	//! exception
	template< class functionCall, class exception, class optionType >
		bool constructorThrowsException(
			optionType option, std::string failMessage = "");


	template< class functionCall, class optionType >
		bool constructorDoesNotThrowException(
			optionType option, std::string failMessage="");

	//! \breif Print the results
	inline void printResult() const;
	
private:

	int passed; //! \brief Number of passed tests
	int numberOfTests; //! \brief If numberOfTestsDefined then
					   //! this is the expected number of tests otherwise it is
					   //! equal to actualTestCount
	int actualTestCount; //! \brief The actual number of test done so far
	const bool numberOfTestsDefined; //! true if numberOfTests is set by the constructer
	std::string name; //! The name of the test

	//! Updates numberOfTests and actualTestCount when a new test is done
	inline void newTest();

	//! Prints the failMessage if it exists otherwise prints a newline
	inline void printFailMessage(std::string failMessage = "");
};

template< class testType > bool test::equal(
		const testType& value, const testType& expected,
		std::string failMessage) {
		
	newTest();
			
	if (value == expected) {
		passed++;
		return true;
	}
	else {
		std::cout << name << ": Values not equal in test number: " << actualTestCount;
		printFailMessage(failMessage);
		return false;
	}
}

template< class testType > bool test::notEqual(
		const testType& value, const testType& expected,
		std::string failMessage) {
		
	newTest();
			
	if (value != expected) {
		passed++;
		return true;
	}
	else {
		std::cout << name << ": Values equal in test number: " << actualTestCount;
		printFailMessage(failMessage);
		return false;
	}
}

bool test::failedTest(std::string failMessage) {

	newTest();
	
	std::cout << name << ": Failed test in test number: " << actualTestCount;
	printFailMessage(failMessage);
	return false;
}
	


template< class functionCall, class exception, class optionType >
	bool test::constructorThrowsException(optionType option, std::string failMessage) {

	newTest();
	
	try {
		functionCall function(option);
	}
	catch ( exception testException ) {
		passed++;
		return true;
	}
	catch ( ... ) {
		std::cout << name << ": The wrong type of execption got thrown in test number: " << actualTestCount;
		printFailMessage(failMessage);
		return false;
	}

	std::cout << name << ": No execption thrown in test number: " << actualTestCount;
	printFailMessage(failMessage);
	return false;
}

template< class functionCall, class optionType >
	bool test::constructorDoesNotThrowException(optionType option, std::string failMessage) {

	newTest();
	
	try {
		functionCall function(option);
	}
	catch ( ... ) {
		std::cout << name << ": The wrong type of execption got thrown in test number: " << actualTestCount;
		printFailMessage(failMessage);
		return false;
	}

	passed++;
	return true;
}

void test::printResult() const {

	bool error = false;
	if (passed == actualTestCount) {
		std::cout << name << ": All " << passed << " tests passed.";
	}
	else {
		std::cout << name << ": " << passed << " out of " << actualTestCount << " passed.";
		error = true;
	}
		
	if (numberOfTestsDefined) {
		if (actualTestCount == numberOfTests) {
			std::cout << " All tests have run";
		}
		else if (actualTestCount < numberOfTests) {
			std::cout << " Only " << actualTestCount << " out of " << numberOfTests << " have been run";
			error = true;
		}
		else {
			std::cout << " Ran " << actualTestCount << " tests but only expected to run " << numberOfTests << " tests";
			error = true;
		}
	}
	if (error) {
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void test::newTest() {

	actualTestCount++;
	if (numberOfTestsDefined) {
		if (actualTestCount > numberOfTests) {
			std::cerr << name <<": Test warning: There are more tests than expected." << std::endl;
			std::cerr << name << ":\tExpected: " << numberOfTests << " got: " << actualTestCount << " so far" << std::endl;
		}
	}
	else {
		numberOfTests++;
	}
}

void test::printFailMessage(std::string failMessage) {
	if (failMessage != "") {
		std::cout << " " << failMessage;
	}
	std::cout << std::endl;
}

#endif
