#ifndef TESTREADFILE
#define TESTREADFILE

#include "test.cxx"
#include "../../src/readfile.cxx"
#include "../../src/exception.cxx"

class testReadfile {
public:
	testReadfile(int& passed, int& ran, int& expected, std::string& messages) 
		: tester(21, "readfile test"), filename("test/auxfiles/readfile.cxx.aux") {
	
		testStdin();
				
		testNamedFile();
				
		testGet_line();
		testGet_line_failEmpty();
		testSkip_to_empty();
		
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	}
private:

	test tester;
	const std::string filename;

	//! \brief Check that the name is set correct when stdin is used
	inline void testStdin();

	//! \brief Test that the named files can be open and that an error is given
	//! if the file can not be opened.
	inline void testNamedFile();

	//! \brief Test the get_line member function
	inline void testGet_line();

	//! \brief Test the get_line_failEmpty member function
	inline void testGet_line_failEmpty();

	//! \brief Test the skip_to_empty member function
	inline void testSkip_to_empty();

	//! \brief Opens the named file (filename) and returns a pointer to the
	//! readfile object. The pointer must be destroyed with delete
	inline readfile* openFile();
};

void testReadfile::testStdin() {

	readfile stdin;
	
	tester.equal(stdin.name(), std::string("<STDIN>"));
	
}

void testReadfile::testNamedFile() {

	tester.constructorThrowsException<readfile, exception, std::string>(
		std::string("hopefullydoesnotexists"));
	
	readfile* file = openFile();
	tester.equal(filename, file->name());
	delete file;
}

void testReadfile::testGet_line() {

	readfile* file = openFile();

	std::string line;
	tester.isTrue(file->get_line(line));
	tester.equal(line, std::string("1"));
	tester.isTrue(file->get_line(line));
	tester.equal(line, std::string("2"));
//	tester.isTrue(file->get_line(line), std::string("Read line with \\r"));
//	tester.equal(line, std::string(""), std::string("line with \\r"));
	tester.isFalse(file->get_line(line));
	tester.equal(line, std::string(""));

	delete file;
}	

void testReadfile::testGet_line_failEmpty() {

	readfile* file = openFile();

	std::string line;
	tester.isTrue(file->get_line_failEmpty(line), std::string("A"));
	tester.equal(line, std::string("1"), std::string("B"));

	tester.isFalse(file->get_line_failEmpty(line), std::string("C"));
	tester.equal(line, std::string(""), std::string("D"));

	tester.isTrue(file->get_line_failEmpty(line), std::string("E"));
	tester.equal(line, std::string("2"), std::string("F"));

//	tester.isTrue(file->get_line(line), std::string("G"));
//	tester.equal(line, std::string("line with \\r"), std::string("H"));

	tester.isFalse(file->get_line_failEmpty(line), std::string("I"));
	tester.equal(line, std::string(""), std::string("J"));

	delete file;
}	

void testReadfile::testSkip_to_empty() {

	readfile* file = openFile();

	std::string line;
	
	file->skip_to_empty();
	tester.isTrue(file->get_line(line));
	tester.equal(line, std::string("2"));

	file->skip_to_empty();
	tester.isFalse(file->get_line(line));
	tester.equal(line, std::string(""));

	delete file;
}	

readfile* testReadfile::openFile() {

	readfile* file;

	try {
		file = new readfile(filename);
	}
	catch (...) {
		std::string error = "Could not open file: " + filename;
		tester.failedTest(error);
		throw;
	}

	return file;
}
#endif
