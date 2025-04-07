#ifndef TESTREADCONSTRAINTS
#define TESTREADCONSTRAINTS

#include "test.cxx"
#include "../../src/readConstraints.cxx"
#include "../../src/foldalign.hxx"
#include "../../src/jl.cxx"

class testReadConstraints {
public:
	typedef std::multimap< positionType, jl > mmap_k;
	typedef std::map< positionType, mmap_k > map_i;

	testReadConstraints(int& passed, int& ran, int& expected, std::string& messages)
	 : tester(test(180, "readConstraints.cxx test")), 
	   auxFilePath("test/auxfiles/") {
		
		testMissingFileException();
		
		std::string filename = auxFilePath + "readConstraints.cxx.noSwap.aux";
		testFile(filename, false);
		
		filename = auxFilePath + "readConstraints.cxx.swap.aux";
		testFile(filename, true);
		
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	};

private:

	test tester;

	const std::string auxFilePath;

	void testMissingFileException();
	void testFile(std::string filename, bool swap);

};

void testReadConstraints::testMissingFileException() {

	std::string filename = "missingFile";

	map_i testMap_i;
	scoreType minScore = 0;
	std::vector<lengthType> testMaxLengths;	

	try {
		readconstraints<jl> read(filename, testMap_i, false, minScore, testMaxLengths);
	}
	catch (exception& error) {
		std::string expectedErrorMessage = "It is not possible to open the constraint file: " + filename;
		tester.equal(error.getMessage(), expectedErrorMessage);
		tester.equal(error.getFatal(), false);
		return;
	}
	catch (...) {
		tester.failedTest("exception test: Wrong exception");
		return;
	}
	tester.failedTest("exception test: No exception");
}

void testReadConstraints::testFile(std::string filename, bool swap) {

	map_i cons;
	scoreType minScore = 0;
	std::vector<lengthType> maxLengths;	

	try {
		readconstraints<jl> read(filename, cons, swap, minScore, maxLengths);
	}
	catch (...) {
		std::string error = "Unexpected exception trying to read constraints file: " + filename;
		tester.failedTest(error);
		return;
	}

	// Test minScore
	tester.equal(minScore, scoreType(-13), "minScore");

	// Test testMaxLengths
	tester.equal(int(maxLengths.size()), 3);
	tester.equal(maxLengths[0], lengthType(30), "maxLength 0");
	tester.equal(maxLengths[1], lengthType(100), "maxLength 1");
	tester.equal(maxLengths[2], lengthType(200), "maxLength 2");

	// Test map_i
	tester.equal(int(cons.size()), 6);

	tester.equal(int(cons[486].size()), 1, "486");
	tester.equal(int(cons[482].size()), 3, "482");
	tester.equal(int(cons[479].size()), 1, "479");
	tester.equal(int(cons[476].size()), 4, "476");
	tester.equal(int(cons[478].size()), 1, "478");
	tester.equal(int(cons[477].size()), 1, "477");
	
	int count = 0;
	for(map_i::iterator itI = cons.begin(); itI != cons.end(); itI++) {
		for(mmap_k::iterator itK = itI->second.begin(); itK != itI->second.end(); itK++) {
			switch (count) {
				case 0:
					tester.equal(itI->first, positionType(476));
					tester.equal(itK->second.j, positionType(489));
					tester.equal(itK->first, positionType(1));
					tester.equal(itK->second.l, positionType(14));
					tester.equal(itK->second.similarityScore, scoreType(28));
					tester.equal(itK->second.energyScore, scoreType(63));
					tester.equal(itK->second.state, stateType(88));
					break;
				case 1:
					tester.equal(itI->first, positionType(476));
					tester.equal(itK->second.j, positionType(489));
					tester.equal(itK->first, positionType(44));
					tester.equal(itK->second.l, positionType(57));
					tester.equal(itK->second.similarityScore, scoreType(-1));
					tester.equal(itK->second.energyScore, scoreType(33));
					tester.equal(itK->second.state, stateType(1));
					break;
				case 2:
					tester.equal(itI->first, positionType(476));
					tester.equal(itK->second.j, positionType(476));
					tester.equal(itK->first, positionType(44));
					tester.equal(itK->second.l, positionType(57));
					tester.equal(itK->second.similarityScore, scoreType(-16));
					tester.equal(itK->second.energyScore, scoreType(3));
					tester.equal(itK->second.state, stateType(2));
					break;
				case 3:
					tester.equal(itI->first, positionType(476));
					tester.equal(itK->second.j, positionType(489));
					tester.equal(itK->first, positionType(44));
					tester.equal(itK->second.l, positionType(54));
					tester.equal(itK->second.similarityScore, scoreType(16));
					tester.equal(itK->second.energyScore, scoreType(33));
					tester.equal(itK->second.state, stateType(8));
					break;
				case 4:
					tester.equal(itI->first, positionType(477));
					tester.equal(itK->second.j, positionType(488));
					tester.equal(itK->first, positionType(2));
					tester.equal(itK->second.l, positionType(13));
					tester.equal(itK->second.similarityScore, scoreType(9));
					tester.equal(itK->second.energyScore, scoreType(6));
					tester.equal(itK->second.state, stateType(88));
					break;
				case 5:
					tester.equal(itI->first, positionType(478));
					tester.equal(itK->second.j, positionType(495));
					tester.equal(itK->first, positionType(24));
					tester.equal(itK->second.l, positionType(41));
					tester.equal(itK->second.similarityScore, scoreType(71));
					tester.equal(itK->second.energyScore, scoreType(-5));
					tester.equal(itK->second.state, stateType(88));
					break;
				case 6:
					tester.equal(itI->first, positionType(479));
					tester.equal(itK->second.j, positionType(494));
					tester.equal(itK->first, positionType(25));
					tester.equal(itK->second.l, positionType(40));
					tester.equal(itK->second.similarityScore, scoreType(64));
					tester.equal(itK->second.energyScore, scoreType(-54));
					tester.equal(itK->second.state, stateType(88));
					break;
				case 7:
					tester.equal(itI->first, positionType(482));
					tester.equal(itK->second.j, positionType(488));
					tester.equal(itK->first, positionType(17));
					tester.equal(itK->second.l, positionType(28));
					tester.equal(itK->second.similarityScore, scoreType(-1));
					tester.equal(itK->second.energyScore, scoreType(7), "482 488 17 28 -1 7 88");
					tester.equal(itK->second.state, stateType(88));
					break;
				case 8:
					tester.equal(itI->first, positionType(482));
					tester.equal(itK->second.j, positionType(489));
					tester.equal(itK->first, positionType(28));
					tester.equal(itK->second.l, positionType(35));
					tester.equal(itK->second.similarityScore, scoreType(59));
					tester.equal(itK->second.energyScore, scoreType(-21), "-21");
					tester.equal(itK->second.state, stateType(88));
					break;
				case 9:
					tester.equal(itI->first, positionType(482));
					tester.equal(itK->second.j, positionType(489));
					tester.equal(itK->first, positionType(32));
					tester.equal(itK->second.l, positionType(39));
					tester.equal(itK->second.similarityScore, scoreType(35));
					tester.equal(itK->second.energyScore, scoreType(-24), "-24");
					tester.equal(itK->second.state, stateType(88));
					break;
				case 10:
					tester.equal(itI->first, positionType(486));
					tester.equal(itK->second.j, positionType(498));
					tester.equal(itK->first, positionType(21));
					tester.equal(itK->second.l, positionType(33));
					tester.equal(itK->second.similarityScore, scoreType(45));
					tester.equal(itK->second.energyScore, scoreType(-34), "-34");
					tester.equal(itK->second.state, stateType(88));
					break;
				default:
					tester.failedTest("Unexpected number of constraints");
			}
			count++;
		}
	}
}
#endif
