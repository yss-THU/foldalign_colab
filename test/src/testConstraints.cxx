#ifndef TESTCONSTRAINTS
#define TESTCONSTRAINTS

#include "../../src/constraints.cxx"
#include "test.cxx"

class testConstraints {
public:

	testConstraints(int& passed, int& ran, int& expected, std::string& messages)
	 : tester(test(20, "constraints.cxx test")), auxFilePath("test/auxfiles/") {
	
		constraints cons(auxFilePath + "constraints.cxx.aux", false);
		
		// Test getMinScore
		tester.equal(cons.getMinScore(), scoreType(-13));

		testExistsEnclosedConstraints();

		testMinExpandLengthMany();
		testMinExpandLength2();
		testMinExpandLength1();
		testMinExpandLength0();


		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	}
private:

	test tester;

	const std::string auxFilePath;
	void testExistsEnclosedConstraints();

	void testMinExpandLengthMany();
	void testMinExpandLength2();
	void testMinExpandLength1();
	void testMinExpandLength0();
};

void testConstraints::testExistsEnclosedConstraints() {

	constraints cons(auxFilePath + "constraints.cxx.aux", false);

	tester.equal(cons.existsEnclosedConstraints(1,1,1,1), false);
	tester.equal(cons.existsEnclosedConstraints(600,600,600,600), false);

	tester.equal(cons.existsEnclosedConstraints(486,21,12,12), true);
	tester.equal(cons.existsEnclosedConstraints(486,21,12,13), true);
	tester.equal(cons.existsEnclosedConstraints(486,20,12,13), true);
	tester.equal(cons.existsEnclosedConstraints(486,20,12,14), true);


//	tester.equal(cons.existsEnclosedConstraints(476,44,13,0), true);
//	tester.equal(cons.existsEnclosedConstraints(476,44,13,11), true);
}

void testConstraints::testMinExpandLengthMany() {

	constraints cons(auxFilePath + "constraints.cxx.aux", false);

	tester.equal(cons.getSeedMinExpandLength(201), lengthType(100));
	tester.equal(cons.getSeedMinExpandLength(200), lengthType(30));
	tester.equal(cons.getSeedMinExpandLength(101), lengthType(30));
	tester.equal(cons.getSeedMinExpandLength(100), lengthType(0));
	tester.equal(cons.getSeedMinExpandLength(31), lengthType(0));
	tester.equal(cons.getSeedMinExpandLength(30), lengthType(0));
}

void testConstraints::testMinExpandLength2() {

	constraints cons(auxFilePath + "constraints.cxx.maxLength.2.aux", false);

	tester.equal(cons.getSeedMinExpandLength(101), lengthType(30));
	tester.equal(cons.getSeedMinExpandLength(100), lengthType(0));
	tester.equal(cons.getSeedMinExpandLength(31), lengthType(0));
	tester.equal(cons.getSeedMinExpandLength(30), lengthType(0));
}

void testConstraints::testMinExpandLength1() {

	constraints cons(auxFilePath + "constraints.cxx.maxLength.1.aux", false);

	tester.equal(cons.getSeedMinExpandLength(31), lengthType(0));
	tester.equal(cons.getSeedMinExpandLength(30), lengthType(0));
}

void testConstraints::testMinExpandLength0() {

	constraints cons(auxFilePath + "constraints.cxx.maxLength.0.aux", false);

	tester.equal(cons.getSeedMinExpandLength(30), lengthType(0));
}

#endif
