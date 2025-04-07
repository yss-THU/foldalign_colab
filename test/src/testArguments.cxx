#ifndef TESTARGUMENTS
#define TESTARGUMENTS

#include "../../src/arguments.cxx"
#include "test.cxx"
#include "../../src/exception.cxx"

#include <string>

class testArguments {
public:

	testArguments(int& passed, int& ran, int& expected, std::string& messages)
	 : tester(test(138, "arguments.cxx test")), 
		numberOfOptions(29),
		noOptionWithName("DOESNOTEXIST"),
		scoretypeTestName("-min_LS_score"),
		positiontypeTestName("-i"),
		lengthtypeTestName("-max_diff"),
		booltypeTestName("-global"),
		stringtypeTestName("-score_matrix"),
		orgScoreValue(0),
		orgPositionValue(1),
		orgLengthValue(25),
		orgStringValue("<default>"),
		orgBoolValue(false),
		newScoreValue(-10),
		secondNewScoreValue(42),
		newPositionValue(100),
		newLengthValue(54),
		newStringValue("test test"),
		newBoolValue(true),
		addScoreName("scoretype_name"),
		addPositionName("positiontype_name"),
		addLengthName("lengthtype_name"),
		addStringName("stringtype_name"),
		addBoolName("booltype_name"),
		addScoreValue(-37),
		addPositionValue(67),
		addLengthValue(454),
		addStringValue("test to the west"),
		addBoolValue(true)	{

		messages += "Constructing arguments.cxx with other than default values, and/or filenames is not tested";
	
		basicTest();
	
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();

	}

private:
	test tester;
	
	const int numberOfOptions;

	const std::string noOptionWithName;
	const std::string scoretypeTestName;
	const std::string positiontypeTestName;
	const std::string lengthtypeTestName;
	const std::string booltypeTestName;
	const std::string stringtypeTestName;

	scoreType orgScoreValue;
	positionType orgPositionValue;
	lengthType orgLengthValue;
	std::string orgStringValue;
	bool orgBoolValue;

	scoreType newScoreValue;
	scoreType secondNewScoreValue;
	positionType newPositionValue;
	lengthType newLengthValue;
	std::string newStringValue;
	bool newBoolValue;

	std::string addScoreName;
	std::string addPositionName;
	std::string addLengthName;
	std::string addStringName;
	std::string addBoolName;

	scoreType addScoreValue;
	positionType addPositionValue;
	lengthType addLengthValue;
	std::string addStringValue;
	bool addBoolValue;


	void basicTest();

	void testExistOpt(arguments& arg);
	void testGetters(arguments& arg);
	template< class returnType >
	void testGet(arguments& arg, std::string optName, 
			returnType expected, returnType (arguments::*func)(std::string optName) );
	void testDefaultValues(arguments& arg);

	void testSetters(arguments& arg);
	template< class newValueType >
	void testSet(arguments& arg, std::string optName, newValueType newValue,
		void (arguments::*set)(const std::string& optName, const newValueType& newValue),
		newValueType (arguments::*get)(std::string optName));

	void testAdders(arguments& arg);
	template< class newValueType >
	void testAdder(arguments& arg, std::string newParameter, newValueType value,
		void (arguments::*add)(const std::string&, const newValueType&),
		bool (arguments::*exists)(const std::string&) const,
		newValueType (arguments::*get)(std::string));

	void testAssignment(arguments& arg);
	void testCopy(arguments& arg);
	void testChangedValues(arguments& arg);
	void makeChangeTestDifference(arguments& arg, arguments& arg2);
};


void testArguments::basicTest() {

	// Test without any input (i.e. default values)

	const int argc = 1;
	char* argv[argc];
	char name[25] = "testArguments Assignment";
	argv[0] = name;
	
	arguments arg(argc, argv);
	
	tester.equal(arg.numberOfOptions(), numberOfOptions, "Basic numberOfOptions");
	tester.equal(arg.numberOfArguments(), 0, "Basic numberOfArguments");

	testExistOpt(arg);
	testGetters(arg);
	
	testDefaultValues(arg);
	
	tester.equal(arg.getDefault(stringtypeTestName), orgStringValue, "Basic getDefault");
	tester.equal(arg.optionNumber(3), stringtypeTestName, "Basic optionNumber");
	tester.equal(arg.descriptionNumber(3), std::string("Path and filename of a score matrix file. There are two buildt in matrices one for local alignment (the default) and one for global alignment. The global alignment matrix is used when option -global is used."), "Basic descriptionNumber");

	testSetters(arg);
	testAdders(arg);

	testAssignment(arg);
	testCopy(arg);
}

void testArguments::testExistOpt(arguments& arg) {

	tester.equal(arg.existStOpt(scoretypeTestName), true, "Basic existStOpt");
	tester.equal(arg.existStOpt(noOptionWithName), false, "Basic existStOpt");
	tester.equal(arg.existPtOpt(positiontypeTestName), true, "Basic existPtOpt");
	tester.equal(arg.existPtOpt(noOptionWithName), false, "Basic existPtOpt");
	tester.equal(arg.existLtOpt(lengthtypeTestName), true, "Basic existLtOpt");
	tester.equal(arg.existLtOpt(noOptionWithName), false, "Basic existLtOpt");
	tester.equal(arg.existStringOpt(stringtypeTestName), true, "Basic existStringOpt");
	tester.equal(arg.existStringOpt(noOptionWithName), false, "Basic existStringOpt");
	tester.equal(arg.existBoolOpt(booltypeTestName), true, "Basic existBoolOpt");
	tester.equal(arg.existBoolOpt(noOptionWithName), false, "Basic existBoolOpt");

}

void testArguments::testGetters (arguments& arg) {

	testGet(arg, std::string(booltypeTestName), orgBoolValue, &arguments::boolOpt);
	testGet(arg, std::string(positiontypeTestName), orgPositionValue, &arguments::ptOpt);
	testGet(arg, std::string(lengthtypeTestName), orgLengthValue, &arguments::ltOpt);
	testGet(arg, std::string(scoretypeTestName), orgScoreValue, &arguments::stOpt);
	testGet(arg, std::string(stringtypeTestName), orgStringValue, &arguments::stringOpt);

}

template< class returnType >
void testArguments::testGet(arguments& arg, std::string optName, 
			returnType expected, returnType (arguments::*func)(std::string optName) ) {

	tester.equal((arg.*func)(optName), expected, "Basic get boolOpt");
	try {
		(arg.*func)(noOptionWithName);
		tester.failedTest("Basic getOpt did not throw an exception");
	}
	catch (exception& excp) {

		std::string expected = "Program error! Options.get. Option " + noOptionWithName + " does not exists.";
		tester.equal(excp.getMessage(), expected, "Basic getOpt exception" + excp.getMessage());
		tester.equal(excp.getFatal(), true, "Basic getOpt fatal");
	}
	catch (...) {
		tester.failedTest("Basic getOpt throws the wrong exception");
	}

}

void testArguments::testDefaultValues(arguments& arg) {

	try {
		tester.equal(arg.ltOpt("-max_diff"), lengthType(25));
		tester.equal(arg.ltOpt("-max_length"), lengthType(-1));
		tester.equal(arg.ltOpt("-min_loop"), lengthType(3));
		tester.equal(arg.stringOpt("-score_matrix"), std::string("<default>"));
//		tester.equal(arg.stringOpt("-seed_constraints"), std::string("<none>"));
//		tester.equal(arg.ltOpt("-min_seed_expand_length"), lengthType(-1));
//		tester.equal(arg.stringOpt("-backtrack_seed_constraints"), std::string("<none>"));
		tester.equal(arg.ltOpt("-chunk_size"), lengthType(-1));
		tester.equal(arg.boolOpt("-no_pruning"), false);
		tester.equal(arg.stringOpt("-ID"), std::string("n.a."));
		tester.equal(arg.boolOpt("-nobranch"), false);
		tester.equal(arg.boolOpt("-no_backtrack"), false);
		tester.equal(arg.boolOpt("-global"), false);
		tester.equal(arg.boolOpt("-use_global_pruning"), false);
		tester.equal(arg.ptOpt("-i"), positionType(1));
		tester.equal(arg.ptOpt("-j"), positionType(-1));
		tester.equal(arg.ptOpt("-k"), positionType(1));
		tester.equal(arg.ptOpt("-l"), positionType(-1));
		tester.equal(arg.stringOpt("-format"), std::string("fasta"));
		tester.equal(arg.stringOpt("-output_format"), std::string("column"));
		tester.equal(arg.boolOpt("-plot_score"), false);
//		tester.equal(arg.boolOpt("-print_seed_constraints"), false);
		tester.equal(arg.stOpt("-min_LS_score"), scoreType(0));
		tester.equal(arg.boolOpt("-print_all_LS_scores"), false);
		tester.equal(arg.ltOpt("-number_of_processors"), lengthType(1));
		tester.equal(arg.boolOpt("-align_self"), false);
		tester.equal(arg.boolOpt("-backtrack_info"), false);
		tester.equal(arg.boolOpt("-all_scores"), false);
		tester.equal(arg.stringOpt("-memory_roof"), std::string("-1"));
		tester.equal(arg.boolOpt("-memory_info"), false);
		tester.equal(arg.boolOpt("-version"), false);
		tester.equal(arg.boolOpt("-help"), false);
		tester.equal(arg.boolOpt("-h"), false);
	}
	catch (exception& exp) {
		std::cerr << exp.getMessage() << std::endl;
		throw;
	}
}

void testArguments::testSetters(arguments& arg) {

	testSet(arg, std::string(scoretypeTestName), newScoreValue,
			&arguments::setSt, &arguments::stOpt);
	testSet(arg, positiontypeTestName, newPositionValue,
			&arguments::setPt, &arguments::ptOpt);
	testSet(arg, lengthtypeTestName, newLengthValue,
			&arguments::setLt, &arguments::ltOpt);
	testSet(arg, booltypeTestName, newBoolValue,
			&arguments::setBool, &arguments::boolOpt);
	testSet(arg, stringtypeTestName, newStringValue,
			&arguments::setString, &arguments::stringOpt);
}

template< class newValueType >
void testArguments::testSet(arguments& arg, std::string optName,
	newValueType newValue,
	void (arguments::*set)(const std::string& optName, const newValueType& newValue),
	newValueType (arguments::*get)(std::string optName)) {

	(arg.*set)(optName, newValue);
	tester.equal((arg.*get)(optName), newValue);

	try {
		(arg.*set)(noOptionWithName, newValue);
		tester.failedTest("Basic set opt did not throw an exception");
	
	}
	catch (exception& excp) {

		std::string expected = "Program error! Options.set. Option " + noOptionWithName + " does not exists.";
		tester.equal(excp.getMessage(), expected, "Basic setOpt exception: " + excp.getMessage());
		tester.equal(excp.getFatal(), true, "Basic setOpt fatal");
	}
	catch (...) {
		tester.failedTest("Basic setOpt throws the wrong exception");
	}
}

void testArguments::testAdders(arguments& arg) {

	testAdder(arg, addScoreName, addScoreValue,
			&arguments::addSt, &arguments::existStOpt, &arguments::stOpt);
	testAdder(arg, addPositionName, addPositionValue,
			&arguments::addPt, &arguments::existPtOpt, &arguments::ptOpt);
	testAdder(arg, addLengthName, addLengthValue,
			&arguments::addLt, &arguments::existLtOpt, &arguments::ltOpt);
	testAdder(arg, addBoolName, addBoolValue,
			&arguments::addBool, &arguments::existBoolOpt, &arguments::boolOpt);
	testAdder(arg, addStringName, addStringValue,
			&arguments::addString, &arguments::existStringOpt,
			&arguments::stringOpt);

}

template< class newValueType >
void testArguments::testAdder(arguments& arg, std::string newParameter,
	newValueType value,
	void (arguments::*add)(const std::string&, const newValueType&),
	bool (arguments::*exists)(const std::string&) const,
	newValueType (arguments::*get)(std::string)) {

	(arg.*add)(newParameter, value);
	bool newParameterExists = (arg.*exists)(newParameter);
	tester.isTrue(newParameterExists);
	if (newParameterExists) {
		tester.equal((arg.*get)(newParameter), value);
	}
	else {
		std::string error = "Could not add new parameter: " + newParameter;
		throw exception(newParameter, true);
	}
}

void testArguments::testAssignment(arguments& arg) {
	
	// Test assignment operator
//	const int argc = 0;
//	char* argv[argc];

	const int argc = 1;
	char* argv[argc];
	char name[25] = "testArguments Assignment";
	argv[0] = name;

	arguments arg2(argc, argv);
	testDefaultValues(arg2);
	arg2 = arg;

	testChangedValues(arg2);

	makeChangeTestDifference(arg, arg2);
		
	// This part is not run in the basic testing.
	tester.equal(arg.numberOfArguments(), arg2.numberOfArguments(), "Assign numberOfArguments");
	for(int i = 0; i < arg2.numberOfArguments(); i++) {
		tester.equal(arg2.argument(i), arg.argument(i), "Assign arguments");
	}
}

void testArguments::testCopy(arguments& arg) {

	arguments copy = arg;
	
	testChangedValues(copy);

	makeChangeTestDifference(arg, copy);	
}

void testArguments::testChangedValues(arguments& arg) {

	// Check that the values changed in testSetters are now also changed for arg2
	tester.equal(arg.stOpt(scoretypeTestName), newScoreValue);
	tester.equal(arg.ptOpt(positiontypeTestName), newPositionValue);
	tester.equal(arg.ltOpt(lengthtypeTestName), newLengthValue);
	tester.equal(arg.boolOpt(booltypeTestName), newBoolValue);
	tester.equal(arg.stringOpt(stringtypeTestName), newStringValue);

	// Check that the values added in testAdders are now in arg
	tester.equal(arg.stOpt(addScoreName), addScoreValue);
	tester.equal(arg.ptOpt(addPositionName), addPositionValue);
	tester.equal(arg.ltOpt(addLengthName), addLengthValue);
	tester.equal(arg.boolOpt(addBoolName), addBoolValue);
	tester.equal(arg.stringOpt(addStringName), addStringValue);
}

void testArguments::makeChangeTestDifference(arguments& arg, arguments& arg2) {

	// Change an option in arg2 and check that the option in arg was not changed
	arg2.setSt(scoretypeTestName, secondNewScoreValue);
	tester.equal(arg2.stOpt(scoretypeTestName), secondNewScoreValue, "Change scoretype");
	tester.equal(arg.stOpt(scoretypeTestName), newScoreValue);

}
#endif
