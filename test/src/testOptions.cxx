#ifndef TESTOPTIONS
#define TESTOPTIONS

#include "../../src/options.cxx"
#include "../../src/convertString.cxx"
#include "test.cxx"

#include <string>
#include <sstream>
class testOptions {
public:

	testOptions(int& passed, int& ran, int& expected, std::string& messages)
	 : tester(41, "options.cxx"), size(6) {
		
		names = new std::string[size];
		values = new std::string[size];
		descriptions = new std::string[size];
		types = new std::string[size];

		names[0] = "name 1";
		names[1] = "name 2";
		names[2] = "name 3";
		names[3] = "name 4";
		names[4] = "name 5";
		names[5] = "name 6";
		values[0] = "Value 1";
		values[1] = "Value 2";
		values[2] = "2";
		values[3] = "-10";
		values[4] = "26";
		values[5] = "545";
		descriptions[0] = "description 1";
		descriptions[1] = "description 2";
		descriptions[2] = "description 3";
		descriptions[3] = "description 4";
		descriptions[4] = "description 5";
		descriptions[5] = "description 6";
		types[0] = "string";
		types[1] = "string";
		types[2] = "int";
		types[3] = "int";
		types[4] = "int";
		types[5] = "int";
		
		floatOpt = new options<float>(std::string("float"), size, names, values,
				types, descriptions, &convertString::toNumeric<float>);
				
		delete floatOpt;
		
		testgetSetExistsAddWithString();

		options<int> intOpt(std::string("int"), size, names, 
				values, types, descriptions, &convertString::toNumeric<int>);

		tester.equal(intOpt.get(names[2]), 2);
		tester.equal(intOpt.get(names[3]), -10);
		tester.equal(intOpt.get(names[5]), 545);
		tester.equal(intOpt.get(names[4]), 26);

		tester.equal(intOpt.getSize(), 4);

		// Test that a name with the wrong type is not in the list
		testException(intOpt, &options<int>::get, names[0],
			std::string("Program error! Options.get. Option name 1 does not exists."),
			true,
			std::string("get with name of different type"));
		
		
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();

		delete[] names;
		delete[] values;
		delete[] descriptions;
		delete[] types;
	}

private:
	test tester;

	const int size;

	std::string* names;
	std::string* values;
	std::string* descriptions;
	std::string* types;
		options<float>* floatOpt;

	void testgetSetExistsAddWithString();

	template< class objectClass, class methodReturnType, class arg1type>
	void testException(objectClass& opt,
		methodReturnType (objectClass::*method)(const arg1type&) const,
		arg1type& arg1,
		std::string exceptionMessage,
		bool exceptionFatal,
		std::string name);

	template< class objectClass, class methodReturnType, class arg1type, class arg2type>
	void testException(objectClass& opt,
		methodReturnType (objectClass::*method)(arg1type&, arg2type&),
		arg1type& arg1,
		arg2type& arg2,
		std::string exceptionMessage,
		bool exceptionFatal,
		std::string name);

	template< class objectClass, class methodReturnType,
		class arg1type, class arg2type, class arg3type>
	void testException(objectClass& opt,
		methodReturnType (objectClass::*method)(arg1type&, arg2type, arg3type),
		arg1type& arg1,
		arg2type& arg2,
		arg3type arg3,
		std::string exceptionMessage,
		bool exceptionFatal,
		std::string name);


	template< typename type >
	static type converter(std::string number) {
		type result;
		
		std::stringstream ss(number);
		
		return ss >> result ? result : throw exception(std::string("Not a number: "+number), true);
	}

	static std::string returnInput(std::string number) {
		return number;
	}
};

void testOptions::testgetSetExistsAddWithString() {

	options<std::string> stringOpt(std::string("string"), size, names, 
								values, types, descriptions, &convertString::toString);

	tester.equal(stringOpt.getSize(), 2);
	std::string noOptionWithName = "DOESNOTEXIST";

	// Test get
	tester.equal(stringOpt.get(names[1]), values[1], values[1] + "!" + stringOpt.get(names[1]));
	tester.equal(stringOpt.get(names[0]), values[0]);

	testException(stringOpt, &options<std::string>::get, noOptionWithName,
		std::string("Program error! Options.get. Option DOESNOTEXIST does not exists."),
		true,
		std::string("get"));

	// Test set
	std::string newValue = "new value";
	stringOpt.set(names[0], newValue);
	tester.equal(stringOpt.get(names[0]), newValue, "get after set returns the wrong value");

	try {
		stringOpt.set(noOptionWithName, newValue);
		tester.failedTest("set did not throw exception as expected");
	}
	catch (exception& excep) {
		
		std::string expect = "Program error! Options.set. Option DOESNOTEXIST does not exists.";
		tester.equal(excep.getMessage(), expect, "set threw exception with the wrong message");
		tester.equal(excep.getFatal(), true);
	}
	catch (...) {
		tester.failedTest("set threw the wrong exception");
		throw;
	}

	// Test exists
	tester.isTrue(stringOpt.exists(names[0]), "exists returns false");
	tester.isFalse(stringOpt.exists(noOptionWithName), "exists returns true");
	
	// Add values
	std::string addOption = "add option";
	std::string addValue = "add option value";
	stringOpt.add(addOption, addValue); // No description
	tester.equal(stringOpt.get(addOption), addValue, "Adding");
	tester.equal(stringOpt.getDescription(addOption), std::string(""));
	

	try {
		stringOpt.add(names[0], values[0]);
		tester.failedTest("Adding an already existing option did not throw an exception");
	}
	catch (	exception& excep) {
		
		std::string expect = "Program error! Options.add. Option name 1 already exists.";
		tester.equal(excep.getMessage(), expect, "add threw exception with the wrong message");
		tester.equal(excep.getFatal(), true);
	}
	catch (...) {
		tester.failedTest("add threw the wrong exception");
		throw;
	}

	// Test get description
	std::string addOption2 = "add option 2";
	std::string addValue2 = "add option value 2";
	std::string addDescription2 = "Description 2"; // There is addDescription 1
	stringOpt.add(addOption2, addValue2, addDescription2);
	tester.equal(stringOpt.get(addOption2), addValue2);
	tester.equal(stringOpt.getDescription(addOption2), addDescription2);
	

	testException(stringOpt, &options<std::string>::getDescription, noOptionWithName,
		std::string("Program error! Options.getDescription. Option DOESNOTEXIST does not exists."),
		true,
		std::string("getDescription"));
	
	tester.equal(stringOpt.get(names[0]), newValue, "New value default");
	tester.equal(stringOpt.get(names[1]), values[1], "Not changed value default");

	testException(stringOpt, &options<std::string>::getDefault, noOptionWithName,
		std::string("Program error! Options.getDefault. Option DOESNOTEXIST does not exists."),
		true,
		std::string("getDefault"));
	
	// Make sure that all options are still correct
	tester.equal(stringOpt.get(names[0]), newValue);
	tester.equal(stringOpt.getDescription(names[0]), descriptions[0]);
	tester.equal(stringOpt.getDefault(names[0]), values[0]);

	tester.equal(stringOpt.get(names[1]), values[1]);
	tester.equal(stringOpt.getDescription(names[1]), descriptions[1], "Description");
	tester.equal(stringOpt.getDefault(names[1]), values[1]);

	tester.equal(stringOpt.get(addOption), addValue);
	tester.equal(stringOpt.getDescription(addOption), std::string(""));
	tester.equal(stringOpt.getDefault(addOption), addValue);

	tester.equal(stringOpt.get(addOption2), addValue2);
	tester.equal(stringOpt.getDescription(addOption2), addDescription2);
	tester.equal(stringOpt.getDefault(addOption2), addValue2);

}

template< class objectClass, class methodReturnType, class arg1type>
void testOptions::testException(objectClass& opt,
	methodReturnType (objectClass::*method)(const arg1type&) const,
	arg1type& arg1,
	std::string exceptionMessage,
	bool exceptionFatal,
	std::string name) {

	try {
		(opt.*method)(arg1);
		tester.failedTest(name + ": did not throw an exception");
	}
	catch (	exception& excep) {
		tester.equal(excep.getMessage(), exceptionMessage, name + ": threw an exception with the wrong message:\n" + excep.getMessage() + "\n" + exceptionMessage);
		tester.equal(excep.getFatal(), exceptionFatal, name + ": threw a exception with the wrong fatality");
	}
	catch (...) {
		tester.failedTest(name + ": threw the wrong exception");
		throw;
	}
}

template< class objectClass, class methodReturnType, class arg1type, class arg2type>
void testOptions::testException(objectClass& opt,
	methodReturnType (objectClass::*method)(arg1type&, arg2type&),
	arg1type& arg1,
	arg2type& arg2,
	std::string exceptionMessage,
	bool exceptionFatal,
	std::string name) {

	try {
		(opt.*method)(arg1, arg2);
		tester.failedTest(name + ": did not throw an exception");
	}
	catch (	exception& excep) {
		tester.equal(excep.getMessage(), exceptionMessage, name + ": threw an exception with the wrong message:\n" + excep.getMessage() + "\n" + exceptionMessage);
		tester.equal(excep.getFatal(), exceptionFatal, name + ": threw a exception with the wrong fatality");
	}
	catch (...) {
		tester.failedTest(name + ": threw the wrong exception");
		throw;
	}
}

template< class objectClass, class methodReturnType,
	class arg1type, class arg2type, class arg3type>
void testOptions::testException(objectClass& opt,
	methodReturnType (objectClass::*method)(arg1type&, arg2type, arg3type),
	arg1type& arg1,
	arg2type& arg2,
	arg3type arg3,
	std::string exceptionMessage,
	bool exceptionFatal,
	std::string name) {

	try {
		(opt.*method)(arg1, arg2,arg3);
		tester.failedTest(name + ": did not throw an exception");
	}
	catch (	exception& excep) {
		tester.equal(excep.getMessage(), exceptionMessage, name + ": threw an exception with the wrong message:\n" + excep.getMessage() + "\n" + exceptionMessage);
		tester.equal(excep.getFatal(), exceptionFatal, name + ": threw a exception with the wrong fatality");
	}
	catch (...) {
		tester.failedTest(name + ": threw the wrong exception");
		throw;
	}
}
	
#endif /* TESTOPTIONS */
