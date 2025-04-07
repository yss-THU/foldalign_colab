#ifndef TESTMATRIX2D
#define TESTMATRIX2D

#include "test.cxx"
#include "../../src/matrix2d.cxx"
#include <iostream>
#include <string>

class testMatrix2d {
public:
	testMatrix2d(int& passed, int& ran, int& expected, std::string& messages) {
	
		test tester(22, "matrix2d.cxx test");
		
		matrix2d<int> matrix(1);

		// Test set and get
		matrix.set(0,0,1);
		
		tester.equal(matrix.get(0,0), 1);

		// Test that a second matrix does not interfere with the first matrix
		matrix2d<int> matrix2(2);
		tester.notEqual(matrix2.get(0,0), 1, std::string("This test can in rare cases fail"));
		
		tester.equal(matrix.getSize(), 1);
		tester.equal(matrix2.getSize(), 2);
		
		// Test that set and get works for sizes > 1
		matrix2.set(0,0,3);
		matrix2.set(0,1,5);
		matrix2.set(1,0,6);
		matrix2.set(1,1,7);

		tester.equal(matrix2.get(0,0), 3);
		tester.equal(matrix2.get(0,1), 5);
		tester.equal(matrix2.get(1,0), 6);
		tester.equal(matrix2.get(1,1), 7);
		
		matrix2.init(57);
		tester.equal(matrix2.get(0,0), 57);
		tester.equal(matrix2.get(0,1), 57);
		tester.equal(matrix2.get(1,0), 57);
		tester.equal(matrix2.get(1,1), 57);
		
		// Test assignment operator
		matrix2 = matrix;
		tester.equal(matrix2.getSize(), 1);
		tester.equal(matrix.get(0,0), 1);
		
		// Test that it is still two separate matrices
		matrix.set(0,0,2);
		tester.equal(matrix.get(0,0), 2);
		tester.equal(matrix2.get(0,0), 1);
		
		// Test copy constructor
		matrix2d<int> matrix3 = matrix2;
		tester.equal(matrix3.get(0,0), 1);
		tester.equal(matrix2.get(0,0), 1);
		
		tester.equal(matrix3.getSize(), 1);
		tester.equal(matrix2.getSize(), 1);
		
		matrix3.set(0,0,-34);
		tester.equal(matrix3.get(0,0), -34);
		tester.equal(matrix2.get(0,0), 1);
		
		
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();

	}

private:
};


#endif
