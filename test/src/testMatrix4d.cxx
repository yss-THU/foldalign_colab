#ifndef TESTMATRIX4D
#define TESTMATRIX4D

#include "test.cxx"
#include "../../src/matrix4d.cxx"
#include <iostream>
#include <string>

class testMatrix4d {
public:
	testMatrix4d(int& passed, int& ran, int& expected, std::string& messages) {
	
		test tester(46, "matrix4d.cxx test");
		
		matrix4d<int> matrix(1);

		// Test set and get
		matrix.set(0,0,0,0,1);
		
		tester.equal(matrix.get(0,0,0,0), 1);

		// Test that a second matrix does not interfere with the first matrix
		matrix4d<int> matrix2(2);
		tester.notEqual(matrix2.get(0,0,0,0), 1, std::string("This test can in rare cases fail"));
		
		tester.equal(matrix.getSize(), 1);
		tester.equal(matrix2.getSize(), 2);
		
		// Test that set and get works for sizes > 1
		matrix2.set(0,0,0,0,3);
		matrix2.set(0,0,0,1,5);
		matrix2.set(0,0,1,0,6);
		matrix2.set(0,0,1,1,7);
		matrix2.set(0,1,0,0,9);
		matrix2.set(0,1,0,1,10);
		matrix2.set(0,1,1,0,50);
		matrix2.set(0,1,1,1,11);
		matrix2.set(1,0,0,0,13);
		matrix2.set(1,0,0,1,17);
		matrix2.set(1,0,1,0,17);
		matrix2.set(1,0,1,1,20);
		matrix2.set(1,1,0,0,26);
		matrix2.set(1,1,0,1,-2);
		matrix2.set(1,1,1,0,100);
		matrix2.set(1,1,1,1,0);

		tester.equal(matrix2.get(0,0,0,0), 3);
		tester.equal(matrix2.get(0,0,0,1), 5);
		tester.equal(matrix2.get(0,0,1,0), 6);
		tester.equal(matrix2.get(0,0,1,1), 7);
		tester.equal(matrix2.get(0,1,0,0), 9);
		tester.equal(matrix2.get(0,1,0,1), 10);
		tester.equal(matrix2.get(0,1,1,0), 50);
		tester.equal(matrix2.get(0,1,1,1), 11);
		tester.equal(matrix2.get(1,0,0,0), 13);
		tester.equal(matrix2.get(1,0,0,1), 17);
		tester.equal(matrix2.get(1,0,1,0), 17);
		tester.equal(matrix2.get(1,0,1,1), 20);
		tester.equal(matrix2.get(1,1,0,0), 26);
		tester.equal(matrix2.get(1,1,0,1), -2);
		tester.equal(matrix2.get(1,1,1,0), 100);
		tester.equal(matrix2.get(1,1,1,1), 0);
		
		matrix2.init(57);
		tester.equal(matrix2.get(0,0,0,0), 57);
		tester.equal(matrix2.get(0,0,0,1), 57);
		tester.equal(matrix2.get(0,0,1,0), 57);
		tester.equal(matrix2.get(0,0,1,1), 57);
		tester.equal(matrix2.get(0,1,0,0), 57);
		tester.equal(matrix2.get(0,1,0,1), 57);
		tester.equal(matrix2.get(0,1,1,0), 57);
		tester.equal(matrix2.get(0,1,1,1), 57);
		tester.equal(matrix2.get(1,0,0,0), 57);
		tester.equal(matrix2.get(1,0,0,1), 57);
		tester.equal(matrix2.get(1,0,1,0), 57);
		tester.equal(matrix2.get(1,0,1,1), 57);
		tester.equal(matrix2.get(1,1,0,0), 57);
		tester.equal(matrix2.get(1,1,0,1), 57);
		tester.equal(matrix2.get(1,1,1,0), 57);
		tester.equal(matrix2.get(1,1,1,1), 57);
		
		// Test assignment operator
		matrix2 = matrix;
		tester.equal(matrix2.getSize(), 1);
		tester.equal(matrix.get(0,0,0,0), 1);
		
		// Test that it is still two separate matrices
		matrix.set(0,0,0,0,2);
		tester.equal(matrix.get(0,0,0,0), 2);
		tester.equal(matrix2.get(0,0,0,0), 1);
		
		// Test copy constructor
		matrix4d<int> matrix3 = matrix2;
		tester.equal(matrix3.get(0,0,0,0), 1);
		tester.equal(matrix2.get(0,0,0,0), 1);
		
		tester.equal(matrix3.getSize(), 1);
		tester.equal(matrix2.getSize(), 1);
		
		matrix3.set(0,0,0,0,-34);
		tester.equal(matrix3.get(0,0,0,0), -34);
		tester.equal(matrix2.get(0,0,0,0), 1);
		
		
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	}

private:
};


#endif
