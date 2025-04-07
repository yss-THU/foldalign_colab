#ifndef TESTSCOREMATRIX
#define TESTSCOREMATRIX

#include "test.cxx"
#include "../../src/scorematrix.cxx"
#include <iostream>
#include <string>

class testScorematrix {
public:
	testScorematrix(int& passed, int& ran, int& expected, std::string& messages)
	 : tester(test(136, "scorematrix.cxx test:")), 
	   sm(scorematrix(false)),
	   auxFilePath("test/auxfiles/") {
//	 : tester(test(2010, "scorematrix.cxx test:")), 
	
		/*****************************************************************
		Note: Only basepair substitution score of the basepairs:
		AU, CG, GC, GU, and AU are checked in the testBpSubstitutionMatrix
		function
		*****************************************************************/

		messages += "The scorematrix test is not yet complete. Only parts of the code are tested\n";
		//messages += "The large number of tests which tests the actual value of the scorematrix has been deactivated since the values are about to be updated\n";

		// Test the default (67 tests)
		testDefault();

		// Test reading from file
		testReadingFile();
		
		// 1943 tests
//		testSetGCparameters();
		
		tester.printResult();

		passed += tester.getTestsPassed();
		ran += tester.getTestsRan();
		expected += tester.getTestsExpected();
	};
private:

	test tester;
	scorematrix sm;
	const std::string auxFilePath;

	void testDefault();
	void testReadingFile();

	void testSetGCparameters();
	void testGC_0_2_0_1();
	void testGC_0_2_0_2();
	void testGC_0_3_0_0();
	void testGC_0_3_0_1();
	void testGC_0_3_0_2();
	void testGC_0_3_0_3();
	void testGC_0_4_0_0();
	void testGC_0_4_0_1();
	void testGC_0_4_0_2();
	void testGC_0_4_0_3();
	void testGC_0_4_0_4();
	void testGC_0_5_0_0();
	void testGC_0_5_0_1();
	void testGC_0_5_0_2();
	void testGC_0_5_0_3();
	void testGC_0_5_0_4();
	void testGC_0_5_0_5();
	void testGC_0_6_0_0();
	void testGC_0_6_0_1();
	void testGC_0_6_0_2();
	void testGC_0_6_0_3();
	void testGC_0_6_0_4();
	void testGC_0_6_0_5();
	void testGC_0_6_0_6();
	void testGC_0_7_0_2();
	void testGC_0_7_0_3();
	void testGC_0_7_0_4();
	void testGC_0_7_0_5();
	void testGC_0_7_0_6();
	
	void testSubScale_0_5(const scoreType gap, std::string prefix = "");
	void testSubScale_1(const scoreType gap, std::string prefix = "");
	void testSubScale_2(const scoreType gap, std::string prefix = "");

	void testGapPrune(
		const scoreType gapOpen,
		const scoreType gapElongation,
		const scoreType gapBpOpen,
		const scoreType gapBpElongation,
		scoreType pruneScore);

	void testSingleSubstitutionMatrix(
		const scoreType aa,
		const scoreType ac,
		const scoreType ag,
		const scoreType au,
		const scoreType cc,
		const scoreType cg,
		const scoreType cu,
		const scoreType gg,
		const scoreType gu,
		const scoreType uu,
		const scoreType gap,
		std::string prefix);
								 
	void testBpSubstitutionMatrix(
		const scoreType auau,
		const scoreType cgau,
		const scoreType gcau,
		const scoreType guau,
		const scoreType uaau,
		const scoreType ugau,
		const scoreType cgcg,
		const scoreType gccg,
		const scoreType gucg,
		const scoreType uacg,
		const scoreType ugcg,
		const scoreType gcgc,
		const scoreType gugc,
		const scoreType uagc,
		const scoreType uggc,
		const scoreType gugu,
		const scoreType uagu,
		const scoreType uggu,
		const scoreType uaua,
		const scoreType ugua,
		const scoreType ugug,
		std::string prefix);

	void testBpSubOneBp(
		const int nuc1,
		const int nuc2,
		const scoreType au,
		const scoreType cg,
		const scoreType gc,
		const scoreType gu,
		const scoreType ua,
		const scoreType ug,
		std::string namePrefix="");
};

void testScorematrix::testDefault() {

	// Test the default
	tester.equal(sm.getName(), std::string("default_local"));
	testGapPrune(-110, -55, -220, -110, -400);
	testSubScale_1(-110);
}


void testScorematrix::testReadingFile() {

	std::string filename = auxFilePath + "scorematrix.cxx.aux";
	sm = scorematrix(filename, false);

	sm.checkSize(10);

	// Test the default
	tester.equal(sm.getName(), filename, "filename: " + sm.getName() + " != " + filename);
	testGapPrune(-11, -4, -2020, -7, -40);
	testSubScale_1(-11);
}

void testScorematrix::testSetGCparameters() {

		scorematrix defaultSM(false);

		sm = defaultSM;
		testGC_0_2_0_1();

		sm = defaultSM;
		testGC_0_2_0_2();

		sm = defaultSM;
		testGC_0_3_0_0();

		sm = defaultSM;
		testGC_0_3_0_1();

		sm = defaultSM;
		testGC_0_3_0_2();

		sm = defaultSM;
		testGC_0_3_0_3();

		sm = defaultSM;
		testGC_0_4_0_0();

		sm = defaultSM;
		testGC_0_4_0_1();

		sm = defaultSM;
		testGC_0_4_0_2();

		sm = defaultSM;
		testGC_0_4_0_3();

		sm = defaultSM;
		testGC_0_4_0_4();

		sm = defaultSM;
		testGC_0_5_0_0();

		sm = defaultSM;
		testGC_0_5_0_1();

		sm = defaultSM;
		testGC_0_5_0_2();

		sm = defaultSM;
		testGC_0_5_0_3();

		sm = defaultSM;
		testGC_0_5_0_4();

		sm = defaultSM;
		testGC_0_5_0_5();

		sm = defaultSM;
		testGC_0_6_0_0();

		sm = defaultSM;
		testGC_0_6_0_1();

		sm = defaultSM;
		testGC_0_6_0_2();

		sm = defaultSM;
		testGC_0_6_0_3();

		sm = defaultSM;
		testGC_0_6_0_4();

		sm = defaultSM;
		testGC_0_6_0_5();

		sm = defaultSM;
		testGC_0_6_0_6();

		sm = defaultSM;
		testGC_0_7_0_2();

		sm = defaultSM;
		testGC_0_7_0_3();

		sm = defaultSM;
		testGC_0_7_0_4();

		sm = defaultSM;
		testGC_0_7_0_5();

		sm = defaultSM;
		testGC_0_7_0_6();

}


void testScorematrix::testGC_0_2_0_1() {

	std::string prefix = "0.2 0.1";
	sm.setGCparameters(0.25,0.15, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.2_0.1"), prefix);
	scoreType gapOpen = -110;
	testGapPrune(gapOpen, -55, -220, -110, -411);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_2_0_2() {

	std::string prefix = "0.2 0.2";
	sm.setGCparameters(0.25,0.25, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.2_0.2"), prefix);
	scoreType gapOpen = -55;
	testGapPrune(gapOpen, -28, -110, -55, -282);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_3_0_0() {

	std::string prefix = "0.3 0.0";
	sm.setGCparameters(0.35,0.05, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.3_0.0"), prefix);
	scoreType gapOpen = -55;
	testGapPrune(gapOpen, -28, -110, -55, -247);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_3_0_1() {

	std::string prefix = "0.3 0.1";
	sm.setGCparameters(0.35,0.15, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.3_0.1"), prefix);
	scoreType gapOpen = -55;
	testGapPrune(gapOpen, -28, -110, -55, -192);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_3_0_2() {

	std::string prefix = "0.3 0.2";
	sm.setGCparameters(0.35,0.25, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.3_0.2"), prefix);
	scoreType gapOpen = -55;
	testGapPrune(gapOpen, -28, -110, -55, -236);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_3_0_3() {

	std::string prefix = "0.3 0.3";
	sm.setGCparameters(0.35,0.35, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.3_0.3"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -341);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_4_0_0() {

	std::string prefix = "0.4 0.0";
	sm.setGCparameters(0.45,0.05, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.4_0.0"), prefix);
	scoreType gapOpen = -110;
	testGapPrune(gapOpen, -55, -220, -110, -247);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_4_0_1() {

	std::string prefix = "0.4 0.1";
	sm.setGCparameters(0.45,0.15, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.4_0.1"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -253);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_4_0_2() {

	std::string prefix = "0.4 0.2";
	sm.setGCparameters(0.45,0.25, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.4_0.2"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -200);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_4_0_3() {

	std::string prefix = "0.4 0.3";
	sm.setGCparameters(0.45,0.35, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.4_0.3"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -272);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_4_0_4() {

	std::string prefix = "0.4 0.4";
	sm.setGCparameters(0.45,0.45, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.4_0.4"), prefix);
	scoreType gapOpen = -110;
	testGapPrune(gapOpen, -55, -220, -110, -611);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_5_0_0() {

	std::string prefix = "0.5 0.0";
	sm.setGCparameters(0.55,0.05, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.5_0.0"), prefix);
	scoreType gapOpen = -55;
	testGapPrune(gapOpen, -28, -110, -55, -193);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_5_0_1() {

	std::string prefix = "0.5 0.1";
	sm.setGCparameters(0.55,0.15, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.5_0.1"), prefix);
	scoreType gapOpen = -110;
	testGapPrune(gapOpen, -55, -220, -110, -299);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_5_0_2() {

	std::string prefix = "0.5 0.2";
	sm.setGCparameters(0.55,0.24, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.5_0.2"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -233);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_5_0_3() {

	std::string prefix = "0.5 0.3";
	sm.setGCparameters(0.55,0.35, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.5_0.3"), prefix);
	scoreType gapOpen = -110;
	testGapPrune(gapOpen, -55, -220, -110, -425);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_5_0_4() {

	std::string prefix = "0.5 0.4";
	sm.setGCparameters(0.55,0.45, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.5_0.4"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -510);
	testSubScale_0_5(gapOpen, prefix + " ");
}


void testScorematrix::testGC_0_5_0_5() {

	std::string prefix = "0.5 0.5";
	sm.setGCparameters(0.55,0.55, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.5_0.5"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -1050);
	testSubScale_2(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_6_0_0() {

	std::string prefix = "0.6 0.0";
	sm.setGCparameters(0.65,0.05, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.6_0.0"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -132);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_6_0_1() {

	std::string prefix = "0.6 0.1";
	sm.setGCparameters(0.65,0.15, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.6_0.1"), prefix);
	scoreType gapOpen = -110;
	testGapPrune(gapOpen, -55, -220, -110, -280);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_6_0_2() {

	std::string prefix = "0.6 0.2";
	sm.setGCparameters(0.65,0.25, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.6_0.2"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -332);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_6_0_3() {

	std::string prefix = "0.6 0.3";
	sm.setGCparameters(0.65,0.35, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.6_0.3"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -201);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_6_0_4() {

	std::string prefix = "0.6 0.4";
	sm.setGCparameters(0.65,0.45, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.6_0.4"), prefix);
	scoreType gapOpen = -110;
	testGapPrune(gapOpen, -55, -220, -110, -337);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_6_0_5() {

	std::string prefix = "0.6 0.5";
	sm.setGCparameters(0.65,0.55, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.6_0.5"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -334);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_6_0_6() {

	std::string prefix = "0.6 0.6";
	sm.setGCparameters(0.65,0.65, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.6_0.6"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -286);
	testSubScale_1(gapOpen, prefix + " ");
}


void testScorematrix::testGC_0_7_0_2() {

	std::string prefix = "0.7 0.2";
	sm.setGCparameters(0.75,0.25, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.7_0.2"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -171);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_7_0_3() {

	std::string prefix = "0.7 0.3";
	sm.setGCparameters(0.75,0.35, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.7_0.3"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -420);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_7_0_4() {

	std::string prefix = "0.7 0.4";
	sm.setGCparameters(0.75,0.45, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.7_0.4"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -251);
	testSubScale_0_5(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_7_0_5() {

	std::string prefix = "0.7 0.5";
	sm.setGCparameters(0.75,0.55, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.7_0.5"), prefix);
	scoreType gapOpen = -220;
	testGapPrune(gapOpen, -110, -440, -220, -147);
	testSubScale_1(gapOpen, prefix + " ");
}

void testScorematrix::testGC_0_7_0_6() {

	std::string prefix = "0.7 0.6";
	sm.setGCparameters(0.75,0.65, 0, 0, false);	// Last two parameters was picked randomly
	tester.equal(sm.getName(), std::string("default_local_gc_0.7_0.6"), prefix);
	scoreType gapOpen = -55;
	testGapPrune(gapOpen, -28, -110, -55, -222);
	testSubScale_0_5(gapOpen, prefix + " ");
}



void testScorematrix::testSubScale_0_5(const scoreType gap, std::string prefix) {
		testBpSubstitutionMatrix(21, 6, 10, 0, 4, -7,
		                         27, 10, -2, 10, 2,
		                         28, 4, 5, -3,
		                         17, -3, -10,
		                         24, 1,
		                         16, prefix);
		testSingleSubstitutionMatrix( 10, -11, -9, -10, 6, -13, -8, 5, -10, 7, gap, prefix);
}

void testScorematrix::testSubScale_1(const scoreType gap, std::string prefix) {

	testSingleSubstitutionMatrix(19, -22, -18, -19, 11, -25, -15, 9, -20, 13, gap, prefix);
	testBpSubstitutionMatrix(42, 11, 19,  0,  8, -14, 
	                         53, 20, -4, 19,  4,
							 56,  7,  9, -6,
							 34, -6, -19,
							 48,  2,
							 32, prefix);
}

void testScorematrix::testSubScale_2(const scoreType gap, std::string prefix) {
		testBpSubstitutionMatrix(84, 22, 38, 0, 16, -28,
		                         106, 40, -8, 38, 8,
		                         112, 14, 18, -12,
		                         68, -12, -38,
		                         96, 4,
		                         64, prefix);
		testSingleSubstitutionMatrix(38, -44, -36, -38, 22, -50, -30, 18, -40, 26, gap, prefix);
}

void testScorematrix::testGapPrune(
		const scoreType gapOpen,
		const scoreType gapElongation,
		const scoreType gapBpOpen,
		const scoreType gapBpElongation,
		scoreType pruneScore) {

	tester.equal(sm.getPruneScore(0,0), scoreType(pruneScore), "Prun 0");
	tester.equal(sm.getPruneScore(1,1), scoreType(pruneScore+1), "Prun 1");
	tester.equal(sm.getPruneScore(2,2), scoreType(pruneScore+2), "Prun 2");
	tester.equal(sm.getGapOpen(), gapOpen, "gap open");
	tester.equal(sm.getGap(), scoreType(gapElongation - gapOpen), "gap elongation");
	tester.equal(sm.getBpGapOpen(), gapBpOpen, "bp gap open");
	tester.equal(sm.getBpAffineGap(), scoreType(gapBpElongation - gapBpOpen), "bp gap elongation");
}	

void testScorematrix::testSingleSubstitutionMatrix(
								 const scoreType aa,
								 const scoreType ac,
								 const scoreType ag,
								 const scoreType au,
								 const scoreType cc,
								 const scoreType cg,
								 const scoreType cu,
								 const scoreType gg,
								 const scoreType gu,
								 const scoreType uu,
								 const scoreType gap,
								 std::string prefix) {

	tester.equal(sm.getInit(1,1), aa, prefix + "aa");
	tester.equal(sm.getInit(1,2), ac, prefix + "ac");
	tester.equal(sm.getInit(1,3), ag, prefix + "ag");
	tester.equal(sm.getInit(1,4), au, prefix + "au");
	tester.equal(sm.getInit(1,0), gap, prefix + "a gap");
	tester.equal(sm.getInit(0,1), gap, prefix + "gap a");
	
	tester.equal(sm.getInit(2,1), ac, prefix + "ca");
	tester.equal(sm.getInit(2,2), cc, prefix + "cc");
	tester.equal(sm.getInit(2,3), cg, prefix + "cg");
	tester.equal(sm.getInit(2,4), cu, prefix + "cu");
	tester.equal(sm.getInit(2,0), gap, prefix + "c gap");
	tester.equal(sm.getInit(0,2), gap, prefix + "gap c");
	
	tester.equal(sm.getInit(3,1), ag, prefix + "ga");
	tester.equal(sm.getInit(3,2), cg, prefix + "gc");
	tester.equal(sm.getInit(3,3), gg, prefix + "gg");
	tester.equal(sm.getInit(3,4), gu, prefix + "gu");
	tester.equal(sm.getInit(2,0), gap, prefix + "g gap");
	tester.equal(sm.getInit(0,2), gap, prefix + "gap g");

	tester.equal(sm.getInit(4,1), au, prefix + "ua");
	tester.equal(sm.getInit(4,2), cu, prefix + "uc");
	tester.equal(sm.getInit(4,3), gu, prefix + "ug");
	tester.equal(sm.getInit(4,4), uu, prefix + "uu");
	tester.equal(sm.getInit(2,0), gap, prefix + "u gap");
	tester.equal(sm.getInit(0,2), gap, prefix + "gap u");
}

void testScorematrix::testBpSubstitutionMatrix(
		const scoreType auau,
		const scoreType cgau,
		const scoreType gcau,
		const scoreType guau,
		const scoreType uaau,
		const scoreType ugau,
		const scoreType cgcg,
		const scoreType gccg,
		const scoreType gucg,
		const scoreType uacg,
		const scoreType ugcg,
		const scoreType gcgc,
		const scoreType gugc,
		const scoreType uagc,
		const scoreType uggc,
		const scoreType gugu,
		const scoreType uagu,
		const scoreType uggu,
		const scoreType uaua,
		const scoreType ugua,
		const scoreType ugug,
		std::string prefix) {

	/*****************************************************************
	Note: Only basepair substitution score of the basepairs:
	AU, CG, GC, GU, and AU are checked in this function
	*****************************************************************/

	testBpSubOneBp(1, 4, auau, cgau, gcau, guau, uaau, ugau, prefix + "au");
	testBpSubOneBp(2, 3, cgau, cgcg, gccg, gucg, uacg, ugcg, prefix + "cg");
	testBpSubOneBp(3, 2, gcau, gccg, gcgc, gugc, uagc, uggc, prefix + "gc");
	testBpSubOneBp(3, 4, guau, gucg, gugc, gugu, uagu, uggu, prefix + "gu");
	testBpSubOneBp(4, 1, uaau, uacg, uagc, uagu, uaua, ugua, prefix + "ua");
	testBpSubOneBp(4, 3, ugau, ugcg, uggc, uggu, ugua, ugug, prefix + "ug");

}

void testScorematrix::testBpSubOneBp(
		const int nuc1,
		const int nuc2,
		const scoreType au,
		const scoreType cg,
		const scoreType gc,
		const scoreType gu,
		const scoreType ua,
		const scoreType ug,
		std::string namePrefix) {
		
	tester.equal(sm.getScore(1,4,nuc1,nuc2), au, namePrefix + "au");
	tester.equal(sm.getScore(2,3,nuc1,nuc2), cg, namePrefix + "cg");
	tester.equal(sm.getScore(3,2,nuc1,nuc2), gc, namePrefix + "gc");
	tester.equal(sm.getScore(3,4,nuc1,nuc2), gu, namePrefix + "gu");
	tester.equal(sm.getScore(4,1,nuc1,nuc2), ua, namePrefix + "ua");
	tester.equal(sm.getScore(4,3,nuc1,nuc2), ug, namePrefix + "ug");
}
#endif
