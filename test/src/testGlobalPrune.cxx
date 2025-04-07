#ifndef TESTGLOBALPRUNE
#define TESTGLOBALPRUNE

#include "../../src/globalPruneTable.cxx"
#include "test.cxx"

class testGlobalPrune {
public:
  testGlobalPrune(int& passed, int& ran, int& expected, std::string& messages) : tester(test(0, "globalPrune.cxx test:")) {

    positionType seqLength1 = 20;
    positionType seqLength2 = 23;
    lengthType delta = 10;
    scoreType affineGapCost = 4;
    bool noPrune = false;
    positionType lengthDifference = seqLength1 > seqLength2 ? seqLength1 - seqLength2 : seqLength2 - seqLength1;
    
    globalPruneTable prune(seqLength1, seqLength2, delta, affineGapCost, noPrune);
    testGet(prune, affineGapCost, lengthDifference, delta);

    globalPruneTable prune2 = prune;
    testGet(prune2, affineGapCost, lengthDifference, delta);

    positionType newSeqLength1 = 31;
    positionType newSeqLength2 = 26;
    scoreType newAffineGapCost = 7;
    positionType newLengthDifference = newSeqLength1- newSeqLength2;
    prune2.setNewSeqLengths(newSeqLength1, newSeqLength2, newAffineGapCost, delta);
    testGet(prune2, newAffineGapCost, newLengthDifference, delta); // Test setNewSeqLengths
    testGet(prune, affineGapCost, lengthDifference, delta); // Make sure prune is unchanged by the changes to prune2

    prune = prune2;
    prune2.setNewSeqLengths(seqLength1, seqLength2, affineGapCost, delta);
    testGet(prune, newAffineGapCost, newLengthDifference, delta);
    
    
    // Test noPrune
    globalPruneTable noPruneTest(seqLength1, seqLength2, delta, affineGapCost, true);
    testGetNoPrune(noPruneTest, lengthDifference, delta);
    
    noPruneTest.setNewSeqLengths(newSeqLength1, newSeqLength2, newAffineGapCost, delta);
    testGetNoPrune(noPruneTest, lengthDifference, delta);

    globalPruneTable noPruneTest2 = noPruneTest;
    testGetNoPrune(noPruneTest2, lengthDifference, delta);
    
    noPruneTest = noPruneTest2;
    testGetNoPrune(noPruneTest, lengthDifference, delta);
    
    passed += tester.getTestsPassed();
    ran += tester.getTestsRan();
    expected += tester.getTestsExpected();

    tester.printResult();

  }
private:
  
  test tester;
  
  void testGet(globalPruneTable& prune, scoreType affineGapCost, positionType lengthDifference, lengthType delta) {

    tester.equal(prune.get(0,0), scoreType(0), "0");
    tester.equal(prune.get(0,1), scoreType(-affineGapCost), "1");
    tester.equal(prune.get(lengthDifference-1,0), scoreType(-(lengthDifference-1)*affineGapCost), "lengthDifference-1");
    tester.equal(prune.get(lengthDifference,0), scoreType(-lengthDifference*affineGapCost), "lengthDifference");
    tester.equal(prune.get(lengthDifference+1,0), scoreType(-lengthDifference*affineGapCost), "lengthDifference+1");
    tester.equal(prune.get(delta,0), scoreType(-lengthDifference*affineGapCost), "delta");
    tester.equal(prune.get(1,0), scoreType(-affineGapCost), "1");
    tester.equal(prune.get(0,lengthDifference-1), scoreType(-(lengthDifference-1)*affineGapCost), "lengthDifference-1");
    tester.equal(prune.get(0,lengthDifference), scoreType(-lengthDifference*affineGapCost), "lengthDifference");
    tester.equal(prune.get(0,lengthDifference+1), scoreType(-lengthDifference*affineGapCost), "lengthDifference+1");
    tester.equal(prune.get(0,delta), scoreType(-lengthDifference*affineGapCost), "delta");
    tester.equal(prune.get(7,delta+7), scoreType(-lengthDifference*affineGapCost), "delta no zero Wi & Wk");
  }
  
  void testGetNoPrune (globalPruneTable& noPruneTest, positionType lengthDifference, lengthType delta) {
    tester.equal(noPruneTest.get(0,0), scoreType(0), "0");
    tester.equal(noPruneTest.get(0,1), scoreType(0), "1");
    tester.equal(noPruneTest.get(lengthDifference-1,0), scoreType(0), "lengthDifference-1");
    tester.equal(noPruneTest.get(lengthDifference,0), scoreType(0), "lengthDifference");
    tester.equal(noPruneTest.get(lengthDifference+1,0), scoreType(0), "lengthDifference+1");
    tester.equal(noPruneTest.get(delta,0), scoreType(0), "delta");
    tester.equal(noPruneTest.get(1,0), scoreType(0), "1");
    tester.equal(noPruneTest.get(0,lengthDifference-1), scoreType(0), "lengthDifference-1");
    tester.equal(noPruneTest.get(0,lengthDifference), scoreType(0), "lengthDifference");
    tester.equal(noPruneTest.get(0,lengthDifference+1), scoreType(0), "lengthDifference+1");
    tester.equal(noPruneTest.get(0,delta), scoreType(0), "delta");
    tester.equal(noPruneTest.get(7,delta+7), scoreType(0), "delta no zero Wi & Wk");
  }
};


#endif
