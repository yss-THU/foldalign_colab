#ifndef TESTPRUNE
#define TESTPRUNE

#include "../../src/prune.cxx"
#include "test.cxx"

class testPrune {
public:
  testPrune(int& passed, int& ran, int& expected, std::string& messages) : tester(test(0, "prune.cxx test:")) {

    scoreType pruneStart = -400;
    scoreType linearPrune = 1;
    lengthType size = 10;
    positionType seqLength1 = 34;
    positionType seqLength2 = 43;
    lengthType delta = 17;
    scoreType affineGapCost = 4;
    bool global = false;
    bool noPrune = false;
    prune localTable(pruneStart, linearPrune, size, seqLength1, seqLength2, delta, affineGapCost, global, noPrune);
    tester.equal(localTable.getLength(), lengthType(size+1), "getLength()");

    testGetLocal(localTable, pruneStart, linearPrune, size);
    
    global = true;
    prune globalTable(pruneStart, linearPrune, size, seqLength1, seqLength2, delta, affineGapCost, global, noPrune);
    tester.equal(globalTable.get(0, 0), pruneStart, "Global get 0");
    tester.equal(globalTable.get(1, 0), scoreType(pruneStart + linearPrune - affineGapCost), "Global get 1,0");
    tester.equal(globalTable.get(0, 1), scoreType(pruneStart + linearPrune - affineGapCost), "Global get 0,1");

    // Test updatePruneTable
    pruneStart = -500;
    linearPrune = 3;
    size = 20;
    seqLength1 = 31;
    seqLength2 = 29;
    delta = 7;
    affineGapCost = 7;
    globalTable.updatePruneTable(pruneStart, linearPrune, size, seqLength1, seqLength2, affineGapCost, delta);
    tester.equal(globalTable.get(0, 0), pruneStart, "Update get 0");
    tester.equal(globalTable.get(1, 0), scoreType(pruneStart + linearPrune - affineGapCost), "Update get 1,0");
    tester.equal(globalTable.get(0, 1), scoreType(pruneStart + linearPrune - affineGapCost), "Update get 0,1");
    tester.equal(globalTable.get(2, 0), scoreType(pruneStart + 2*linearPrune - 2*affineGapCost), "Update get 2,0");
    tester.equal(globalTable.get(0, 2), scoreType(pruneStart + 2*linearPrune - 2*affineGapCost), "Update get 0,2");
    tester.equal(globalTable.get(4, 2), scoreType(pruneStart + 4*linearPrune - 2*affineGapCost), "Update get 4,2");
    tester.equal(globalTable.get(2, 4), scoreType(pruneStart + 4*linearPrune - 2*affineGapCost), "Update get 2,4");
    
    
    passed += tester.getTestsPassed();
    ran += tester.getTestsRan();
    expected += tester.getTestsExpected();

    tester.printResult();

  }
private:

  void testGetLocal(prune& localTable, scoreType pruneStart, scoreType linearPrune, lengthType size);
  
  test tester;
};

inline void testPrune::testGetLocal(prune& localTable, scoreType pruneStart, scoreType linearPrune, lengthType size) {
    tester.equal(localTable.get(0, 0), pruneStart, "Get 0");
    tester.equal(localTable.get(1, 0), scoreType(pruneStart+linearPrune), "Get 1,0"); // Check that Wi==0 and Wi==1 are the same
    tester.equal(localTable.get(0, 1), scoreType(pruneStart+linearPrune), "Get 0,1"); // (1,0) and (0,1) should also be the same
    tester.equal(localTable.get(2, 0), scoreType(pruneStart+2*linearPrune), "Get 2,0"); // Should still return Wi==0 value
    tester.equal(localTable.get(0, 2), scoreType(pruneStart+2*linearPrune), "Get 0,2");
    tester.equal(localTable.get(2, 1), scoreType(pruneStart+2*linearPrune), "Get 2,1"); // Should still return Wi==0 value
    tester.equal(localTable.get(1, 2), scoreType(pruneStart+2*linearPrune), "Get 1,2");
    tester.equal(localTable.get(2, 3), scoreType(pruneStart+3*linearPrune), "Get 2,3"); // Should still return Wi==0 value
    tester.equal(localTable.get(3, 2), scoreType(pruneStart+3*linearPrune), "Get 3,2");
    tester.equal(localTable.get(size, size), scoreType(pruneStart+(size)*linearPrune), "Get 10,10");

}

#endif
