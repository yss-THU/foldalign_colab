#ifndef PRUNE
#define PRUNE

#include "foldalign.hxx"
#include "globalPruneTable.cxx"
#include "pruneTable.cxx"

class prune {
public:
  prune(scoreType pruneStart, scoreType linearPrune, lengthType localPruneSize,
	const positionType seqLength1, const positionType seqLength2,
	const lengthType delta, scoreType affineGapCost, bool Global, bool NoPrune = false)
  : global(Global), noPrune(NoPrune),
    localPrune(pruneStart, linearPrune, localPruneSize, noPrune),
    globalPrune(seqLength1, seqLength2, delta, affineGapCost, noPrune) {};

  scoreType get(const lengthType Wi, const lengthType Wk) {

    scoreType score = localPrune.get(Wi) > localPrune.get(Wk) ? localPrune.get(Wi) : localPrune.get(Wk);
    if (global) {
      return score + globalPrune.get(Wi, Wk);
    }
    else {
      return score;
    }
  }

  void updatePruneTable(scoreType pruneStart, scoreType linearPrune, lengthType localPruneSize,
		   const positionType seqLength1, const positionType seqLength2, const scoreType affineGapCost, const lengthType delta) {

    localPrune = pruneTable(pruneStart, linearPrune, localPruneSize, noPrune);
    setNewSeqLengths(seqLength1, seqLength2, affineGapCost, delta);

  }

  lengthType getLength() const {
      return localPrune.getLength();
  }

  void setPruneCoefficient(scoreType value) {
    localPrune.setPruneCoefficient(value);
  }

  void buildLocalPruneTableFromFile(readfile*& file) {
    localPrune.buildTableFromFile(file);
  }

  inline void reSizeTable(const lengthType size) {
    localPrune.reSizeTable(size);
  }

  void setNewSeqLengths(const positionType seqLength1, const positionType seqLength2, const scoreType affineGapCost, lengthType delta) {
    globalPrune.setNewSeqLengths(seqLength1, seqLength2, affineGapCost, delta);
  }

private:
  prune();

  bool global;
  bool noPrune;

  pruneTable localPrune;
  globalPruneTable globalPrune;
};

#endif
