#ifndef GLOBALPRUN
#define GLOBALPRUN

#include "foldalign.hxx"
#include "exception.cxx"
#include <string>
#include <iostream>

class globalPruneTable {
public:
  globalPruneTable(const positionType seqLength1, const positionType seqLength2,
		   const lengthType Delta, scoreType affineGapCost, bool NoPrune = false)
    : delta(Delta), noPrune(NoPrune) {

    globalDeltaPrun = new scoreType[2*delta+1];
    if (globalDeltaPrun == 0) {
      std::string error = "Could not allocate the delta pruning table. Most likely cause: Out of memory.";
      throw exception(error, false);
    }

    buildTable(seqLength1, seqLength2, affineGapCost);
  }

  inline scoreType get(const lengthType Wi, const lengthType Wk) const {
    return globalDeltaPrun[Wi - Wk + delta];
  }

  void setNewSeqLengths(const positionType seqLength1, const positionType seqLength2, scoreType affineGapCost, const lengthType Delta) {
    delta = Delta;
    buildTable(seqLength1, seqLength2, affineGapCost);
  }

  ~globalPruneTable() {
     delete[] globalDeltaPrun;
  }

  globalPruneTable& operator=(const globalPruneTable& p) {
    if (this != &p) {
      delete[] globalDeltaPrun;
      delta = p.delta;
      noPrune = p.noPrune;
      globalDeltaPrun = new scoreType[2*delta+1];
      for(lengthType d=0; d < 2*delta+1; d++) {
	       globalDeltaPrun[d] = p.globalDeltaPrun[d];
      }
    }
    return *this;
  };

  globalPruneTable(const globalPruneTable& p) {
    delta = p.delta;
    noPrune = p.noPrune;
    globalDeltaPrun = new scoreType[2*delta+1];
    for(lengthType d=0; d < 2*delta+1; d++) {
      globalDeltaPrun[d] = p.globalDeltaPrun[d];
    }
  }

private:
  scoreType* globalDeltaPrun;
  lengthType delta;
  bool noPrune;

  globalPruneTable(); // No default constructor

  inline void buildTable(const positionType seqLength1, const positionType seqLength2, const scoreType affineGapCost) {

    delete[] globalDeltaPrun;
    globalDeltaPrun = new scoreType[2*delta+1];

    if ( !noPrune ) {

      // Setup the Wi -Wk depended cost (0 in the non global case)
      const positionType deltaLength = seqLength1 > seqLength2 ? seqLength1 - seqLength2 : seqLength2 - seqLength1;

      for(lengthType d = 0; d < delta; d++) {
         // The maximum compensation is the length difference between the
	       // two input sequences.
	       scoreType pscore;
	       if (delta -d < deltaLength) {
           pscore = -(delta  - d)*affineGapCost;
	       }
	       else {
	         pscore = -deltaLength*affineGapCost;
	       }
	       globalDeltaPrun[d] = pscore;
	       globalDeltaPrun[2*delta - d] = pscore;
      }
      globalDeltaPrun[delta] = 0;
    }
    else {
      for(lengthType d = 0; d < 2*delta+1; d++) {
	       globalDeltaPrun[d] = 0;
      }
    }
  }
};

#endif
