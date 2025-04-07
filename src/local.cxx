#ifndef LOCAL
#define LOCAL

#include "backtrack.cxx"
#include "fold.cxx"
#include "cell.cxx"
#include "longCell.cxx"
#include "exception.cxx"
#include "constraints.cxx"

#define LOCAL_TEMPLATE class stmcell, class ltmcell, class startCell, \
                      bool global, bool realigning, bool mblrealign

#define LOCAL_TEMPLATE_PARAMETERS stmcell, ltmcell, startCell, \
                                 global, realigning, mblrealign
											
											

template< LOCAL_TEMPLATE >
class local {
public:

  local(positionType begin_1, positionType end_1, 
	positionType begin_2, positionType end_2, 
	const sequence* const seq_1, const sequence* const seq_2,
	arguments& arg, 
	scorematrix& score,
	output*& out,
	constraints*& cons,
	constraints*& btCons) {

    const lengthType lambda = arg.ltOpt("-max_length");
    const lengthType chunk_size = arg.ltOpt("-chunk_size");
    const bool plot_score = arg.boolOpt("-plot_score");
    const bool all_scores = arg.boolOpt("-all_scores");
    const lengthType n_threads = arg.ltOpt("-number_of_processors");
    const bool globalPrune = arg.boolOpt("-use_global_pruning");
    const bool printSeedConstraints = false; //arg.boolOpt("-print_seed_constraints");
    const lengthType len_1 = seq_1 ->getLength(); // Short names for the length of the sequences
    const lengthType len_2 = seq_2 ->getLength();

    std::string error = "Unknown error";

    // The score and coordinates of the very best alignment is stored in r
    // r (result) is initialized to the big negative number.
    results r(big_neg, double(-1000.0), noState, 0, 0, 0, 0);
    
    longTermMemory< longCell >* ltm = 0;
    if ( !arg.boolOpt("-nobranch") ) {

      error = "Could not allocate long term memory";
      try {
	ltm = new longTermMemory< longCell >(chunk_size, lambda+n_threads, len_1, 
					     len_2, false, n_threads);
      }
      catch ( exception ) {throw;}
      catch ( ... ) {throw exception(error, false);}
      if (ltm == 0) {throw exception(error, false);}
			
    }
	
    //*****************************************
    //
    // Output the headings


    // When local scores or all scores are printed an extra header section
    // is printed.
    if (plot_score || all_scores || printSeedConstraints) {
      out->localscorehead();
    }
		
    
    //*******************************************
    //
    // Split sequence two into chuncks of size 2*lambda and call align.
    //
		
    // Stop is the end position for the final chunk
    positionType stop = end_2 - lambda;

    // Make sure the stop point is not placed before the starting point.
    if (stop <= begin_2) {stop = begin_2+1;}

    lengthType grow = chunk_size - lambda;
    if (grow <= 0) {grow=1;}

    bool lastRun = false;
    for(positionType begin = begin_2; begin < stop; begin+=grow) {
      positionType end = begin +chunk_size-1;
      if (end >= end_2) {end = end_2; lastRun = true;}
      try {
	if (n_threads == 1) {
		if (globalPrune) {
			fold<cell, longCell, longCell, false, true, false, false, false>(
					begin_1, end_1,
				    begin, end, r, seq_1, seq_2, arg, score, ltm, out, lastRun,
				    cons);
		}
		else {
			fold<cell, longCell, longCell, false, false, false, false, false>(
					begin_1, end_1,
				    begin, end, r, seq_1, seq_2, arg, score, ltm, out, lastRun,
				    cons);
		}
	}
	else {
		if (globalPrune) {
			fold<cell, longCell, longCell, false, true, false, false, true>(
				begin_1, end_1,
			    begin, end, r, seq_1, seq_2, arg, score, ltm, out, lastRun,
			    cons);
		}
		else {
			fold<cell, longCell, longCell, false, false, false, false, true>(
				begin_1, end_1,
			    begin, end, r, seq_1, seq_2, arg, score, ltm, out, lastRun,
			    cons);
		}
	}
      }
      catch ( exception& exc ) {
	if ( !plot_score ) {out->head();}
	out->saveOutputError(plot_score, exc.getMessage());
	if (ltm != 0) {delete ltm;}
	throw;
      }
      catch ( ... ) {
	error = "Could not align sequences.";
	out->saveOutputError(plot_score);
	if (ltm != 0) {delete ltm;}
	throw exception(error, false);
      }
      
    }


    //******************************************
    //
    // Finish the scan output
    //

    // After print local scores or all scores print an end line. If a structure
    // is to be backtracked the print the standard header.
    if (plot_score || all_scores) {out->plotscoreSep();}
    
    // This longTermMemory is no longer needed. (The backtrack object will
    // make its own)
    if (ltm != 0) {delete ltm; ltm = 0;}

    //******************************************
    //
    // Do the backtrack
    //
		
    try {
      if ( !arg.boolOpt("-no_backtrack") ) {
				
	// Backtrack time
	scoreType topscore = r.getScore();
								
	if (topscore == big_neg) {
	  // No alignment found
					
	  std::cerr << "No structural alignment was found between the sequence: ";
	  if (arg.boolOpt("switch")) {
	    std::cerr << seq_2->getName() << " and " << seq_1->getName() << std::endl;
	  }
	  else {
	    std::cerr << seq_1->getName() << " and " << seq_2->getName() << std::endl;
	  }
	  
	  out->head();
	  out->parameters();
	  out->errorNoLocal();		
	  out->saveOutput(false);
	}
	else {
	  positionType best_i;
	  positionType best_k;
	  lengthType   best_Wi;
	  lengthType   best_Wk;
	  r.getPos(best_i, best_k, best_Wi, best_Wk);
	  scoreType similarityScore = r.getSimilarityScore();
	  scoreType energyScore = r.getEnergyScore();
	  
	  positionType j =  best_i+best_Wi;
	  positionType l =  best_k+best_Wk;

//std::cout << "Found that: best_i-" << best_i << " best_k-" << best_k << " best_Wi-" << best_Wi << " best_Wk-" << best_Wk << " (j-" << j << " l-" << l << ")\n";
	  backtrack(similarityScore, energyScore, best_i, j, best_k, l, seq_1, seq_2, arg,
		    score, out, cons);
	}
      }
    }
    catch ( exception ) {throw;}
    catch (...) {
      error = "A serious error occurred during backtrack of the alignment.";
      throw exception(error, false);
    }
  }
};


# endif
