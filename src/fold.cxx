#ifndef FOLD
#define FOLD
#include "foldalign.hxx"
#include "foldK.cxx"
#include "sequence.cxx"
#include "output.cxx"
#include "helper.cxx"
#include "arguments.cxx"
#include "results.cxx"
#include "longTermMemory.cxx"
#include "shortTermMemory.cxx"
#include "exception.cxx"
#include "stack_ssl.cxx"
#include "mbllist.cxx"
#include "longCell.cxx"
#include "constraints.cxx"
#include "jl.cxx"
#include "cell.cxx"
#include "threadMaster.cxx"

extern "C" {
#include <pthread.h>
}


#define FOLD_TEMPLATE class stmcell, class ltmcell, class startCell, \
                      bool global, bool globalPrune, bool realigning, \
					  bool mblrealign, bool multiThreaded


#define FOLD_TEMPLATE_PARAMETERS stmcell, ltmcell, startCell, \
                                 global, globalPrune, \
								 realigning, mblrealign, multiThreaded

template< FOLD_TEMPLATE >
class fold {
public:
  //*************************************
  // The constructor intialize the stm matrix,
  //
  // The subsequence 
  // begin_1 --> end_1
  // in the first sequence is aligned to subsequence
  // begin_2 --> end_2
  // in the second sequence

  fold(const positionType begin_I, const positionType end_I,
       const positionType begin_K, const positionType end_K,
       results& res,
       const sequence* const seq_1, const sequence* const seq_2,
       arguments& arg, 
       scorematrix& score,
       longTermMemory< ltmcell >* const ltm,
       output*& out,
       const bool lastRun,
       constraints* cons = 0,
       longTermMemory< startCell >* const startCoordinates = 0,
       stack_ssl< mbllist, int >* const mblMem = 0);



private:
  //*************************************
  // The global variables


  // The local names for the parameters stored in the argument object
  const lengthType lambda;
  const lengthType delta;
  const lengthType chunk_size;
  const positionType memroof;
  const bool mem_info;
  const lengthType n_threads;
  const bool printSeedConstraints;

	
  // This numbers are added to i to find the range of the k coordinate during
  // global, realigning, mblrealign.
  const positionType k_offSet_i;
  const positionType k_low;
  const positionType k_high;

  shortTermMemory< stmcell >* newSTMandInit(const positionType begin_I, const positionType end_I,
				        const positionType begin_K, const positionType end_K,
				        arguments& arg,
				        scorematrix& s_matrix);
  

  void initLTM(longTermMemory< ltmcell >* const ltm, const positionType begin_I, const positionType end_I,
	       const positionType begin_K, const positionType end_K,
	       longTermMemory< startCell >* const startCoordinates);
};

//******************
// The constructor *
//******************


template< FOLD_TEMPLATE >
fold< FOLD_TEMPLATE_PARAMETERS >::fold(
		     const positionType begin_I, const positionType end_I,
		     const positionType begin_K, const positionType end_K,
		     results& res,
		     const sequence* const seq_1, const sequence* const seq_2,
		     arguments& arg, 
		     scorematrix& s_matrix,
		     longTermMemory< ltmcell >* const ltm,
		     output*& out,
		     const bool lastRun,
		     constraints* cons,
		     longTermMemory< startCell >* const startCoordinates,
		     stack_ssl< mbllist, int >* const mblMem)
: lambda(arg.ltOpt("-max_length")),
  delta(arg.ltOpt("-max_diff")),
  chunk_size(arg.ltOpt("-chunk_size")),
  memroof(arg.ptOpt("memory_roof")),
  mem_info(arg.boolOpt("-memory_info")),
  n_threads(arg.ltOpt("-number_of_processors") + 1), //+1 because one worker thread is always Idle (see Suplementary Material)
  printSeedConstraints(false),//arg.boolOpt("-print_seed_constraints")),
  k_offSet_i(begin_K - begin_I),
  k_low(k_offSet_i - 2*delta),
  k_high(k_offSet_i + 2*delta) {

	shortTermMemory< stmcell >* stm = newSTMandInit(begin_I, end_I, begin_K, end_K, arg, s_matrix);

	initLTM(ltm, begin_I, end_I, begin_K, end_K, startCoordinates);

	// For all coordinates along the I-sequence
	positionType i_start = end_I;
	if (startCoordinates != 0) {
		if (startCoordinates->getIStart() > -1) {
			i_start = startCoordinates->getIStart();
			stm->set_I_Position(i_start);
		}
	}

	long n_cons = 0;
	long n_pruns = 0;

	if (cons != 0) {cons->reset();}

	if (!realigning && !mblrealign && printSeedConstraints) {
		std::cout << "MAX_LENGTH: " << lambda << " Cons." << std::endl;
	}
	
	threadMaster< FOLD_TEMPLATE_PARAMETERS >* threader = 
			new threadMaster< FOLD_TEMPLATE_PARAMETERS >(
				n_threads, begin_I, end_I, begin_K, end_K, res, seq_1, seq_2,
				arg, s_matrix, ltm, stm, out, lastRun, cons, startCoordinates,
				mblMem);

	positionType i = i_start;
	while(i >= begin_I) {

		// The j coordinate (j = i +Wi) must be less then the min_top_I if the
		// alignment is to be expanded in the j direction.
		// min_top_I must have the smallest value of i + lambda or the end of the
		// sequence.
		positionType min_top_I = helper::min(positionType(i + lambda), end_I);

		// Get the nucleotide for this position
		const positionType c_i = seq_1->getPos(i);

		// Set the k range
		const positionType k_start = global || realigning || mblrealign ?
			i + k_low < begin_K ? begin_K :
			i + k_low : begin_K;
		const positionType k_stop	= global || realigning || mblrealign ?
			i + k_high > end_K ? end_K : 
			i+k_high : end_K;
		

//std::cout << "Thread start: i-" << i << " k_start-" << k_start << " k_stop-" << k_stop << std::endl;
		threader->startNewThread(i, k_start, k_stop, min_top_I, c_i, cons,
				n_pruns, n_cons);
		// Ready for the next i
		i--;
	}

	delete threader;
	
//	if (cons != 0) {
//		std::cerr << "n cons: " << n_cons << " n pruns: " << n_pruns << std::endl;
//	}

	delete stm;
//	std::cerr << "It is done" << std::endl;
}

template< FOLD_TEMPLATE >
shortTermMemory< stmcell >* fold< FOLD_TEMPLATE_PARAMETERS >::newSTMandInit(
					const positionType begin_I, const positionType end_I,
				        const positionType begin_K, const positionType end_K,
				        arguments& arg,
				        scorematrix& s_matrix) {

	//***********************************************
	//
	// Allocating memory and initalizing
	//
	// Notice: You must delete stm when it is no longer needed.
	
	std::string error = "This is not an error, ok maybe it is.";
	
	shortTermMemory< stmcell >* stm;
	
	try {
		// Functions in this object handles most of the printing
		error = "Could not allocate memory for out object. Most likely cause: Out of memory.";

		// Allocate the short term memory
		if ( global || realigning || mblrealign ) {
			// Global alignment and realignment
			error = "Could not allocate short term memory needed for backtrack or global alignment.";
			// Set up the short term memory

			// Dimension along the I-sequence
			const positionType i_dimension = 1 + n_threads;
			// Dimension along the K-sequence. Only a narrow band is needed. See
			// comment below where the range of k is defined.
			const positionType k_dimension = 4*delta+3;
			// The Wi dimension. An alignment cannot be longer than this
			const lengthType Wi_dimension = end_I - begin_I +1;

//std::cout << "Allocating STM" << std::endl;
			stm = new shortTermMemory< stmcell >(i_dimension, k_dimension,
						 Wi_dimension, arg, s_matrix,
						 n_threads,
						 (begin_I - begin_K + 2*delta +1));
//std::cout << "STM ready"<< std::endl;
		}
		else {
			// The normal scan (local) alignment
			error = "Could not allocate short term memory. Most likely cause: Out of memory.";

			// Dimension along the I-sequence
			const positionType i_dimension = 1 + n_threads;
			// Dimension along the K-sequence. The full chunk is needed.
			const positionType k_dimension = chunk_size+3;
			// The Wi dimension. An alignment cannot be longer than this
			const lengthType Wi_dimension = lambda+1;

			stm = new shortTermMemory< stmcell >(i_dimension, k_dimension, 
						 Wi_dimension, arg, s_matrix,
						 n_threads);
		}
	}
	catch ( exception ) {throw;}
	catch ( ... ) {throw exception(error, false);}

	if (stm == 0) {throw exception(error, false);}

	// Set the position of the memory matrices
	stm->setPositions(end_I, end_K);
	
	return stm;
}

template< FOLD_TEMPLATE >
void fold< FOLD_TEMPLATE_PARAMETERS >::initLTM(longTermMemory< ltmcell >* const ltm,
					       const positionType begin_I, const positionType end_I,
					       const positionType begin_K, const positionType end_K,
					       longTermMemory< startCell >* const startCoordinates) {


	if (ltm != 0) {
		if ( global || realigning || mblrealign ) {
			ltm->setPositions(begin_I, begin_K);
			if (startCoordinates != 0) {startCoordinates->setPositions(begin_I, begin_K);}
		}
		else {
			ltm->setPositions(end_I+1, end_K);
		}
	}
}

#endif /*FOLD */
