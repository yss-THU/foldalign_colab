#ifndef FOLDK
#define FOLDK
#include <sstream>
#include <unistd.h>
#include "foldalign.hxx"
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
#include "thread.cxx"
#include "runKargs.cxx"
#include "foldThreadHandler.cxx"
#include "stmSubMatrix.cxx"
#include "tables.cxx"
#include "combineSimilarityEnergy.cxx"
//#include "fold.cxx"

extern "C" {
#include <pthread.h>
}


// ENERGY_FUNC defines the parameters of the energy calculation functions.
// It is only used in the hope that it makes the code a little more readable
#define ENERGY_FUNC const positionType i, const positionType k, \
    const lengthType Wi, const lengthType Wk,				\
    const positionType c_i, const positionType c_k,			\
    const positionType c_j, const positionType c_l,			\
    stmcell*& cCell

#define FOLDK_TEMPLATE class stmcell, class ltmcell, class startCell, \
                      bool global, bool globalPrune, bool realigning, bool mblrealign, \
											bool multiThreaded


#define FOLDK_TEMPLATE_PARAMETERS stmcell, ltmcell, startCell, \
                                 global, globalPrune, realigning, mblrealign, multiThreaded


template< FOLDK_TEMPLATE >
class foldK {
public:
  //*************************************
  // The constructor intialize the stm matrix,
  //
  // The subsequence
  // begin_1 --> end_1
  // in the first sequence is aligned to subsequence
  // begin_2 --> end_2
  // in the second sequence

  foldK(const positionType begin_I, const positionType end_I,
				const positionType begin_K, const positionType end_K,
				results& res,
				const sequence one, const sequence two,
				arguments argu,
				scorematrix score,
				shortTermMemory< stmcell >* const short_D,
				longTermMemory< ltmcell >* const long_D,
				output*& outp,
				positionType& k_pos,
				foldThreadHandler*& fThreadHandler,
				const bool lstRun,
				constraints* cons = 0,
				longTermMemory< startCell >* const startCoordinates = 0,
				stack_ssl< mbllist, int >* const mblMem = 0);


  // For a given i make the calculations for k,Wi and Wk.
  void* runK(void* paras);

private:
  //*************************************
  // The global variables

  // The current coordinates
  const positionType curr_begin_I;
  const positionType curr_end_I;
  const positionType curr_begin_K;
  const positionType curr_end_K;


  // The coordinates and score for the best alignment.
  results& r;

  // The output class controls most of the printing
  output* out;

	// The k position
  positionType& k;

	foldThreadHandler*& threadHandler;

  const int thread_number;
  // True if this is the last run used to make sure there is no doubble
  // printing of coordinates when -plot_score is used
  const bool lastRun;

  // This the heading of plotscore lines
  const std::string lshead;

  // The local names for the parameters stored in the argument object
  const bool plot_score;
  const bool all_scores;
  const bool no_backtrack;
  const bool flip;
  const bool nobranch;
  const bool print_all_scores;
  const lengthType lambda;
  const lengthType delta;
  const lengthType min_loop;
  const lengthType chunk_size;
  const bool noprune;
  const bool hpstart;
  const positionType memroof;
  const bool mem_info;
	const bool printSeedConstraints;
	const scoreType min_LS_score;
	const lengthType minExpandLength;

  // Given a Window size and a position can this position be extanded to the
  // left without violating the window size < lambda
  positionType max_bottom_I;
  positionType max_bottom_K;

  // The current minimum of i+lambda, end_I. Ditto K
  positionType min_top_I;
  positionType min_top_K;

  // The sequences
  const sequence seq_1;
  const sequence seq_2;

  // The argument object
  arguments arg;

  // The scorematrix
  scorematrix s_matrix;

  // Some of score from the score matrix
  const scoreType mblNuc;
  const scoreType mblHelix;
  const scoreType affineGapCost;
  const scoreType bpOpenGapCost;
  const scoreType bpAffineGapCost;

  // The memory handler object
  shortTermMemory< stmcell >* const stm;
  stmSubMatrix< stmcell >* subStm;
  longTermMemory< ltmcell >* const ltm;

  // An object which holds start coordinates.
  longTermMemory< startCell >* const startCoord;

  // Stack of branch points. Not needed during local alignment
  stack_ssl< mbllist, int >* const mblm;
  stack_ssl< mbllist, int >* local_mblm;

  // The length of the two sequences
  const positionType length_I;
  const positionType length_K;

  // This numbers are added to i to find the range of the k coordinate during
  // global, realigning, mblrealign.
  const positionType k_offSet_i;
  const positionType k_low;
  const positionType k_high;

  // This is true if the high or low score overflow alert has been triggert
  bool highLowAlert;
	const scoreType constraintPruneScore;
  // Program flow control arrays flow_size and the index states are defined
  // in foldalign.hxx
  typedef void (foldK< FOLDK_TEMPLATE_PARAMETERS >::*flowPtr) (
			   const positionType i, const positionType k,
			   const lengthType Wi, const lengthType Wk,
			   const positionType c_i, const positionType c_k,
			   const positionType c_j, const positionType c_l,
			   stmcell*& cCell);

  flowPtr p2calc_ikWiWk[flow_size];
  flowPtr p2calc_iWi[flow_size];
  flowPtr p2calc_kWk[flow_size];
  flowPtr p2calc_ik[flow_size];
  flowPtr p2calc_WiWk[flow_size];
  flowPtr p2calc_i[flow_size];
  flowPtr p2calc_k[flow_size];
  flowPtr p2calc_Wi[flow_size];
  flowPtr p2calc_Wk[flow_size];
  flowPtr p2calc_mbl[flow_size];

  // The functions in this array corrects the external unpaired sequence scores
  void (foldK< FOLDK_TEMPLATE_PARAMETERS >::*p2end[flow_size]) (
				 scoreType& score,
				 const positionType i, const positionType k,
				 const lengthType Wi, const lengthType Wk,
				 stmcell*& cCell);

  // True for all states which cannot be the left part
  // of a branch point
//  bool right_branch[flow_size];

  // True for states which should be stored.
  // equal to the right_branch during scan but eqaul all states during
  // backtrack or global alignment
//  bool right_store[flow_size];

	const tables<realigning, mblrealign> array;

//	std::map<lengthType, std::map<lengthType, mbllist*> >* mblConstraints;

  //*************************************
  // Functions

	// Init the alignments
	inline void initAlignment(const positionType i, const positionType k,
						constraints* cons,
						const positionType c_i, const positionType c_k,
						lengthType& Wi_start);


	inline void initNucleotidePair(const positionType i, const positionType k,
				const positionType c_i, const positionType c_k);

	inline void initStartCoord(const positionType i, const positionType k,
				lengthType& Wi_start);

	inline void initSeedConstraints(const positionType i, const positionType k,
				lengthType& Wi_start, constraints* cons);

  // Makes sure that the alignment is only expanded in the possible directions
  inline void expandAlignment(positionType i, positionType k,
			      lengthType Wi, lengthType Wk,
			      const positionType c_i, const positionType c_k,
			      const positionType c_j, const positionType c_l,
			      stmcell*& cCell);

  // Calls the different score functions.
  inline void expandAlignment_ikWiWk(positionType i, positionType k,
			      lengthType Wi, lengthType Wk,
			      const positionType c_i, const positionType c_k,
			      const positionType c_j, const positionType c_l,
			      const bool WkWi_m1, const bool WkWi_p1,
			      const bool ki_m1, const bool ki_p1,
			      stmcell*& cCell);

  inline void expandAlignment_iWi(positionType i, positionType k,
				 lengthType Wi, lengthType Wk,
				 const positionType c_i, const positionType c_k,
				 const positionType c_j, const positionType c_l,
				 const bool WkWi_m1, const bool ki_m1,
				 const bool ki_p1, const bool WkWi_p1,
				  stmcell*& cCell);

  inline void expandAlignment_kWk(positionType i, positionType k,
				 lengthType Wi, lengthType Wk,
				 const positionType c_i, const positionType c_k,
				 const positionType c_j, const positionType c_l,
				 const bool WkWi_m1, const bool ki_m1,
				 const bool ki_p1, const bool WkWi_p1,
				 stmcell*& cCell);


  // Some checks on Wimax when it is updated.
  inline void checkWimax(lengthType& Wi_max, const positionType i);


  // In which directions can an alignment be expanded.
  inline int check_coord_local(const positionType i, const positionType k,
			       const lengthType Wi, const lengthType Wk) const;

  // Is there room for one more base-pair in this sequence?
  inline bool extra_bp_coord_IK(const positionType i, const positionType k,
				const lengthType Wi, const lengthType Wk) const;
  inline bool extra_bp_coord_I(const positionType i, const positionType k,
			       const lengthType Wi, const lengthType Wk) const;
  inline bool extra_bp_coord_K(const positionType i, const positionType k,
			       const lengthType Wi, const lengthType Wk) const;


  //**************************************************************************
  // Energy calculation functions
  inline void hp2end(scoreType& score, const positionType i, const positionType k, const lengthType Wi, const lengthType Wk, stmcell*& cCell);
  inline void stem2end(scoreType& score, const positionType i, const positionType k, const lengthType Wi, const lengthType Wk, stmcell*& cCell);
  inline void bulge_ik2end(scoreType& score, const positionType i, const positionType k, const lengthType Wi, const lengthType Wk, stmcell*& cCell);
  inline void bulge_WiWk2end(scoreType& score, const positionType i, const positionType k, const lengthType Wi, const lengthType Wk, stmcell*& cCell);
  inline void ilIK2end(scoreType& score, const positionType i, const positionType k, const lengthType Wi, const lengthType Wk, stmcell*& cCell);

  // Hairpin functions
  template< stateType return_state, lengthType grow1, lengthType grow3 >
  inline void hp_align_one_pos(const positionType i, const positionType k,
			       const lengthType Wi, const lengthType Wk,
			       const int& nuc_I, const int& nuc_K,
			       stmcell*& cCell, const scoreType& affine = 0);
  template< stateType return_state > inline void hp_init_pot_bp_IK( ENERGY_FUNC );
  template< stateType return_state > inline void hp_align_i( ENERGY_FUNC )
  {hp_align_one_pos<return_state, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i), 0, cCell);}
  template< stateType return_state > inline void hp_align_gap_k_i( ENERGY_FUNC )
  {hp_align_one_pos<return_state, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i), 0, cCell, affineGapCost);}
  template< stateType return_state > inline void hp_align_k( ENERGY_FUNC )
  {hp_align_one_pos<return_state, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k), cCell);}
  template< stateType return_state > inline void hp_align_gap_I_k( ENERGY_FUNC )
  {hp_align_one_pos<return_state, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k), cCell, affineGapCost);}
  template< stateType return_state > inline void hp_align_Wi( ENERGY_FUNC )
  {hp_align_one_pos<return_state, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i+Wi), 0, cCell);}
  template< stateType return_state > inline void hp_align_gap_K_Wi( ENERGY_FUNC )
  {hp_align_one_pos<return_state, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i+Wi), 0, cCell, affineGapCost);}
  template< stateType return_state > inline void hp_align_Wk( ENERGY_FUNC )
  {hp_align_one_pos<return_state, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k+Wk), cCell);}
  template< stateType return_state > inline void hp_align_gap_I_Wk( ENERGY_FUNC )
  {hp_align_one_pos<return_state, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k+Wk), cCell, affineGapCost);}
  template< stateType return_state > inline void hp_align_ik( ENERGY_FUNC )
  {hp_align_one_pos<return_state, 1, 1>(i, k, Wi, Wk, seq_1.getPos(i), seq_2.getPos(k), cCell);}
  template< stateType return_state > inline void hp_align_WiWk( ENERGY_FUNC )
  {hp_align_one_pos<return_state, 1, 1>(i, k, Wi, Wk, seq_1.getPos(i+Wi), seq_2.getPos(k+Wk), cCell);}

  // Stem functions
  template< stateType return_state > inline void stemIK2ikWiWk( ENERGY_FUNC );
  template< stateType return_state > inline void stemIK2onepair(
			 const positionType i,
			 const positionType k,
			 const lengthType Wi,
			 const lengthType Wk,
			 const int& nuc_start,
			 const int& nuc_end,
			 const positionType& c_start,
			 const positionType& c_end,
			 stmcell*& cCell,
			 const scoreType& affine = 0);
  template< stateType return_state > inline void stemIK2iWi( ENERGY_FUNC )
  {stemIK2onepair<return_state>(i, k, Wi, Wk, seq_1.getPos(i), seq_1.getPos(i+Wi), c_i, c_j, cCell);}
  template< stateType return_state > inline void stemIstemgkWk2iWi( ENERGY_FUNC )
  {stemIK2onepair<return_state>(i, k, Wi, Wk, seq_1.getPos(i), seq_1.getPos(i+Wi), c_i, c_j, cCell, bpAffineGapCost);}
  template< stateType return_state > inline void stemIK2kWk( ENERGY_FUNC )
  {stemIK2onepair<return_state>(i, k, Wi, Wk, seq_2.getPos(k), seq_2.getPos(k+Wk), c_k, c_l, cCell);}
  template< stateType return_state > inline void stemgiWistemK2kWk( ENERGY_FUNC )
  {stemIK2onepair<return_state>(i, k, Wi, Wk, seq_2.getPos(k), seq_2.getPos(k+Wk), c_k, c_l, cCell, bpAffineGapCost);}

  // Bulge functions
  template< stateType return_state, lengthType grow_I, lengthType grow_K >
  inline void growBulge_ik(const positionType i, const positionType k,
			   const lengthType Wi, const lengthType Wk,
			   const int& nuc_I, const int& nuc_K,
			   stmcell*& cCell, const scoreType affine = 0);
  template< stateType return_state, lengthType grow_I, lengthType grow_K >
  inline void growBulge_WiWk(const positionType i, const positionType k,
			     const lengthType Wi, const lengthType Wk,
			     const int& nuc_I, const int& nuc_K,
			     stmcell*& cCell, const scoreType affine = 0);
  template< stateType return_state > inline void bWiIbWkK2ikWiWkCore(
			  const positionType i, const positionType k,
			  const lengthType Wi, const lengthType Wk,
			  const positionType c_i,
			  const positionType c_k,
			  const positionType c_j,
			  const positionType c_l,
			  stmcell*& cCell,
			  const lengthType& len_I,
			  const lengthType& len_K,
			  const lengthType& len1,
			  const lengthType& len2,
			  const lengthType& len3,
			  const lengthType& len4 );
  template< stateType return_state, lengthType grow_i, lengthType grow_k,
	    lengthType grow_Wi, lengthType grow_Wk >
  inline void bulge2intCore(const positionType i,
			    const positionType k, const lengthType Wi,
			    const lengthType Wk, const int& nuc_I,
			    const int& nuc_K,
			    const lengthType& len_I, const lengthType& len_K,
			    stmcell*& cCell, const scoreType& affine = 0 );
  template< stateType return_state > inline void stemIK2ik( ENERGY_FUNC )
  {growBulge_ik<return_state, 1, 1>(i, k, Wi, Wk, seq_1.getPos(i), seq_2.getPos(k), cCell);}
  template< stateType return_state > inline void stemIK2i( ENERGY_FUNC )
  {growBulge_ik<return_state, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i), 0, cCell);}
  template< stateType return_state > inline void stemIKgap2i( ENERGY_FUNC )
  {growBulge_ik<return_state, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i), 0, cCell, affineGapCost);}
  template< stateType return_state > inline void stemIK2k( ENERGY_FUNC )
  {growBulge_ik<return_state, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k), cCell);}
  template< stateType return_state > inline void stemIKgap2k( ENERGY_FUNC )
  {growBulge_ik<return_state, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k), cCell, affineGapCost);}
  template< stateType return_state > inline void stemIK2WiWk( ENERGY_FUNC )
  {growBulge_WiWk<return_state, 1, 1>(i, k, Wi, Wk, seq_1.getPos(i+Wi), seq_2.getPos(k+Wk), cCell);}
  template< stateType return_state > inline void stemIK2Wi( ENERGY_FUNC )
  {growBulge_WiWk<return_state, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i+Wi), 0, cCell);}
  template< stateType return_state > inline void stemIKgap2Wi( ENERGY_FUNC )
  {growBulge_WiWk<return_state, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i+Wi), 0, cCell, affineGapCost);}
  template< stateType return_state > inline void stemIK2Wk( ENERGY_FUNC )
  {growBulge_WiWk<return_state, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k+Wk), cCell);}
  template< stateType return_state > inline void stemIKgap2Wk( ENERGY_FUNC )
  {growBulge_WiWk<return_state, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k+Wk), cCell, affineGapCost);}
  template< stateType return_state > inline void bWiIbWkK2ikWiWk( ENERGY_FUNC )
  {bWiIbWkK2ikWiWkCore<return_state>(i, k, Wi, Wk, c_i, c_k, c_j, c_l, cCell, cCell->getLength2(), cCell->getLength4(), 0, cCell->getLength2(), 0, cCell->getLength4());}
  template< stateType return_state > inline void biIbkK2ikWiWk( ENERGY_FUNC )
  {bWiIbWkK2ikWiWkCore<return_state>(i, k, Wi, Wk, c_i, c_k, c_j, c_l, cCell, cCell->getLength1(), cCell->getLength3(), cCell->getLength1(), 0, cCell->getLength3(), 0);}
  template< stateType return_state > inline void biIbkK2WiWk( ENERGY_FUNC )
  {bulge2intCore<return_state, 0, 0, 1, 1>(i, k, Wi, Wk, seq_1.getPos(i+Wi), seq_2.getPos(k+Wk), cCell->getLength1(), cCell->getLength3(), cCell);}
  template< stateType return_state > inline void biIbkK2Wk( ENERGY_FUNC )
  {bulge2intCore<return_state, 0, 0, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k+Wk), cCell->getLength1(), cCell->getLength3(), cCell);}
  template< stateType return_state > inline void biIbkKgap2Wk( ENERGY_FUNC )
  {bulge2intCore<return_state, 0, 0, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k+Wk), cCell->getLength1(), cCell->getLength3(), cCell, affineGapCost);}
  template< stateType return_state > inline void biIbkK2Wi( ENERGY_FUNC )
  {bulge2intCore<return_state, 0, 0, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i+Wi), 0, cCell->getLength1(), cCell->getLength3(), cCell);}
  template< stateType return_state > inline void biIbkKgap2Wi( ENERGY_FUNC )
  {bulge2intCore<return_state, 0, 0, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i+Wi), 0, cCell->getLength1(), cCell->getLength3(), cCell, affineGapCost);}
  template< stateType return_state > inline void bWiIbWkK2ik( ENERGY_FUNC )
  {bulge2intCore<return_state, 1, 1, 0, 0>(i, k, Wi, Wk, seq_1.getPos(i), seq_2.getPos(k), cCell->getLength2(), cCell->getLength4(), cCell);}
  template< stateType return_state > inline void bWiIbWkK2i( ENERGY_FUNC )
  {bulge2intCore<return_state, 1, 0, 0, 0>(i, k, Wi, Wk, seq_1.getPos(i), 0, cCell->getLength2(), cCell->getLength4(), cCell);}
  template< stateType return_state > inline void bWiIbWkKgap2i( ENERGY_FUNC )
  {bulge2intCore<return_state, 1, 0, 0, 0>(i, k, Wi, Wk, seq_1.getPos(i), 0, cCell->getLength2(), cCell->getLength4(), cCell, affineGapCost);}
  template< stateType return_state > inline void bWiIbWkK2k( ENERGY_FUNC )
  {bulge2intCore<return_state, 0, 1, 0, 0>(i, k, Wi, Wk, 0, seq_2.getPos(k), cCell->getLength2(), cCell->getLength4(), cCell);}
  template< stateType return_state > inline void bWiIbWkKgap2k( ENERGY_FUNC )
  {bulge2intCore<return_state, 0, 1, 0, 0>(i, k, Wi, Wk, 0, seq_2.getPos(k), cCell->getLength2(), cCell->getLength4(), cCell, affineGapCost);}

  // Internal loop functions
  template< stateType return_state > inline void ilIilK2ikWiWk( ENERGY_FUNC );
  template< stateType return_state, lengthType grow_i, lengthType grow_k,
	    lengthType grow_Wi, lengthType grow_Wk >
  inline void growInterLoop(const positionType i,
			    const positionType k, const lengthType Wi,
			    const lengthType Wk, const int& nuc_I,
			    const int& nuc_K,
			    stmcell*& cCell, const scoreType& affine = 0);
  template< stateType return_state > inline void ilIilK2ik( ENERGY_FUNC )
  {growInterLoop<return_state, 1, 1, 0, 0>(i, k, Wi, Wk, seq_1.getPos(i), seq_2.getPos(k), cCell);}
  template< stateType return_state > inline void ilIilK2i( ENERGY_FUNC )
  {growInterLoop<return_state, 1, 0, 0, 0>(i, k, Wi, Wk, seq_1.getPos(i), 0, cCell);}
  template< stateType return_state > inline void ilIilKgap2i( ENERGY_FUNC )
  {growInterLoop<return_state, 1, 0, 0, 0>(i, k, Wi, Wk, seq_1.getPos(i), 0, cCell, affineGapCost);}
  template< stateType return_state > inline void ilIilK2k( ENERGY_FUNC )
  {growInterLoop<return_state, 0, 1, 0, 0>(i, k, Wi, Wk, 0, seq_2.getPos(k), cCell);}
  template< stateType return_state > inline void ilIilKgap2k( ENERGY_FUNC )
  {growInterLoop<return_state, 0, 1, 0, 0>(i, k, Wi, Wk, 0, seq_2.getPos(k), cCell, affineGapCost);}
  template< stateType return_state > inline void ilIilK2WiWk( ENERGY_FUNC )
  {growInterLoop<return_state, 0, 0, 1, 1>(i, k, Wi, Wk, seq_1.getPos(i+Wi), seq_2.getPos(k+Wk), cCell);}
  template< stateType return_state > inline void ilIilK2Wi( ENERGY_FUNC )
  {growInterLoop<return_state, 0, 0, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i+Wi), 0, cCell);}
  template< stateType return_state > inline void ilIilKgap2Wi( ENERGY_FUNC )
  {growInterLoop<return_state, 0, 0, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i+Wi), 0, cCell, affineGapCost);}
  template< stateType return_state > inline void ilIilK2Wk( ENERGY_FUNC )
  {growInterLoop<return_state, 0, 0, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k+Wk), cCell);}
  template< stateType return_state > inline void ilIilKgap2Wk( ENERGY_FUNC )
  {growInterLoop<return_state, 0, 0, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k+Wk), cCell, affineGapCost);}

  // MBL functions
  template< stateType return_state > inline void mIK2ikWiWk( ENERGY_FUNC );
  template< stateType return_state, lengthType grow_i, lengthType grow_k,
	    lengthType grow_Wi, lengthType grow_Wk >
  inline void mIK2unpaired(const positionType i,
			   const positionType k, const lengthType Wi,
			   const lengthType Wk, const int& nuc_I,
			   const int& nuc_K,
			   stmcell*& cCell, const scoreType& affine = 0);
  template< stateType return_state > inline void mIK2ik( ENERGY_FUNC )
  {mIK2unpaired<return_state, 1, 1, 0, 0>(i, k, Wi, Wk, seq_1.getPos(i), seq_2.getPos(k), cCell);}
  template< stateType return_state > inline void mIK2i( ENERGY_FUNC )
  {mIK2unpaired<return_state, 1, 0, 0, 0>(i, k, Wi, Wk, seq_1.getPos(i), 0, cCell);}
  template< stateType return_state > inline void mIKgap2i( ENERGY_FUNC )
  {mIK2unpaired<return_state, 1, 0, 0, 0>(i, k, Wi, Wk, seq_1.getPos(i), 0, cCell, affineGapCost);}
  template< stateType return_state > inline void mIK2k( ENERGY_FUNC )
  {mIK2unpaired<return_state, 0, 1, 0, 0>(i, k, Wi, Wk, 0, seq_2.getPos(k), cCell);}
  template< stateType return_state > inline void mIKgap2k( ENERGY_FUNC )
  {mIK2unpaired<return_state, 0, 1, 0, 0>(i, k, Wi, Wk, 0, seq_2.getPos(k), cCell, affineGapCost);}
  template< stateType return_state > inline void mIK2WiWk( ENERGY_FUNC )
  {mIK2unpaired<return_state, 0, 0, 1, 1>(i, k, Wi, Wk, seq_1.getPos(i+Wi), seq_2.getPos(k+Wk), cCell);}
  template< stateType return_state > inline void mIK2Wi( ENERGY_FUNC )
  {mIK2unpaired<return_state, 0, 0, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i+Wi), 0, cCell);}
  template< stateType return_state > inline void mIKgap2Wi( ENERGY_FUNC )
  {mIK2unpaired<return_state, 0, 0, 1, 0>(i, k, Wi, Wk, seq_1.getPos(i+Wi), 0, cCell, affineGapCost);}
  template< stateType return_state > inline void mIK2Wk( ENERGY_FUNC )
  {mIK2unpaired<return_state, 0, 0, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k+Wk), cCell);}
  template< stateType return_state > inline void mIKgap2Wk( ENERGY_FUNC )
  {mIK2unpaired<return_state, 0, 0, 0, 1>(i, k, Wi, Wk, 0, seq_2.getPos(k+Wk), cCell, affineGapCost);}

  template< stateType return_state > inline void bWiWkIK2mbl( ENERGY_FUNC );
  template< stateType return_state > inline void stemIK2mbl( ENERGY_FUNC );
  inline void mblCore(const scoreType left_similarity_score,
			  const scoreType left_energy_score,
  			  const positionType i,
		      const positionType k, const lengthType Wi,
		      const lengthType len2, const lengthType len4,
		      const lengthType Wk, const stateType return_state);

  // A method use by the constructor to initialize the flow arrays
  inline void buildFlowArrays();

  // This function should never be called
  inline void errorFunction( ENERGY_FUNC ) {
    if (highLowAlert && cCell->getState() == 0 &&
		(cCell->getSimilarityScore() <= big_neg || cCell->getEnergyScore() <= big_neg)) {
      std::cerr << "Score underflow error detected" << std::endl;
    }
    else {
      std::cerr << "Program error! Error function called with arguments: i = "
		<< i << " k = " << k << " Wi = " << Wi << " Wk = " << Wk << std::flush;
      std::cerr << " Old cell score: " << cCell->getSimilarityScore() << " " << cCell->getEnergyScore() << " Old cell state: " << int(cCell->getState()) << std::flush;
      std::cerr << " Name 1: " << seq_1.getName() << " Name 2: " << seq_2.getName() << std::endl;
    }
  };

  // The no calculation function
  inline void noCalc( ENERGY_FUNC ) {return;}
  inline void noCalc(scoreType& score, const positionType i,
		     const positionType k, const lengthType Wi,
		     const lengthType Wk, stmcell*& cCell) {return;}

  inline void stmstore(positionType i, positionType k,
		       lengthType Wi, lengthType Wk,
		       scoreType similarityScore, scoreType energyScore,
			   stateType state,
		       lengthType len1, lengthType len2,
		       lengthType len3, lengthType len4,
		       mbllist* mblptr);

	inline void printPlotScore(const positionType k, results& r_local);

  inline void dump_all_scores(const positionType& i, const positionType& k,
			      const lengthType& Wi, const lengthType& Wk,
			      const scoreType& similarityScore, const scoreType& energyScore,
				  const stateType& state, stmcell*& cCell);

};

//******************
// The constructor *
//******************


template< FOLDK_TEMPLATE >
foldK< FOLDK_TEMPLATE_PARAMETERS >::foldK(
								const positionType begin_I, const positionType end_I,
					      const positionType begin_K, const positionType end_K,
								results& res,
								const sequence one, const sequence two,
						  	arguments argu,
					     	scorematrix score,
								shortTermMemory< stmcell >* const short_D,
					     	longTermMemory< ltmcell >* const long_D,
						  	output*& outp,
								positionType& k_pos,
								foldThreadHandler*& fThreadHandler,
								const bool lstRun,
								constraints* cons,
								longTermMemory< startCell >* const startCoordinates,
					     	stack_ssl< mbllist, int >* const mblMem)
: curr_begin_I(begin_I),
  curr_end_I(end_I),
  curr_begin_K(begin_K),
  curr_end_K(end_K),
  r(res),
  out(outp),
  k(k_pos),
  threadHandler(fThreadHandler),
  thread_number( multiThreaded ? threadHandler->getThreadNumber() : 0),
  lastRun(lstRun),
  lshead(out->getLShead()),
  plot_score(argu.boolOpt("-plot_score")),
  all_scores(argu.boolOpt("-all_scores")),
  no_backtrack(argu.boolOpt("-no_backtrack")),
  flip(argu.boolOpt("switch")),
  nobranch(argu.boolOpt("-nobranch")),
  print_all_scores(argu.boolOpt("-print_all_LS_scores")),
  lambda(argu.ltOpt("-max_length")),
  delta(argu.ltOpt("-max_diff")),
  min_loop(argu.ltOpt("-min_loop")),
  chunk_size(argu.ltOpt("-chunk_size")),
  noprune(argu.boolOpt("-no_pruning")),
  hpstart(argu.boolOpt("hpstart")),
  memroof(argu.ptOpt("memory_roof")),
  mem_info(argu.boolOpt("-memory_info")),
//  printSeedConstraints(argu.boolOpt("-print_seed_constraints")),
  printSeedConstraints(false),
	min_LS_score(argu.stOpt("-min_LS_score")),
  minExpandLength(mblrealign | global | realigning ? cons == 0 ? 0 : cons->getSeedMinExpandLength(curr_end_I - curr_begin_I+1) : -1),//argu.ltOpt("-min_seed_expand_length")),
  seq_1(one),
  seq_2(two),
  arg(argu),
  s_matrix(score),
  mblNuc(s_matrix.getMblNuc()),
  mblHelix(2*s_matrix.getMblAffine()),
  affineGapCost(s_matrix.getGap()),
  bpOpenGapCost(s_matrix.getBpGapOpen()),
  bpAffineGapCost(s_matrix.getBpAffineGap()),
  stm(short_D),
  ltm(long_D),
  startCoord(startCoordinates),
  mblm(mblMem),
  local_mblm(0),
  length_I(seq_1.getLength()),
  length_K(seq_2.getLength()),
  k_offSet_i(curr_begin_K - curr_begin_I),
  k_low(k_offSet_i - 2*delta),
  k_high(k_offSet_i + 2*delta),
  highLowAlert(false),
  constraintPruneScore(cons == 0 ? s_matrix.getPruneScore(0, 0) - delta*affineGapCost : cons->getMinScore() -100)
{

//	pthread_mutex_init(&lockout, NULL);
//std::cout << "HPSTART: " << hpstart << " " << startCoordinates << " " << cons << std::endl;
  // Set up the pruning table

//std::cout << "Min expand: " << minExpandLength << " mblrealign " << mblrealign << " global " << global << " realign " << realigning << " cons " << cons << std::endl;

	// Set up the flow arrays
	buildFlowArrays();
}


template< FOLDK_TEMPLATE >
void* foldK< FOLDK_TEMPLATE_PARAMETERS >::runK(void* paras) {

	runKargs* k_paras = (runKargs*) paras;

	const positionType i = k_paras->i;
	const positionType k_start = k_paras->k_start;
	const positionType k_stop = k_paras->k_stop;
	min_top_I = k_paras->min_top_I;
	const positionType c_i = k_paras->c_i;
	constraints* cons = k_paras->cons;
//	long n_pruns = k_paras->n_pruns;
//	long n_cons = k_paras->n_cons;

	//if (i == 7) sleep(2); //This line forces an old multithread bug. Debug only. REMOVE ME

	// All allowed coordinates along the K-sequence.
	if (multiThreaded) {
		threadHandler->lockAll();
	}

//std::cerr << "I: " << i << " " << std::endl;
	subStm = new stmSubMatrix< stmcell >(stm, i, k_stop);

	for(k = k_stop; k >= k_start; k--) {

		if (multiThreaded) {
			threadHandler->releasePostThread();
		}

		results r_local;

		subStm->newK(k);
		// Get the nucleotide for this position
		const positionType c_k = seq_2.getPos(k);

		lengthType Wi_start = 0;

		// The i coordinate must be larger than the max_bottom_I value before
		// the alignment can be expanded in the i direction.
		max_bottom_I = curr_begin_I;

		if (mblrealign) {
			// When a mbl score is calculated its information is stored in the
			// local mblm list. When an alignment is found to be a branch point
			// then the local information is moved into the global mblm list.

			local_mblm = new stack_ssl< mbllist, int >;

			if (local_mblm == 0) {
				std::string error = "Could not allocate memory for the local mbl list. Most likely cause: Out of memory";
				throw exception(error, false);
			}
		}

		initAlignment(i, k, cons, c_i, c_k, Wi_start);

		// This is the maximum Wi (i,k given) for which a value has been assign.
		// Get the size of the largest Wi window seen so far.
		lengthType Wi_max = subStm->getWimax(i,k, thread_number);


		// Make sure the top Wi size do not overflow the sequence end, or
		// lambda etc.
		checkWimax(Wi_max, i);

		two_link< ltmcell >* subLtm = 0;
		if (ltm != 0) {
			subLtm = ltm->getLockAndResetSubMatrix(i, k, thread_number);
		}


		// Wi window size along the I-sequence (-1).
		// j = i + Wi. Window size 0 is a one nucleotide long alignment
//std::cout << "LOOP i-" << i << " k-" <<k << " Wi:" << Wi_start << "--" << Wi_max <<  " " << hpstart << std::endl;
		for(lengthType Wi = Wi_start; Wi <= Wi_max; Wi++) {


			// If the window size is equal to lambda then it is not possible to
			// expand in the i direction since this would make the resulting
			// alignment longer than lambda.
			if (Wi == lambda) {max_bottom_I = i;}

			// The nucleotide at the i+Wi position.
			const positionType c_j = seq_1.getPos(i+Wi);

			// The range of the window size along the K-sequence: Wk
			const lengthType Wk_start = Wi - delta > 0 ? Wi - delta : 0; // Plus one because the for loop
			lengthType Wk_end	 = Wi + delta+1;	 // runs from = to < not <=.

			if (Wk_end	 > lambda	) {Wk_end	 = lambda;}
			if (Wk_end + k > curr_end_K +1) {Wk_end = curr_end_K - k+1;}

			// Only alignments with k larger than the max_bottom_K can be expanded
			// in the k direction.
			max_bottom_K = curr_begin_K;

			// An alignment must have an l = k + Wk which is less than the
			// min_top_K to be expanded in the l direction.
			min_top_K = (k + lambda < curr_end_K) ? k + lambda : curr_end_K;

			for(lengthType Wk = Wk_start; Wk < Wk_end; Wk++) {

				// Get the current alignment stmcell. It holds score, state, and the
				// four lengths.

				stmcell* cCell = subStm->getPos(i, k, Wi, Wk);

				// If the current alignment stmcell is empty (pruned away) hurry
				// to the next alignment
				if (cCell == 0) {continue;}

				//threadHandler->lockOutput();
				//std::cout << "CELL: i-" << i << " k-" << k << " Wi-" << Wi << " Wk-" << Wk << " Score:"  << cCell->getSimilarityScore() << "/" << cCell->getEnergyScore() << " State:" << int(cCell->getState()) << std::endl;
				//threadHandler->unlockOutput();

				// If Wk is equal to lambda then the current alignment can not
				// be expanded in the k direction. Since Wk == lambda is the last
				// Wk for a given Wi it is not nessecary to reset the max_bottom_K
				// to begin_K for this Wi. For the next Wi it will happen
				// automatically
				if (Wk == lambda) {max_bottom_K = k;}

				// The nucleotide at the k+Wk position
				const positionType c_l = seq_2.getPos(k+Wk);

				// Get the score and state
				scoreType similarityScore = cCell->getSimilarityScore();
				scoreType energyScore = cCell->getEnergyScore();
				const stateType state = cCell->getState();
				mbllist* mblpointer = 0;

				const stateType seed = array.convertState2seed(state);
				if (!realigning && !mblrealign && printSeedConstraints && seed != noState) {
					if (seed == noState) {
						std::cerr << "Seed error: " << i << " " << i+Wi << " " << k << " " << k+Wk << " " << similarityScore << " " << energyScore << " " << int(seed) << " "  << flip << std::endl;
					}
					if (flip) {
						std::cout << k << " " << k+Wk << " " << i << " " << i+Wi << " " << similarityScore << " " << energyScore << " " << int(seed) << " Cons." << std::endl;
					}
					else {
						std::cout << i << " " << i+Wi << " " << k << " " << k+Wk << " " << similarityScore << " " << energyScore << " " << int(seed) << " Cons." << std::endl;
					}
				}

				if ( (similarityScore < warn_low ||
					  similarityScore > warn_high ||
					  energyScore < warn_low ||
					  energyScore > warn_high ) &&
					  !highLowAlert) {

					highLowAlert = true;
					std::cerr << "Warning. The score has a size where";
					if (similarityScore < warn_low || energyScore < warn_low) {
						std::cerr << " underflow ";
					}
					else {
						std::cerr << " overflow ";
					}
					std::cerr << "compared to the big_neg score is possible. ";
					std::cerr << "big_neg is defined in foldalign.hxx. ";
					std::cerr << "Please change it to a lower number and recompile ";
					std::cerr << "(or email me for help)." << std::endl;
					std::cerr << "(" << i << "," << i+Wi << ") -> (" << k << "," << k+ Wk << ") Score: " << similarityScore << " " << energyScore << " State: " << int(state) << std::endl;
				}

				if (mblrealign) {

					mblpointer = cCell->getPointer();

					if ( state == mblIK) {

						// Get the mbllist info from local_mblm
						mbllist mbl = *cCell->getPointer();
						mblpointer = mblm->push_ptr(mbl);
						cCell->setPointer(mblpointer);
					}
				}

				// If the alignment contains unpaired nucleotides which are not
				// closed by any base-pair (external unpaired nucleotides) the
				// score is recalculated to the score of mbl unpaired nucleotides.
				const stateType endState = array.convertSeed2state(state);
				(this->*p2end[endState])(energyScore, i, k, Wi, Wk, cCell);

				// Store the score and state in the long term memory if
				// the structure can form the right part of a bifurcated
				// structure
				if (ltm != 0) {

					// Only alignments with base-pairs are stored.
					if (array.get_right_store(state)) {

						// Retrive a pointer to a long term storage cell
						ltmcell* store = subLtm->putNext(Wi, Wk);

						// an store the values in it.
						store->set(similarityScore, energyScore, state, mblpointer);
					}
				}

				// This prints the all_scores information. (Mainly a debug
				// option)
				if (all_scores) {
					dump_all_scores(i, k, Wi, Wk, similarityScore, energyScore, state, cCell);
				}

				// The algorithm keeps track of the best local alignment found
				// so far and the best local alignment with coordinates i and k
				// found so far. If the current score is better than the previous
				// best local i & k alignment then store it and check if the score
				// is also better than the score of best alignment found so far.
				scoreType score = combineSimilarityEnergy(similarityScore, energyScore);
				const double loglen = Wi > Wk? s_matrix.getLog(Wi): s_matrix.getLog(Wk);
				const double logScore = double(score)/loglen;
				if ((logScore > r_local.getLogScore() ||
					 (logScore == r_local.getLogScore() && score > r_local.getScore())) &&
					 (state >= min_struc_state)) {
					r_local.store(score, logScore, similarityScore, energyScore, state, i, k, Wi, Wk);
					if ((logScore > r.getLogScore() ||
						 (logScore == r.getLogScore() && score > r.getScore())) &&
						!(global || realigning)) {
						r.store(score, logScore, similarityScore, energyScore, state, i, k, Wi, Wk, mblpointer);
					}
				}

				cCell->setState(array.convertSeed2state(state));

				// Finally it is time to expand the alignment
				expandAlignment(i, k, Wi, Wk, c_i, c_k, c_j, c_l, cCell);
			}

			// There might be a new Wi_max, check it
			Wi_max = subStm->getWimax(i, k);
			max_bottom_I = curr_begin_I;
			checkWimax(Wi_max, i);
		} // End Wi loop

		if (ltm != 0) {
			ltm->unlock(i,k,thread_number);
		}

		if (mblrealign) {
			delete local_mblm;
		}

		if (plot_score) {
			printPlotScore(k, r_local);
		}

		if (multiThreaded) {
			threadHandler->lockAll();
		}

	} // End of the k loop

	delete subStm;
	subStm = 0;

	k = -2; // With k== -2 the post thread can always continue (multithreading)

	if (multiThreaded) {
		threadHandler->finishK();
	}
	stm->clear_old_cells(i);
//std::cerr << "I: " << i << " done" << std::endl;
	return 0;

}


template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::initAlignment(
		const positionType i, const positionType k,
		constraints* cons,
		const positionType c_i, const positionType c_k, lengthType& Wi_start) {

	if (hpstart) {
		initNucleotidePair(i, k, c_i, c_k);
	}

	if (startCoord != 0) {
		initStartCoord(i, k, Wi_start);
	}

	if (cons != 0) {
		initSeedConstraints(i, k, Wi_start, cons);
	}
}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::initNucleotidePair(
				const positionType i, const positionType k,
				const positionType c_i, const positionType c_k) {

	// Currently the algortihm starts by aligning one nucleotide to another.																											//

	stateType init_state = hp_init;

	// The inital score is the alignment score plus the length cost.
	scoreType similarityScore = s_matrix.getInit(c_i, c_k);
	scoreType energyScore = s_matrix.getHpLength(1,1);

	mbllist* empty = 0;
	stmstore(i, k, 0, 0, similarityScore, energyScore, init_state, 1, 0, 1, 0, empty);

}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::initStartCoord(
				const positionType i, const positionType k, lengthType& Wi_start) {

	//****************************************************************//
	//																																//
	// Get the initial alignment from the startCoord object					 //
	//																																//
	//****************************************************************//
	lengthType Wi;
	lengthType Wk;

	startCoord->resetCurrent(i,k);
	startCell* start = startCoord->getNextPos(i, k, Wi, Wk);

	if ( !hpstart ) { Wi_start = Wi; }

	while (start != 0) {
		if (start->getSimilarityScore() > big_neg &&
		   start->getEnergyScore() > big_neg) {
			stmstore(i, k , Wi, Wk, start->getSimilarityScore(),
					 start->getEnergyScore(), start->getState(), 0, 0, 0, 0, 0);
		}

		start = startCoord->getNextPos(i, k, Wi, Wk);
	}
	startCoord->unlock(i,k);
}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::initSeedConstraints(
				const positionType i, const positionType k, lengthType& Wi_start,
				constraints* cons) {

	// Initialises alignments with seed constraints

	jl* start = 0;
	positionType start_size = 0;

	cons->get(i, k, start, start_size);

	for(positionType p = 0; p < start_size; p++) {

		positionType j = start[p].j;
		positionType l = start[p].l;
		lengthType Wi = lengthType(j - i);
		lengthType Wk = lengthType(l - k);
		scoreType similarityScore = start[p].similarityScore;
		scoreType energyScore = start[p].energyScore;
		stateType state = start[p].state;

		if (j > curr_end_I || l > curr_end_K) {
			continue;
		}

		positionType dist = Wi - Wk;
		if (dist > delta || dist < -delta) {
			continue;
		}

		if (state == seed_mblIK && i == curr_begin_I && k == curr_begin_K &&
		    j == curr_end_I && l == curr_end_K) {
			continue;
		}

	    mbllist* mblpointer = 0;

		stmstore(i, k, Wi, Wk, similarityScore, energyScore, state, 0, 0, 0, 0, mblpointer);

		if (Wi_start > j -i ) {Wi_start = j -i ;}

	}

}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::stmstore(positionType i,
							 positionType k,
							 lengthType Wi,
							 lengthType Wk,
							 scoreType similarityScore,
							 scoreType energyScore,
							 stateType state,
							 lengthType len1,
							 lengthType len2,
							 lengthType len3,
							 lengthType len4,
							 mbllist* mblptr) {


	scoreType score = combineSimilarityEnergy(similarityScore, energyScore);

	if (score < s_matrix.getPruneScore(Wi, Wk)) {
	    return;
	}

	// Keep the maximum score
	stmcell* cell_ptr = subStm->putPos(i, k, Wi, Wk);

	scoreType oldScore = combineSimilarityEnergy(cell_ptr->getSimilarityScore(),
								 				 cell_ptr->getEnergyScore());

	if (score > oldScore) {
		cell_ptr->set(similarityScore, energyScore, state, len1, len2, len3, len4, mblptr);
	}
}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::checkWimax(lengthType& Wi_max,
							const positionType i) {

	// Do not overstep the end of the sequence.
	// This could probabaly be optimized
	if (i + Wi_max +1 > curr_end_I) {
		Wi_max = curr_end_I - i;
	}

	if (Wi_max >= lambda) {
		Wi_max = lambda;
	}

}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::printPlotScore(
			const positionType k, results& r_local) {

	// Make sure that each line is only printed once.
	if ((k < curr_end_K - lambda +1) || (lastRun)) {
		if ((r_local.getScore() >= min_LS_score) || print_all_scores) {
			positionType pos_i;
			lengthType pos_j;
			positionType pos_k;
			lengthType pos_l;
			r_local.getPos(pos_i, pos_k, pos_j, pos_l);
			if (flip) {
				helper::swap(pos_i, pos_k);
				helper::swap(pos_j, pos_l);
			}

			if (multiThreaded) {threadHandler->lockOutput();}

			std::cout << lshead << pos_i << " " << pos_i+pos_j;
			std::cout << " "	 << pos_k << " " << pos_k+pos_l;
			std::cout << " "	 << r_local.getScore();
			std::cout << " "   << r_local.getState();
			std::cout << " " << r_local.getSimilarityScore();
			std::cout << " " << r_local.getEnergyScore();
			std::cout << std::endl;

			if (multiThreaded) {threadHandler->unlockOutput();}
		}
	}
}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::dump_all_scores(
		     const positionType& i, const positionType& k,
		     const lengthType& Wi, const lengthType& Wk,
		     const scoreType& similarityScore, const scoreType& energyScore,
			 const stateType& state,
		     stmcell*& cCell) {

	// This prints the all_scores information. (Mainly a debug
	// option)

	const int exp_score = check_coord_local(i, k, Wi, Wk);
	std::ostringstream buf;

	if (!flip) {
		buf << "; AS " << int(i) << " " << int(i + Wi);
		buf << " " << int(k) << " " << int(k+Wk);
		buf << " LTM score " << int(similarityScore) << " " << int(energyScore);
		buf << " STM score " << int(cCell->getSimilarityScore()) << " " << int(cCell->getEnergyScore());
		buf << " State " << int(state) << " Lengths: ";
		buf << int(cCell->getLength1()) << " ";
		buf << int(cCell->getLength2()) << " ";
		buf << int(cCell->getLength3()) << " ";
		buf << int(cCell->getLength4()) << " Exp_score: ";
		buf << exp_score << " Prune_score: ";
		buf << int(s_matrix.getPruneScore(Wi, Wk)) << std::endl;
	}
	else {
		buf << "; AS " << int(k) << " " << int(k + Wk);
		buf << " " << int(i) << " " << int(i+Wi);
		buf << " LTM score " << int(similarityScore) << " " << int(energyScore);
		buf << " STM score " << int(cCell->getSimilarityScore()) << " " << int(cCell->getEnergyScore());
		buf << " State " << int(state) << " Lengths: ";
		buf << int(cCell->getLength3()) << " ";
		buf << int(cCell->getLength4()) << " ";
		buf << int(cCell->getLength1()) << " ";
		buf << int(cCell->getLength2()) << " Exp_score: ";
		buf << exp_score << " Prune_score: ";
		buf << int(s_matrix.getPruneScore(Wi, Wk)) << std::endl;
	}

	if (multiThreaded) {
		threadHandler->lockOutput();
		std::cout << buf.str();
		threadHandler->unlockOutput();
	}
	else
		std::cout << buf.str();
}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::expandAlignment(
												positionType i, positionType k,
			lengthType Wi, lengthType Wk,
			const positionType c_i, const positionType c_k,
			const positionType c_j, const positionType c_l,
			stmcell*& cCell) {

	// Calls a function which expands the input alignment depended on its
	// coordinates.


	// Checks if the alignment can be expanded to a longer alignment in only one
	// of the sequences. ie is this true when Wi or Wk is growes without the
	// other: Wi - delta <= Wk <= Wi + delta
	const bool WkWi_m1 = Wk - Wi -1 >= -delta ? true : false;
	const bool WkWi_p1 = Wk - Wi +1 <=	delta ? true : false;


	const bool ki_m1 = k - 1 < i + k_low && (global || realigning || mblrealign)
		? false : true;
	const bool ki_p1 = k + 1 > i + k_high && (global || realigning || mblrealign)
		? false : true;

	// The state is used to index the energy function.
	const stateType state = cCell->getState();

	// The alignment can be expanded in four directions. These can be combined
	// into 16 different cases. The most insterresting is case 15 which expands in
	// all four direction. Case 0 is the case where expansion in any direction is
	// not possible.

	// exp_score masks which direction the alignment can be expanded.
	const int exp_score = check_coord_local(i, k, Wi, Wk);

	switch (exp_score) {
	case 15:
		// Expand in all directions
		expandAlignment_ikWiWk(i, k, Wi, Wk, c_i, c_k, c_j, c_l, WkWi_m1, WkWi_p1, ki_m1, ki_p1, cCell);
		return;
	case 14:
		// Expand in all directions except i

		// Expand WiWk
		(this->*p2calc_WiWk[state])(i, k, Wi+1, Wk+1, c_i, c_k, c_j, c_l, cCell);

		if (WkWi_m1) {
			// Expand Wi
			(this->*p2calc_Wi[state])(i, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
		}

		expandAlignment_kWk(i, k, Wi, Wk, c_i, c_k, c_j, c_l, WkWi_m1, WkWi_p1, ki_m1, ki_p1, cCell);

		if (! nobranch) {
			// Branching
			(this->*p2calc_mbl[state])(i, k, Wi, Wk, c_i, c_k, c_j, c_l, cCell);
		}
		return;
	case 13:
		// Expand in all directions except k

		// Expand WiWk
		(this->*p2calc_WiWk[state])(i, k, Wi+1, Wk+1, c_i, c_k, c_j, c_l, cCell);

		if (WkWi_p1) {
			// Expand Wk
			(this->*p2calc_Wk[state])(i, k, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
		}

		expandAlignment_iWi(i, k, Wi, Wk, c_i, c_k, c_j, c_l, WkWi_m1, WkWi_p1, ki_m1, ki_p1, cCell);

		if (! nobranch) {
			// Branching
			(this->*p2calc_mbl[state])(i, k, Wi, Wk, c_i, c_k, c_j, c_l, cCell);
		}
		return;
	case 11:
		// Expand in all directions except Wi

		if (Wi >= minExpandLength && Wk >= minExpandLength) {
			// Expand ik
			(this->*p2calc_ik[state])(i-1, k-1, Wi+1, Wk+1, c_i, c_k, c_j, c_l, cCell);

			if (WkWi_m1 && ki_m1) {
				// Expand i
				(this->*p2calc_i[state])(i-1, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
			}
		}

		expandAlignment_kWk(i, k, Wi, Wk, c_i, c_k, c_j, c_l, WkWi_m1, WkWi_p1, ki_m1, ki_p1, cCell);
		return;
	case 7:
		// Expand in all directions except Wk

		if (Wi >= minExpandLength && Wk >= minExpandLength) {

			// Expand ik
			(this->*p2calc_ik[state])(i-1, k-1, Wi+1, Wk+1, c_i, c_k, c_j, c_l, cCell);

			if (WkWi_p1 && ki_p1) {
				// Expand k
				(this->*p2calc_k[state])(i, k-1, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
			}
		}


		expandAlignment_iWi(i, k, Wi, Wk, c_i, c_k, c_j, c_l, WkWi_m1, WkWi_p1, ki_m1, ki_p1, cCell);
		return;
	case 3:
		// Expand in the ik directions

		if (Wi >= minExpandLength && Wk >= minExpandLength) {
			// Unpaired on one side
			// Expand ik
			(this->*p2calc_ik[state])(i-1, k-1, Wi+1, Wk+1, c_i, c_k, c_j, c_l, cCell);

			// One side with gap
			if (WkWi_m1 && ki_m1) {
				// Expand i
				(this->*p2calc_i[state])(i-1, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
			}

			if (WkWi_p1 && ki_p1) {
				// Expand k
				(this->*p2calc_k[state])(i, k-1, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
			}
		}
		return;
	case 5:

		// Expand in the iWi directions
		expandAlignment_iWi(i, k, Wi, Wk, c_i, c_k, c_j, c_l, WkWi_m1, WkWi_p1, ki_m1, ki_p1, cCell);

		return;
	case 9:
		// Expand in the iWk directions

		if (Wi >= minExpandLength && Wk >= minExpandLength) {
			if (WkWi_m1 && ki_m1) {
				// Expand i
				(this->*p2calc_i[state])(i-1, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
			}
		}
		if (WkWi_p1) {
			// Expand Wk
			(this->*p2calc_Wk[state])(i, k, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
		}
		return;
	case 6:
		// Expand in the kWi directions

	 	if (Wi >= minExpandLength && Wk >= minExpandLength) {
			if (WkWi_p1 && ki_p1) {
				// Expand k
				(this->*p2calc_k[state])(i, k-1, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
			}
		}
		if (WkWi_m1) {
			// Expand Wi
			(this->*p2calc_Wi[state])(i, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
		}
		return;
	case 10:
		// Expand in the kWk directions
		expandAlignment_kWk(i, k, Wi, Wk, c_i, c_k, c_j, c_l, WkWi_m1, WkWi_p1, ki_m1, ki_p1, cCell);
		return;
	case 12:
		// Expand in the WiWk directions

		(this->*p2calc_WiWk[state])(i, k, Wi+1, Wk+1, c_i, c_k, c_j, c_l, cCell);
		if (WkWi_m1) {
			(this->*p2calc_Wi[state])(i, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
		}
		if (WkWi_p1) {
			// Expand Wk
			(this->*p2calc_Wk[state])(i, k, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
		}
		// Expand MBL
		if (! nobranch) {
			(this->*p2calc_mbl[state])(i, k, Wi, Wk, c_i, c_k, c_j, c_l, cCell);
		}
		return;
	case 1:
		// Expand i
		if (Wi >= minExpandLength && Wk >= minExpandLength) {
			if (WkWi_m1 && ki_m1) {
				(this->*p2calc_i[state])(i-1, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
			}
		}
		return;
	case 2:
		// Expand k
		if (Wi >= minExpandLength && Wk >= minExpandLength) {
			if (WkWi_p1 && ki_p1) {
				(this->*p2calc_k[state])(i, k-1, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
			}
		}
		return;
	case 4:
		// Expand Wi
		if (WkWi_m1) {
			(this->*p2calc_Wi[state])(i, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
		}
		return;
	case 8:
		// Expand Wk
		if (WkWi_p1) {
			(this->*p2calc_Wk[state])(i, k, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
		}
		return;
	case 0:
		// This is the last stmcell no expansion is possible
		if ( global || mblrealign || realigning ) {

			scoreType similarityScore = cCell->getSimilarityScore();
			scoreType energyScore = cCell->getEnergyScore();
			stateType state = cCell->getState();

			// End correction
			(this->*p2end[state])(energyScore, i, k, Wi, Wk, cCell);

			scoreType score = combineSimilarityEnergy(similarityScore, energyScore);
			r.store(score, big_neg, similarityScore, energyScore, state,i, k, Wi, Wk, cCell->getPointer());
		}
		return;
	default:
		std::string error = "Program error! Illegal expand alignment score found.";
		throw exception(error, false);
	}

}


template< FOLDK_TEMPLATE >
inline int foldK< FOLDK_TEMPLATE_PARAMETERS >::check_coord_local(const positionType i, const positionType k, const lengthType Wi, const lengthType Wk) const {

	// This function determines in which directions the alignment can be expanded
	// It should be used for local alignment (scan alignments) since there are
	// extra constrains during global alignment and realignment.

	int exp_score = 0;
	const int i_score = 1;
	const int k_score = 2;
	const int Wi_score = 4;
	const int Wk_score = 8;
	if ( i > max_bottom_I) {
		exp_score += i_score;
	}
	if ( i + Wi < min_top_I ) {
		exp_score += Wi_score;
	}
	if ( k > max_bottom_K) {
		exp_score += k_score;
	}
	if ( k + Wk < min_top_K ) {
		exp_score += Wk_score;
	}
	return exp_score;
}


template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::expandAlignment_ikWiWk(
				  positionType i, positionType k,
				  lengthType Wi, lengthType Wk,
				  const positionType c_i,
				  const positionType c_k,
				  const positionType c_j,
				  const positionType c_l,
				  const bool WkWi_m1, const bool WkWi_p1,
				  const bool ki_m1, const bool ki_p1,
				  stmcell*& cCell) {

	// Call the functions which expand an alignment in all directions

	stateType state = cCell->getState();

	// True if the next set of nucleotides basepairs otherwise false;
	const bool base_pair_I = s_matrix.getBasepair(seq_1.getPos(i-1), seq_1.getPos(i+Wi+1));
	const bool base_pair_K = s_matrix.getBasepair(seq_2.getPos(k-1), seq_2.getPos(k+Wk+1));

	if (Wi >= minExpandLength && Wk >= minExpandLength) {
		if ( base_pair_I && (Wi >= min_loop) && (Wi < lambda -1) ) {

				// The nucleotides at the new coordinates in the I-sequence base-pairs
				// the hairpin loop has the minimum length and does not exeed lambda

			if ( base_pair_K && (Wk >= min_loop) && (Wk < lambda -1) ) {

				// The K-sequence can also base-pair
				// Expand with base-pair in both sequences

				(this->*p2calc_ikWiWk[state])(i-1, k-1, Wi+2, Wk+2, c_i, c_k, c_j, c_l, cCell);

				if (Wk - Wi < delta -1) {

					// Expand with base-pair in the kWk direction
					(this->*p2calc_kWk[state])(i, k-1, Wi, Wk+2, c_i, c_k, c_j, c_l, cCell);
				}
			}

			if (Wi - Wk < delta -1) {

				// Expand with base-pair in the iWi direction
				(this->*p2calc_iWi[state])(i-1, k, Wi+2, Wk, c_i, c_k, c_j, c_l, cCell);
			}

		} else if ( base_pair_K && (Wk >= min_loop) && (Wk < lambda -1) && (Wk - Wi < delta -1) ) {

			// No I-sequence base-pair. Expand with bp in the K-sequence
			(this->*p2calc_kWk[state])(i, k-1, Wi, Wk+2, c_i, c_k, c_j, c_l, cCell);

		}

		// Unpaired on one side
		// Expand ik
		(this->*p2calc_ik[state])(i-1, k-1, Wi+1, Wk+1, c_i, c_k, c_j, c_l, cCell);

	}

	// Expand jl
	(this->*p2calc_WiWk[state])(i, k, Wi+1, Wk+1, c_i, c_k, c_j, c_l, cCell);

	// One side with gap
	// Expand i
	if (WkWi_m1) {
		if (Wi >= minExpandLength && Wk >= minExpandLength) {
			// Expand i
			if (ki_m1) {
				(this->*p2calc_i[state])(i-1, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
			}
		}
		// Expand Wi
		(this->*p2calc_Wi[state])(i, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
	}

	if (WkWi_p1) {
		if (Wi >= minExpandLength && Wk >= minExpandLength) {
			if (ki_p1) {
				// Expand k
				(this->*p2calc_k[state])(i, k-1, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
			}
		}
		// Expand Wk
		(this->*p2calc_Wk[state])(i, k, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
	}

	// Expand with a multibranched loop

	if (! nobranch ) {
		(this->*p2calc_mbl[state])(i, k, Wi, Wk, c_i, c_k, c_j, c_l, cCell);
	}
}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::expandAlignment_iWi(
				 positionType i, positionType k,
				 lengthType Wi, lengthType Wk,
				 const positionType c_i, const positionType c_k,
				 const positionType c_j, const positionType c_l,
				 const bool WkWi_m1, const bool WkWi_p1,
				 const bool ki_m1, const bool ki_p1,
				 stmcell*& cCell) {
	// Expand in the iWi directions

	const stateType state = cCell->getState();
	const bool base_pair_I = s_matrix.getBasepair(seq_1.getPos(i-1), seq_1.getPos(i+Wi+1));

	if (Wi >= minExpandLength && Wk >= minExpandLength) {
		if (base_pair_I && (Wi >= min_loop) && (Wi < lambda -1) && (Wi - Wk < delta -1) ) {
			// Potential insert bp
			(this->*p2calc_iWi[state])(i-1, k, Wi+2, Wk, c_i, c_k, c_j, c_l, cCell);
		}
	}
	if (WkWi_m1) {
		if (Wi >= minExpandLength && Wk >= minExpandLength) {
			if (ki_m1) {
				(this->*p2calc_i[state])(i-1, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
			}
		}
		(this->*p2calc_Wi[state])(i, k, Wi+1, Wk, c_i, c_k, c_j, c_l, cCell);
	}

}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::expandAlignment_kWk(positionType i, positionType k, lengthType Wi, lengthType Wk, const positionType c_i, const positionType c_k, const positionType c_j, const positionType c_l, const bool WkWi_m1, const bool WkWi_p1, const bool ki_m1, const bool ki_p1, stmcell*& cCell) {

	const stateType state = cCell->getState();
	const bool base_pair_K = s_matrix.getBasepair(seq_2.getPos(k-1), seq_2.getPos(k+Wk+1));

	if (Wi >= minExpandLength && Wk >= minExpandLength) {
		if (base_pair_K && (Wk >= min_loop) && (Wk < lambda -1) && (Wk - Wi < delta -1) ) {
			// Expand kl
			(this->*p2calc_kWk[state])(i, k-1, Wi, Wk+2, c_i, c_k, c_j, c_l, cCell);
		}
	}
	if (WkWi_p1) {
		if (Wi >= minExpandLength && Wk >= minExpandLength) {
			if (ki_p1) {
				(this->*p2calc_k[state])(i, k-1, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
			}
		}
		(this->*p2calc_Wk[state])(i, k, Wi, Wk+1, c_i, c_k, c_j, c_l, cCell);
	}
}


/******************************************************************************
*
* These functions check if the alignment can expand with an extra bp
*
*/


template< FOLDK_TEMPLATE >
inline bool foldK< FOLDK_TEMPLATE_PARAMETERS >::extra_bp_coord_IK(const positionType i, const positionType k, const lengthType Wi, const lengthType Wk) const {

	// This function returns true if an alignment with the i,k, Wi,Wk coordinates
	// can be expanded into an alignment i-1, k-1, Wi+2, Wk+2 with a basepair

	// Using the correct borders (max_bottom and min_top) leeds to trouble during
	// backtrack. An "easy" solution has been chosen to solve the problem.
	// Instead of using the actual borders the sequence borders are used (position
	// 1 and length). This allows the algorithm to make structures where there is
	// a single base-pair at the end of a structure. The backtrack algorithm
	// should be able to handle this, but it is still a bit messy.
	// The old "correct" code has been out-commented but is left below.

//	if ((i > max_bottom_I) && (k > max_bottom_K) &&
//	    (Wi < min_top_I) && (Wk < min_top_K) &&
//		 (Wi +2 <= lambda) && (Wk +2 <= lambda)) {
//		return s_matrix.getBasepair(seq_1.getPos(i-1), seq_1.getPos(i + Wi +1)) &&
//		       s_matrix.getBasepair(seq_2.getPos(k-1), seq_2.getPos(k + Wk +1));
//	}

	if ((i > 1) && (k > 1) &&
	    (i+Wi < length_I) && (k+Wk < length_K) &&
		 (Wi +2 <= lambda) && (Wk +2 <= lambda)) {
		return s_matrix.getBasepair(seq_1.getPos(i-1), seq_1.getPos(i + Wi +1)) &&
		       s_matrix.getBasepair(seq_2.getPos(k-1), seq_2.getPos(k + Wk +1));
	}
	return false;
}


template< FOLDK_TEMPLATE >
inline bool foldK< FOLDK_TEMPLATE_PARAMETERS >::extra_bp_coord_I(const positionType i, const positionType k, const lengthType Wi, const lengthType Wk) const {

	// Can the alignment be expanded with a base-pair in the I-sequence?

	// See important comments in the extra_bp_coord_IK function.

//	if ((i > max_bottom_I) && (Wi < min_top_I) && (Wi +2 <= lambda)) {
	if ((i > 1) && (i+Wi < length_I) && (Wi +2 <= lambda)) {
		return s_matrix.getBasepair(seq_1.getPos(i-1), seq_1.getPos(i + Wi +1));
	}

	return false;
}


template< FOLDK_TEMPLATE >
inline bool foldK< FOLDK_TEMPLATE_PARAMETERS >::extra_bp_coord_K(const positionType i, const positionType k, const lengthType Wi, const lengthType Wk) const {

	// Can the alignment be expanded with a base-pair in the K-sequence?

	// See important comments in the extra_bp_coord_IK function.

//	if ((k > max_bottom_K) && (Wk < min_top_K) && (Wk +2 <= lambda)) {
	if ((k > 1) && (k+Wk < length_K) && (Wk +2 <= lambda)) {
		return s_matrix.getBasepair(seq_2.getPos(k-1), seq_2.getPos(k + Wk +1));
	}

	return false;
}


/******************************************************************************
*
* Energy functions
*
*/

/*********************************
*
* The end functions.
* Recalculates the scores of unpaired nucleotides which are not inclosed by a
* base-pair. ie unpaired nucleotides at the ends of the alignment. The score
* is recalculated from hairpin, bulge, or internal loop to mbl unpaired state.
*/


template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end(scoreType& score,
                                  const positionType i, const positionType k,
											 const lengthType Wi, const lengthType Wk,
											 stmcell*& cCell) {

	// Recalculate the score from hairpin state to outside state ie mbl state.

	lengthType len1 = cCell->getLength1();
	lengthType len3 = cCell->getLength3();

	// Remove the hairpin length score
	score -= s_matrix.getHpLength(len1, len3);

	// Add the mblNuc score
	score += (len1 + len3)*mblNuc;

}

template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::stem2end(scoreType& score,
                                  const positionType i, const positionType k,
											 const lengthType Wi, const lengthType Wk,
											 stmcell*& cCell) {

	// Add the nonGC end score
	score += s_matrix.getNonGCEnd(seq_1.getPos(i),
											seq_1.getPos(i+Wi));
	score += s_matrix.getNonGCEnd(seq_2.getPos(k),
											seq_2.getPos(k+Wk));


}


template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_ik2end(scoreType& score,
                                  const positionType i, const positionType k,
											 const lengthType Wi, const lengthType Wk,
											 stmcell*& cCell) {

	// Recalculate the score from hairpin state to outside state ie mbl state.

	lengthType len1 = cCell->getLength1();
	lengthType len3 = cCell->getLength3();

	// Remove the hairpin length score
	score -= s_matrix.getBulgeLength(len1, len3);

	// Add the mblNuc score
	score += (len1 + len3)*mblNuc;

	// Add the nonGC end score
	score += s_matrix.getNonGCEnd(seq_1.getPos(i+len1),
											seq_1.getPos(i+Wi));
	score += s_matrix.getNonGCEnd(seq_2.getPos(k+len3),
											seq_2.getPos(k+Wk));
}


template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_WiWk2end(scoreType& score,
                                  const positionType i, const positionType k,
											 const lengthType Wi, const lengthType Wk,
											 stmcell*& cCell) {

	// Recalculate the score from hairpin state to outside state ie mbl state.

	lengthType len2 = cCell->getLength2();
	lengthType len4 = cCell->getLength4();

	// Remove the hairpin length score
	score -= s_matrix.getBulgeLength(len2, len4);

	// Add the mblNuc score
	score += (len2 + len4)*mblNuc;

	// Add the nonGC end score
	score += s_matrix.getNonGCEnd(seq_1.getPos(i),
											seq_1.getPos(i+Wi-len2));
	score += s_matrix.getNonGCEnd(seq_2.getPos(k),
											seq_2.getPos(k+Wk-len4));
}


template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end(scoreType& score,
                                  const positionType i, const positionType k,
											 const lengthType Wi, const lengthType Wk,
											 stmcell*& cCell) {

	// Recalculate an internal loop state score to an outside (mbl) state score

	lengthType len1 = cCell->getLength1();
	lengthType len2 = cCell->getLength2();
	lengthType len3 = cCell->getLength3();
	lengthType len4 = cCell->getLength4();

	score -= s_matrix.getIntLoopLength(len1, len2);
	score -= s_matrix.getIntLoopLength(len3, len4);

	score += (len1 + len2 + len3 + len4) * mblNuc;

	// Add the nonGC end score
	score += s_matrix.getNonGCEnd(seq_1.getPos(i+len1),
											seq_1.getPos(i+Wi-len2));
	score += s_matrix.getNonGCEnd(seq_2.getPos(k+len3),
											seq_2.getPos(k+Wk-len4));
}

//=========================================================================
// Hairpin functions


template< FOLDK_TEMPLATE >
template<stateType return_state, lengthType grow1, lengthType grow3>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::hp_align_one_pos(const positionType i, const positionType k,
                                   const lengthType Wi, const lengthType Wk,
											  const int& nuc_I, const int& nuc_K,
											  stmcell*& cCell, const scoreType& affine) {

	// Hairpin to pot hp end
	scoreType similarityScore = cCell->getSimilarityScore();
	similarityScore += s_matrix.getInit(nuc_I, nuc_K);
	similarityScore += affine;

	// Get the lengths
	lengthType len1 = cCell->getLength1() + grow1;
	lengthType len3 = cCell->getLength3() + grow3;

	// Add the new length cost and subtract the old length cost
	scoreType energyScore = cCell->getEnergyScore() +
							s_matrix.getHpLength(len1, len3) -
							s_matrix.getHpLength(len1-grow1, len3-grow3);

	// Store
	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state, len1, 0, len3, 0, cCell->getPointer());

}


template< FOLDK_TEMPLATE >
template<stateType return_state>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::hp_init_pot_bp_IK( ENERGY_FUNC ) {

	// If the next set of nucleotides do not base-pair then it is not necessary
	// to do this calculation

	if (!extra_bp_coord_IK(i, k, Wi, Wk)) {return;}

	// Hairpin to pot hp end
	scoreType similarityScore = cCell->getSimilarityScore();

	similarityScore += s_matrix.getScore(seq_1.getPos(i), seq_1.getPos(i+Wi),
	                           seq_2.getPos(k), seq_2.getPos(k+Wk));

	const lengthType len1 = cCell->getLength1();
	const lengthType len3 = cCell->getLength3();

	// Add the closing stacking score if both loop lengths are above three
	scoreType energyScore = cCell->getEnergyScore();
	if ((len1 > 3) && (len3 > 3)) {
		energyScore += s_matrix.getHpClose(seq_1.getPos(i), seq_1.getPos(i+Wi), c_i, c_j);
		energyScore += s_matrix.getHpClose(seq_2.getPos(k), seq_2.getPos(k+Wk), c_k, c_l);
	}

	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state, len1, 0, len3, 0, cCell->getPointer());

}

//=============================================================================
// Stem functions


template< FOLDK_TEMPLATE >
template<stateType return_state>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::stemIK2ikWiWk( ENERGY_FUNC ) {

	// Stem to stem

	// Add the base-pair sub scores
	scoreType similarityScore = cCell->getSimilarityScore();
	similarityScore += s_matrix.getScore(seq_1.getPos(i), seq_1.getPos(i+Wi), seq_2.getPos(k), seq_2.getPos(k+Wk));

	// Add the stackings
	scoreType energyScore = cCell->getEnergyScore();
	energyScore += s_matrix.getStack(seq_1.getPos(i), seq_1.getPos(i+Wi), c_i, c_j);
	energyScore += s_matrix.getStack(seq_2.getPos(k), seq_2.getPos(k+Wk), c_k, c_l);

	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state, 0, 0, 0, 0, cCell->getPointer());

}


template< FOLDK_TEMPLATE >
template<stateType return_state>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::stemIK2onepair(const positionType i, const positionType k, const lengthType Wi, const lengthType Wk, const int& nuc_start, const int& nuc_end, const positionType& c_start, const positionType& c_end, stmcell*& cCell, const scoreType& affine) {

	// Stem to stem_I_stem_gap_kWk

	scoreType similarityScore = cCell->getSimilarityScore();
	similarityScore += bpOpenGapCost + affine;

	// Add the stackings
	scoreType energyScore = cCell->getEnergyScore();
	energyScore += s_matrix.getStack(nuc_start, nuc_end, c_start, c_end);

	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state, 0, 0, 0, 0, cCell->getPointer());

}



//=============================================================================
// Bulge functions


template< FOLDK_TEMPLATE >
template<stateType return_state, lengthType grow_I, lengthType grow_K>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::growBulge_ik(const positionType i, const positionType k, const lengthType Wi, const lengthType Wk, const int& nuc_I, const int& nuc_K, stmcell*& cCell, const scoreType affine) {

	// Substitution score
	scoreType similarityScore = cCell->getSimilarityScore();
	similarityScore += s_matrix.getInit(nuc_I, nuc_K) + affine;

	lengthType len1 = cCell->getLength1() + grow_I;
	lengthType len3 = cCell->getLength3() + grow_K;

	// Subtract the old length cost
	// Add the new length cost
	scoreType energyScore = cCell->getEnergyScore()
							- s_matrix.getBulgeLength(len1-grow_I, len3-grow_K)
							+ s_matrix.getBulgeLength(len1, len3);

	// Store (The ones are due to the potential unpaired nuc. in the bp)
	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state, len1, 0, len3, 0, cCell->getPointer());
}


template< FOLDK_TEMPLATE >
template<stateType return_state, lengthType grow_I, lengthType grow_K>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::growBulge_WiWk(const positionType i, const positionType k, const lengthType Wi, const lengthType Wk, const int& nuc_I, const int& nuc_K, stmcell*& cCell, const scoreType affine) {

	// Substitution score
	scoreType similarityScore = cCell->getSimilarityScore()
								+ s_matrix.getInit(nuc_I, nuc_K) + affine;


	// Energy score
	lengthType len2 = cCell->getLength2() + grow_I;
	lengthType len4 = cCell->getLength4() + grow_K;

	// Subtract the old length cost
	// Add the new length cost
	scoreType energyScore = cCell->getEnergyScore()
							- s_matrix.getBulgeLength(len2-grow_I, len4-grow_K)
							+ s_matrix.getBulgeLength(len2, len4);

	// Store (The ones are due to the potential unpaired nuc. in the bp)
	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state, 0, len2, 0, len4, cCell->getPointer());
}


template< FOLDK_TEMPLATE >
template<stateType return_state>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::bWiIbWkK2ikWiWkCore(const positionType i, const positionType k, const lengthType Wi, const lengthType Wk, const positionType c_i, const positionType c_k, const positionType c_j, const positionType c_l, stmcell*& cCell, const lengthType& len_I, const lengthType& len_K, const lengthType& len1, const lengthType& len2, const lengthType& len3, const lengthType& len4 ) {

	// If the next set of nucleotides do not base-pair then it is not necessary
	// to do this calculation
	if (!extra_bp_coord_IK(i, k, Wi, Wk)) {return;}

	// This alignment can only be expanded into the stem state

	// Pot hp end to stem
	// The potential end is scored as hairpin here, later it will be recalculated
	// to base-pair (if the next nucleotides base-pairs)
	scoreType similarityScore = cCell->getSimilarityScore()
			+ s_matrix.getScore(seq_1.getPos(i), seq_1.getPos(i+Wi),
	                           seq_2.getPos(k), seq_2.getPos(k+Wk));


	// The positions of the opening bp nucleotides
	const positionType po_i = i + len1 + 1; // +1 since it is one step back
	const positionType po_k = k + len3 + 1;
	const positionType po_j = i + Wi - len2 - 1;
	const positionType po_l = k + Wk - len4 - 1;
	// The nucleotides at the base-pair which opens the loop
	const positionType o_i = seq_1.getPos(po_i);
	const positionType o_k = seq_2.getPos(po_k);
	const positionType o_j = seq_1.getPos(po_j);
	const positionType o_l = seq_2.getPos(po_l);

	scoreType energyScore = cCell->getEnergyScore();
	if ((len_I > 1) || (len_K > 1)) {
		// Either of the bulges are longer than one no stacking is allowed
		// Add non-GC cost

		// The closing bp
		energyScore += s_matrix.getNonGCEnd(seq_1.getPos(i), seq_1.getPos(i+Wi));
		energyScore += s_matrix.getNonGCEnd(seq_2.getPos(k), seq_2.getPos(k+Wk));

		// The opening bp
		energyScore += s_matrix.getNonGCEnd(o_i, o_j);
		energyScore += s_matrix.getNonGCEnd(o_k, o_l);

	}
	else {
		// Stacking across bulge is possible
		energyScore += s_matrix.getStack(seq_1.getPos(i), seq_1.getPos(i+Wi), o_i, o_j);
		energyScore += s_matrix.getStack(seq_2.getPos(k), seq_2.getPos(k+Wk), o_k, o_l);

	}

	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state, 0, 0, 0, 0, cCell->getPointer());
}


template< FOLDK_TEMPLATE >
template<stateType return_state, lengthType grow_i, lengthType grow_k, lengthType grow_Wi, lengthType grow_Wk>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge2intCore(const positionType i, const positionType k, const lengthType Wi, const lengthType Wk, const int& nuc_I, const int& nuc_K, const lengthType& len_I, const lengthType& len_K, stmcell*& cCell, const scoreType& affine ) {

	// Substitution score
	scoreType similarityScore = cCell->getSimilarityScore()
								+ s_matrix.getInit(nuc_I, nuc_K) + affine;


	lengthType len1 = cCell->getLength1() + grow_i;
	lengthType len2 = cCell->getLength2() + grow_Wi;
	lengthType len3 = cCell->getLength3() + grow_k;
	lengthType len4 = cCell->getLength4() + grow_Wk;

	scoreType energyScore = cCell->getEnergyScore()
							- s_matrix.getBulgeLength(len_I, len_K)
							+ s_matrix.getIntLoopLength(len1, len2)
							+ s_matrix.getIntLoopLength(len3, len4);

	// Store (The ones are due to the potential unpaired nuc. in the bp)
	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state, len1, len2, len3, len4, cCell->getPointer());

}


//=============================================================================
// Internal loop functions



template< FOLDK_TEMPLATE >
template<stateType return_state>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIilK2ikWiWk( ENERGY_FUNC ) {

	// If the next set of nucleotides do not base-pair then it is not necessary
	// to do this calculation
	if (!extra_bp_coord_IK(i, k, Wi, Wk)) {return;}

	// Pot il end to stem
	// Add the substitution score
	scoreType similarityScore = cCell->getSimilarityScore()
					+ s_matrix.getScore(seq_1.getPos(i), seq_1.getPos(i+Wi),
	                			        seq_2.getPos(k), seq_2.getPos(k+Wk));


	lengthType len1 = cCell->getLength1();
	lengthType len2 = cCell->getLength2();
	lengthType len3 = cCell->getLength3();
	lengthType len4 = cCell->getLength4();

	// Add the opening and closing stacking score
	scoreType energyScore = cCell->getEnergyScore();
	if ((len1 > 0) && (len2 > 0)) {
		const positionType pos_i =  i + len1;
		const positionType pos_Wi = i + Wi - len2;
		energyScore += s_matrix.getIntLoopOpen(seq_1.getPos(pos_i+1),
		                                 seq_1.getPos(pos_Wi-1),
													seq_1.getPos(pos_i),
													seq_1.getPos(pos_Wi));
		energyScore += s_matrix.getIntLoopClose(
								seq_1.getPos(i), seq_1.getPos(i+Wi), c_i, c_j);
	}
	if ((len3 > 0) && (len4 > 0)) {
		const positionType pos_k =  k + len3;
		const positionType pos_Wk = k + Wk - len4;
		energyScore += s_matrix.getIntLoopOpen(seq_2.getPos(pos_k+1),
		                                 seq_2.getPos(pos_Wk-1),
													seq_2.getPos(pos_k),
													seq_2.getPos(pos_Wk));
		energyScore += s_matrix.getIntLoopClose(
								seq_2.getPos(k), seq_2.getPos(k+Wk), c_k, c_l);
	}

	// Store (The ones are due to the potential unpaired nuc. in the bp)
	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state,
			 len1, len2, len3, len4, cCell->getPointer());

}


template< FOLDK_TEMPLATE >
template<stateType return_state, lengthType grow_i, lengthType grow_k,
         lengthType grow_Wi, lengthType grow_Wk >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::growInterLoop(
                                    const positionType i, const positionType k,
												const lengthType Wi, const lengthType Wk,
												const int& nuc_I, const int& nuc_K,
												stmcell*& cCell, const scoreType& affine) {

	// The potential end is scored as hairpin here, later it will be recalculated
	// to base-pair (if the next nucleotides base-pairs)
	// Pot hp end to stem
	scoreType similarityScore = cCell->getSimilarityScore()
								+ s_matrix.getInit(nuc_I, nuc_K) + affine;

	lengthType len1 = cCell->getLength1() + grow_i;
	lengthType len2 = cCell->getLength2() + grow_Wi;
	lengthType len3 = cCell->getLength3() + grow_k;
	lengthType len4 = cCell->getLength4() + grow_Wk;

	// Subtract the old length cost
	// Add the new length cost
	scoreType energyScore = cCell->getEnergyScore()
						- s_matrix.getIntLoopLength(len1-grow_i, len2-grow_Wi)
						- s_matrix.getIntLoopLength(len3-grow_k, len4-grow_Wk)
						+ s_matrix.getIntLoopLength(len1, len2)
						+ s_matrix.getIntLoopLength(len3, len4);

	// Store (The ones are due to the potential unpaired nuc. in the bp)
	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state,
			 len1, len2, len3, len4, cCell->getPointer());

}
//=============================================================================
// These functions handle mbls


template< FOLDK_TEMPLATE >
template<stateType return_state>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::mIK2ikWiWk( ENERGY_FUNC ) {

	// If the next set of nucleotides do not base-pair then it is not necessary
	// to do this calculation
	if (!extra_bp_coord_IK(i, k, Wi, Wk)) {return;}

	// Add the substitution score
	scoreType similarityScore = cCell->getSimilarityScore()
						+ s_matrix.getScore(seq_1.getPos(i), seq_1.getPos(i+Wi),
	                      			     seq_2.getPos(k), seq_2.getPos(k+Wk));

	// Add the mbl closing cost
	// Add the non-GC stem closing cost
	scoreType energyScore = cCell->getEnergyScore()
					+ s_matrix.getMbl()
					+ s_matrix.getNonGCEnd(seq_1.getPos(i), seq_1.getPos(i+Wi))
					+ s_matrix.getNonGCEnd(seq_2.getPos(k), seq_2.getPos(k+Wk));


	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state,
				0, 0, 0, 0, cCell->getPointer());

}


template< FOLDK_TEMPLATE >
template<stateType return_state, lengthType grow_i, lengthType grow_k,
         lengthType grow_Wi, lengthType grow_Wk>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::mIK2unpaired(
                                    const positionType i, const positionType k,
												const lengthType Wi, const lengthType Wk,
												const int& nuc_I, const int& nuc_K,
												stmcell*& cCell, const scoreType& affine) {


	// Pot hp end to stem
	// The potential end is scored as hairpin here, later it will be recalculated
	// to base-pair (if the next nucleotides base-pairs)
	scoreType similarityScore = cCell->getSimilarityScore()
								+ s_matrix.getInit(nuc_I, nuc_K) + affine;

	lengthType len1 = cCell->getLength1() + grow_i;
	lengthType len2 = cCell->getLength2() + grow_Wi;
	lengthType len3 = cCell->getLength3() + grow_k;
	lengthType len4 = cCell->getLength4() + grow_Wk;


	// Add the mbl unpaired score for both nucleotides
	scoreType energyScore = cCell->getEnergyScore()
							+ (grow_i + grow_k + grow_Wi + grow_Wk)*mblNuc;

	// Store (The ones are due to the potential unpaired nuc. in the bp)
	stmstore(i, k, Wi, Wk, similarityScore, energyScore, return_state, len1, len2, len3, len4,
	         cCell->getPointer());


}


template< FOLDK_TEMPLATE >
template<stateType return_state>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::bWiWkIK2mbl( ENERGY_FUNC ) {

	scoreType similarityScore = cCell->getSimilarityScore();

	// Change the bulge score into the mbl/end score
	scoreType energyScore = cCell->getEnergyScore();
	bulge_WiWk2end(energyScore, i, k, Wi, Wk, cCell);

	mblCore(similarityScore, energyScore, i, k, Wi, Wk, cCell->getLength2(),
			cCell->getLength4(), return_state);
}


template< FOLDK_TEMPLATE >
template<stateType return_state>
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::stemIK2mbl( ENERGY_FUNC ) {

	scoreType energyScore = cCell->getEnergyScore();
	const stateType state = cCell->getState();

	if (array.get_right_branch(state)) {

		// Add the non-GC stem close cost in the stem cases
		// but not in the mbl cases because the score does not apply to the mbl
		// case.

		stem2end(energyScore, i, k, Wi, Wk, cCell);
	}

	mblCore(cCell->getSimilarityScore(), energyScore, i, k, Wi, Wk,
			cCell->getLength2(), cCell->getLength4(), return_state);

}


template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::mblCore(
			    const scoreType left_similarity_score,
				const scoreType left_energy_score,
			    const positionType i, const positionType k,
			    const lengthType Wi, const lengthType Wk,
			    const lengthType len2, const lengthType len4,
			    const stateType return_state) {

	if ( (Wi < min_loop +4) || (Wk < min_loop +4) ) {
		return;
	}

	// Calculate the (i,k) coords of the right part
	const positionType right_i = i + Wi + 1;
	const positionType right_k = k + Wk + 1;

	// Calculate the maximum length of the right part
	const lengthType max_right_Wi = (i + lambda < curr_end_I ) ?
		lambda - Wi -1: curr_end_I -i - Wi;

	const lengthType max_right_Wk = (k + lambda < curr_end_K) ?
		lambda - Wk -1: curr_end_K - k - Wk;

	if (( max_right_Wi < min_loop +4) || (max_right_Wk < min_loop +4)) {
		return;
	}

	const lengthType low_left_delta	= Wi - Wk - delta;
	const lengthType high_left_delta = Wi - Wk + delta;

	if ( (global || mblrealign) &&
			 ((right_k < right_i + k_low) || (right_k > right_i + k_high))) {
		return;
	}

//threadHandler->lockOutput();
//std::cout << "mbl " << thread_number << std::endl;
//threadHandler->unlockOutput();
	two_link< ltmcell >* const subLtm = ltm->getLockAndResetSubMatrix(right_i, right_k, thread_number);
//threadHandler->lockOutput();
//std::cout << "mbl done " << thread_number << std::endl;
//threadHandler->unlockOutput();
//	ltm->resetCurrent(right_i, right_k, thread_number);
	while (true) {
		lengthType right_Wi;
		lengthType right_Wk;

//threadHandler->lockOutput();
//std::cout << "mbl getCell " << thread_number << std::endl;
//threadHandler->unlockOutput();
		// Get the right part information from the ltm
		ltmcell* rightCell = subLtm->getNext(right_Wi, right_Wk);
//threadHandler->lockOutput();
//std::cout << "mbl cell done " << thread_number << std::endl;
//threadHandler->unlockOutput();

		// There are not any potential right parts left at this position
		if (rightCell == 0) {
			break;
			//			ltm->unlock(thread_number);
			//			return;
		}

		// Ensure that the right part has a branch state. This is always true
		// during local scan alignment but during global or realignmnet it is
		// not always true.
		if ( (global || mblrealign) &&
		     ( !array.get_right_branch(rightCell->getState()) )) {
			continue;
		}

		// Make sure that the length of the K-subsequence is no more than delta
		// nucleotides shorter than the I-subsequence
		// Keep getting new right parts until the minimum length is reached.
		if (low_left_delta >= right_Wk - right_Wi) {
			continue;
		}

		// Check the length of the right part
		// The right_Wk > max_right_Wk keeps the total K-subsequence length below
		// lambda
		// high_left_delta etc. Keeps the length difference between the two
		// subsequences below delta.
		if ((right_Wk > max_right_Wk) ||
			(high_left_delta <= right_Wk - right_Wi)) {
			// The K right part is to long. Get the next right Wi
//threadHandler->lockOutput();
//std::cout << "mbl last " << thread_number << std::endl;
//threadHandler->unlockOutput();

			subLtm->lastWk();
//threadHandler->lockOutput();
//std::cout << "mbl last done " << thread_number << std::endl;
//threadHandler->unlockOutput();
			continue;
		}

		// Length check of the length of the right part on the I sequence
		// If it is to long then there are not any possible right parts left
		if (right_Wi > max_right_Wi) {
			//			ltm->deletePos(right_i, right_k);
			//			ltm->unlock(thread_number);
			//			return;
			break;
		}

		const scoreType similarityScore = left_similarity_score
										  + rightCell->getSimilarityScore();

		const scoreType energyScore = left_energy_score
									  + rightCell->getEnergyScore() + mblHelix;

		mbllist* mblpointer = 0;
		if ( (mblrealign) ) {
			// Store the extra multibranch information for later use. This is only
			// done during the mblrealignment. (normal realign should not be
			// branched)

			mbllist mbl(i, k, Wi - len2, Wk - len4,
			right_i, right_k, right_Wi, right_Wk, similarityScore, energyScore);
			mblpointer = local_mblm->push_ptr(mbl);
		}

		stmstore(i, k, Wi + right_Wi+1, Wk + right_Wk +1, similarityScore,
					energyScore, return_state, 0,0,0,0, mblpointer);

	}
	ltm->unlock(right_i, right_k, thread_number);
}

//=================================================================
// FLOW CONTROL ARRAYS BELOW


template< FOLDK_TEMPLATE >
inline void foldK< FOLDK_TEMPLATE_PARAMETERS >::buildFlowArrays() {

// Diane, never drink coffee that has been anywhere near a fish.

// These functions is used to pick the function dependent on the previous state
// they controls how a new score is calculated dependend on the previous state

	// Init all to give a big negative score
	helper::init_array(p2calc_ikWiWk, flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::errorFunction);
	helper::init_array(p2calc_iWi,   flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc);
	helper::init_array(p2calc_kWk,   flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc);
	helper::init_array(p2calc_ik,   flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::errorFunction);
	helper::init_array(p2calc_WiWk,   flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::errorFunction);
	helper::init_array(p2calc_i,    flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::errorFunction);
	helper::init_array(p2calc_Wi,    flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::errorFunction);
	helper::init_array(p2calc_k,    flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::errorFunction);
	helper::init_array(p2calc_Wk,    flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::errorFunction);
	helper::init_array(p2calc_mbl,    flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc);
	helper::init_array(p2end,         flow_size, &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc);

	p2calc_ikWiWk[hp_init]      						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_align_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_align_WiWk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_I_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_I_WiWk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_K_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_K_WiWk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_Wi_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_Wi_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_Wk_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_Wk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_iWi_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_iWi_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_kWk_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_kWk_Wi]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_iWk_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_iWk_Wi]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_kWi_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_init_gap_kWi_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_init_pot_bp_IK<hp_pb_IK>;
	p2calc_ikWiWk[hp_pb_IK]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ikWiWk<stem_IK>;				// -> stem_IK
	p2calc_ikWiWk[stem_IK]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ikWiWk<stem_IK>;				//	stem_IK
	p2calc_ikWiWk[stem_I_stem_gap_kWk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ikWiWk<stem_IK>;	// stem_IK
	p2calc_ikWiWk[stem_no_mbl_I_stem_gap_kWk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ikWiWk<stem_IK>;	// stem_IK
	p2calc_ikWiWk[stem_gap_iWi_stem_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ikWiWk<stem_IK>;	// stem_IK
	p2calc_ikWiWk[stem_no_mbl_gap_iWi_stem_K]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ikWiWk<stem_IK>;	// stem_IK
	p2calc_ikWiWk[bi_I_bk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2ikWiWk<bi_pb_I_bk_pb_K>;	// bi_pb_I_bk_pb_K
	p2calc_ikWiWk[bi_I_bk_gap_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2ikWiWk<bi_pb_I_bk_pb_K>;	// bi_pb_I_bk_pb_K
	p2calc_ikWiWk[bi_I_bk_gap_kWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2ikWiWk<bi_pb_I_bk_pb_K>;	// bi_pb_I_bk_pb_K
	p2calc_ikWiWk[bi_I_bk_gap_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2ikWiWk<bi_pb_I_bk_pb_K>;	// bi_pb_I_bk_pb_K
	p2calc_ikWiWk[bi_pb_I_bk_pb_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ikWiWk<stem_IK>;		// stem_IK
	p2calc_ikWiWk[bi_gap_i_bk_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2ikWiWk<bi_pb_I_bk_pb_K>;	// bi_pb_I_bk_pb_K
	p2calc_ikWiWk[bi_gap_i_bk_gap_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2ikWiWk<bi_pb_I_bk_pb_K>;	// bi_pb_I_bk_pb_K
	p2calc_ikWiWk[bi_gap_iWi_bk_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2ikWiWk<bi_pb_I_bk_pb_K>;	// bi_pb_I_bk_pb_K
	p2calc_ikWiWk[bi_gap_Wi_bk_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2ikWiWk<bi_pb_I_bk_pb_K>;	// bi_pb_I_bk_pb_K
	p2calc_ikWiWk[bi_gap_Wi_bk_gap_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2ikWiWk<bi_pb_I_bk_pb_K>;	// bi_pb_I_bk_pb_K
	p2calc_ikWiWk[bWi_I_bWk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ikWiWk<bWi_pb_I_bWk_pb_K>;	// bWi_pb_I_bWk_pb_K
	p2calc_ikWiWk[bWi_I_bWk_gap_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ikWiWk<bWi_pb_I_bWk_pb_K>;	// bWi_pb_I_bWk_pb_K
	p2calc_ikWiWk[bWi_I_bWk_gap_kWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ikWiWk<bWi_pb_I_bWk_pb_K>;	// bWi_pb_I_bWk_pb_K
	p2calc_ikWiWk[bWi_I_bWk_gap_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ikWiWk<bWi_pb_I_bWk_pb_K>;	// bWi_pb_I_bWk_pb_K
	p2calc_ikWiWk[bWi_pb_I_bWk_pb_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ikWiWk<stem_IK>;		// stem_IK
	p2calc_ikWiWk[bWi_gap_i_bWk_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ikWiWk<bWi_pb_I_bWk_pb_K>;	// bWi_pb_I_bWk_pb_K
	p2calc_ikWiWk[bWi_gap_i_bWk_gap_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ikWiWk<bWi_pb_I_bWk_pb_K>;	// bWi_pb_I_bWk_pb_K
	p2calc_ikWiWk[bWi_gap_iWi_bWk_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ikWiWk<bWi_pb_I_bWk_pb_K>;	// bWi_pb_I_bWk_pb_K
	p2calc_ikWiWk[bWi_gap_Wi_bWk_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ikWiWk<bWi_pb_I_bWk_pb_K>;	// bWi_pb_I_bWk_pb_K
	p2calc_ikWiWk[bWi_gap_Wi_bWk_gap_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ikWiWk<bWi_pb_I_bWk_pb_K>;	// bWi_pb_I_bWk_pb_K
	p2calc_ikWiWk[il_I_il_K_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_I_il_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_I_il_gap_k_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_I_il_gap_k_WiWk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_I_il_gap_Wk_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_I_il_gap_Wk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_I_il_gap_kWk_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_I_il_gap_kWk_Wi]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_pb_I_il_pb_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ikWiWk<stem_IK>;		// stemIK
	p2calc_ikWiWk[il_gap_iWi_il_K_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_gap_iWi_il_K_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_gap_i_il_K_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_gap_i_il_K_WiWk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_gap_i_il_gap_Wk_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_gap_i_il_gap_Wk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_gap_Wi_il_K_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_gap_Wi_il_K_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_gap_Wi_il_gap_k_i]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[il_gap_Wi_il_gap_k_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ikWiWk<il_pb_I_il_pb_K>;	// il_pb_I_il_pb_K
	p2calc_ikWiWk[mblIK]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_bWi_I_mbl_bWk_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_bWi_I_mbl_bWk_gap_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_bWi_gap_Wi_mbl_bWk_K]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_I_mbl_il_K_ik]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_I_mbl_il_K_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_I_mbl_il_gap_k_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_I_mbl_il_gap_k_WiWk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_I_mbl_il_gap_kWk_i]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_I_mbl_il_gap_kWk_Wi]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_I_mbl_il_gap_Wk_ik]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_I_mbl_il_gap_Wk_Wi]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_pb_I_mbl_il_pb_K]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ikWiWk<stem_IK>;		// stem_IK
	p2calc_ikWiWk[mbl_il_gap_Wi_mbl_il_K_ik]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_gap_Wi_mbl_il_K_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_gap_Wi_mbl_il_gap_k_i]	= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_gap_Wi_mbl_il_gap_k_Wk]	= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_gap_iWi_mbl_il_K_k]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_gap_iWi_mbl_il_K_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_gap_i_mbl_il_K_k]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_gap_i_mbl_il_K_WiWk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_gap_i_mbl_il_gap_Wk_k]	= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K
	p2calc_ikWiWk[mbl_il_gap_i_mbl_il_gap_Wk_Wi]	= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ikWiWk<mbl_il_pb_I_mbl_il_pb_K>;	// mbl_il_pb_I_mbl_il_pb_K

	p2calc_iWi[hp_pb_IK] 								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2iWi<stem_no_mbl_I_stem_gap_kWk>;					// hp_pb_I_hp_gap_kWk
	p2calc_iWi[stem_IK] 									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2iWi<stem_I_stem_gap_kWk>;					// stem_I_stem_gap_kWk
	p2calc_iWi[stem_I_stem_gap_kWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIstemgkWk2iWi<stem_I_stem_gap_kWk>;				// stem_I_stem_gap_kWk
	p2calc_iWi[stem_no_mbl_I_stem_gap_kWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIstemgkWk2iWi<stem_no_mbl_I_stem_gap_kWk>;				// stem_I_stem_gap_kWk
	p2calc_iWi[stem_gap_iWi_stem_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2iWi<stem_I_stem_gap_kWk>;
	p2calc_iWi[stem_no_mbl_gap_iWi_stem_K]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2iWi<stem_I_stem_gap_kWk>;
	p2calc_iWi[bi_pb_I_bk_pb_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2iWi<stem_no_mbl_I_stem_gap_kWk>;			// stem_I_bk_pb_gap_kWk
	p2calc_iWi[bWi_pb_I_bWk_pb_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2iWi<stem_no_mbl_I_stem_gap_kWk>;			// stem_I_bWk_pb_gap_kWk
	p2calc_iWi[il_pb_I_il_pb_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2iWi<stem_no_mbl_I_stem_gap_kWk>;			// stem_I_il_pb_gap_kWk
	p2calc_iWi[mbl_il_pb_I_mbl_il_pb_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2iWi<stem_no_mbl_I_stem_gap_kWk>;			// stem_I_mbl_il_pb_gap_kWk


	p2calc_kWk[hp_pb_IK]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2kWk<stem_no_mbl_gap_iWi_stem_K>;
	p2calc_kWk[stem_IK]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2kWk<stem_gap_iWi_stem_K>;					// stem_gap_iWi_stem_K
	p2calc_kWk[stem_I_stem_gap_kWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2kWk<stem_gap_iWi_stem_K>;		// stem_gap_iWi_stem_K
	p2calc_kWk[stem_no_mbl_I_stem_gap_kWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2kWk<stem_gap_iWi_stem_K>;		// stem_gap_iWi_stem_K
	p2calc_kWk[stem_gap_iWi_stem_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemgiWistemK2kWk<stem_gap_iWi_stem_K>;		// stem_gap_iWi_stem_K
	p2calc_kWk[stem_no_mbl_gap_iWi_stem_K]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemgiWistemK2kWk<stem_no_mbl_gap_iWi_stem_K>;		// stem_gap_iWi_stem_K
	p2calc_kWk[bi_pb_I_bk_pb_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2kWk<stem_no_mbl_gap_iWi_stem_K>;			// bi_pb_gap_iWi_stem_K
	p2calc_kWk[bWi_pb_I_bWk_pb_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2kWk<stem_no_mbl_gap_iWi_stem_K>;			// bWi_pb_gap_iWi_stem_K
	p2calc_kWk[il_pb_I_il_pb_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2kWk<stem_no_mbl_gap_iWi_stem_K>;			// il_pb_gap_iWi_stem_K
	p2calc_kWk[mbl_il_pb_I_mbl_il_pb_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2kWk<stem_no_mbl_gap_iWi_stem_K>;			// mbl_il_pb_gap_iWi_stem_K


	p2calc_ik[hp_init]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_align_ik>;		 // hp_init_align
	p2calc_ik[hp_init_align_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_align_ik>;		 // hp_init_align
	p2calc_ik[hp_init_align_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_align_ik>;		 // hp_init_align
	p2calc_ik[hp_init_gap_I_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_align_ik>;		 // hp_init_align
	p2calc_ik[hp_init_gap_I_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_align_ik>;		 // hp_init_align
	p2calc_ik[hp_init_gap_K_i]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_align_ik>;		 // hp_init_align
	p2calc_ik[hp_init_gap_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_align_ik>;		 // hp_init_align
	p2calc_ik[hp_init_gap_Wi_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wi_ik>;   // hp_init_gap_Wi
	p2calc_ik[hp_init_gap_Wi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wi_ik>;   // hp_init_gap_Wi
	p2calc_ik[hp_init_gap_Wk_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wk_ik>;   // hp_init_gap_Wk
	p2calc_ik[hp_init_gap_Wk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wk_ik>;   // hp_init_gap_Wk
	p2calc_ik[hp_init_gap_iWi_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wi_ik>;  // hp_init_gap_Wi
	p2calc_ik[hp_init_gap_iWi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wi_ik>;  // hp_init_gap_Wi
	p2calc_ik[hp_init_gap_kWk_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wk_ik>;  // hp_init_gap_Wk
	p2calc_ik[hp_init_gap_kWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wk_ik>;  // hp_init_gap_Wk
	p2calc_ik[hp_init_gap_iWk_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wk_ik>;  // hp_init_gap_Wk
	p2calc_ik[hp_init_gap_iWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wk_ik>;  // hp_init_gap_Wk
	p2calc_ik[hp_init_gap_kWi_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wi_ik>;  // hp_init_gap_Wi
	p2calc_ik[hp_init_gap_kWi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_ik<hp_init_gap_Wi_ik>;  // hp_init_gap_Wi
	p2calc_ik[hp_pb_IK]    								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;		 // hp_init_align
	p2calc_ik[stem_IK]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_I_bk_K>;					// bi_I_bk_K
	p2calc_ik[stem_I_stem_gap_kWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_I_bk_gap_Wk>;			// bi_I_bk_gap_Wk
	p2calc_ik[stem_no_mbl_I_stem_gap_kWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_ik[stem_gap_iWi_stem_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_gap_Wi_bk_K>;			// bi_gap_Wi_bk_K
	p2calc_ik[stem_no_mbl_gap_iWi_stem_K]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_ik[bi_I_bk_K]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_I_bk_K>;					// bi_I_bk_K
	p2calc_ik[bi_I_bk_gap_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_I_bk_K>;					// bi_I_bk_K
	p2calc_ik[bi_I_bk_gap_kWk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_I_bk_gap_Wk>;				 // bi_I_bk_K
	p2calc_ik[bi_I_bk_gap_Wk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_I_bk_gap_Wk>;				 // bi_I_bk_K
	p2calc_ik[bi_pb_I_bk_pb_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;  			 // il_I_il_K
	p2calc_ik[bi_gap_i_bk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_I_bk_K>;  				 // bi_I_bk_K
	p2calc_ik[bi_gap_i_bk_gap_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_I_bk_gap_Wk>;				 // bi_I_bk_K
	p2calc_ik[bi_gap_iWi_bk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_gap_Wi_bk_K>;				 // bi_I_bk_K
	p2calc_ik[bi_gap_Wi_bk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_gap_Wi_bk_K>;				 // bi_I_bk_K
	p2calc_ik[bi_gap_Wi_bk_gap_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2ik<bi_gap_Wi_bk_K>;				 // bi_I_bk_K
	p2calc_ik[bWi_I_bWk_K]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ik<il_I_il_K_ik>;				// il_I_il_K
	p2calc_ik[bWi_I_bWk_gap_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ik<il_I_il_K_ik>;				// il_I_il_K
	p2calc_ik[bWi_I_bWk_gap_kWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ik<il_I_il_gap_Wk_ik>;				// il_I_il_K
	p2calc_ik[bWi_I_bWk_gap_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ik<il_I_il_gap_Wk_ik>;				// il_I_il_K
	p2calc_ik[bWi_pb_I_bWk_pb_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;				// il_I_il_K
	p2calc_ik[bWi_gap_i_bWk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ik<il_I_il_K_ik>;				// il_I_il_K
	p2calc_ik[bWi_gap_i_bWk_gap_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ik<il_I_il_gap_Wk_ik>;				// il_I_il_K
	p2calc_ik[bWi_gap_iWi_bWk_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ik<il_gap_Wi_il_K_ik>;				// il_I_il_K
	p2calc_ik[bWi_gap_Wi_bWk_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ik<il_gap_Wi_il_K_ik>;				// il_I_il_K
	p2calc_ik[bWi_gap_Wi_bWk_gap_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2ik<il_gap_Wi_il_K_ik>;				// il_I_il_K
	p2calc_ik[il_I_il_K_ik]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_K_ik>;					// il_I_il_K
	p2calc_ik[il_I_il_K_WiWk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_K_ik>;					// il_I_il_K
	p2calc_ik[il_I_il_gap_k_i]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_K_ik>;					// il_I_il
	p2calc_ik[il_I_il_gap_k_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_K_ik>;					// il_I_il
	p2calc_ik[il_I_il_gap_Wk_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_gap_Wk_ik>;				// il_I_il_gap_Wk
	p2calc_ik[il_I_il_gap_Wk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_gap_Wk_ik>;				// il_I_il_gap_Wk
	p2calc_ik[il_I_il_gap_kWk_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_gap_Wk_ik>;				// il_I_il_gap_Wk
	p2calc_ik[il_I_il_gap_kWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_gap_Wk_ik>;				// il_I_il_gap_Wk
	p2calc_ik[il_pb_I_il_pb_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;				// il_I_il_K
	p2calc_ik[il_gap_iWi_il_K_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_gap_Wi_il_K_ik>;				// il_gap_Wi_il_K
	p2calc_ik[il_gap_iWi_il_K_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_gap_Wi_il_K_ik>;				// il_gap_Wi_il_K
	p2calc_ik[il_gap_i_il_K_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_K_ik>;					// il_I_il_K
	p2calc_ik[il_gap_i_il_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_K_ik>;					// il_I_il_K
	p2calc_ik[il_gap_i_il_gap_Wk_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_gap_Wk_ik>;				// il_I_il_K
	p2calc_ik[il_gap_i_il_gap_Wk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_I_il_gap_Wk_ik>;				// il_I_il_K
	p2calc_ik[il_gap_Wi_il_K_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_gap_Wi_il_K_ik>;				// il_gap_Wi_il_K
	p2calc_ik[il_gap_Wi_il_K_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_gap_Wi_il_K_ik>;				// il_gap_Wi_il_K
	p2calc_ik[il_gap_Wi_il_gap_k_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_gap_Wi_il_K_ik>;				// il_gap_Wi_il_K
	p2calc_ik[il_gap_Wi_il_gap_k_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2ik<il_gap_Wi_il_K_ik>;				// il_gap_Wi_il_K
	p2calc_ik[mblIK]  									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_K_ik>;						// mbl_bi_I_mbl_bk_K
	p2calc_ik[mbl_bWi_I_mbl_bWk_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_K_ik>;						// mbl_bi_I_mbl_bk_K
	p2calc_ik[mbl_bWi_I_mbl_bWk_gap_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_gap_Wk_ik>;						// mbl_bi_I_mbl_bk_K
	p2calc_ik[mbl_bWi_gap_Wi_mbl_bWk_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_gap_Wi_mbl_il_K_ik>;						// mbl_bi_I_mbl_bk_K
	p2calc_ik[mbl_il_I_mbl_il_K_ik]  				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_K_ik>;				// mbl_il_I_mbl_il_K
	p2calc_ik[mbl_il_I_mbl_il_K_WiWk]  				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_K_ik>;				// mbl_il_I_mbl_il_K
	p2calc_ik[mbl_il_I_mbl_il_gap_k_i]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_K_ik>;				// mbl_il_I_mbl_il_K
	p2calc_ik[mbl_il_I_mbl_il_gap_k_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_K_ik>;				// mbl_il_I_mbl_il_K
	p2calc_ik[mbl_il_I_mbl_il_gap_kWk_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_gap_Wk_ik>;			// mbl_il_I_mbl_il_gap_Wk
	p2calc_ik[mbl_il_I_mbl_il_gap_kWk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_gap_Wk_ik>;			// mbl_il_I_mbl_il_gap_Wk
	p2calc_ik[mbl_il_I_mbl_il_gap_Wk_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_gap_Wk_ik>;
	p2calc_ik[mbl_il_I_mbl_il_gap_Wk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_gap_Wk_ik>;
	p2calc_ik[mbl_il_pb_I_mbl_il_pb_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_ik[mbl_il_gap_Wi_mbl_il_K_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_gap_Wi_mbl_il_K_ik>;
	p2calc_ik[mbl_il_gap_Wi_mbl_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_gap_Wi_mbl_il_K_ik>;
	p2calc_ik[mbl_il_gap_Wi_mbl_il_gap_k_i]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_gap_Wi_mbl_il_K_ik>;
	p2calc_ik[mbl_il_gap_Wi_mbl_il_gap_k_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_gap_Wi_mbl_il_K_ik>;
	p2calc_ik[mbl_il_gap_iWi_mbl_il_K_k]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_gap_Wi_mbl_il_K_ik>;
	p2calc_ik[mbl_il_gap_iWi_mbl_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_gap_Wi_mbl_il_K_ik>;
	p2calc_ik[mbl_il_gap_i_mbl_il_K_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_K_ik>;
	p2calc_ik[mbl_il_gap_i_mbl_il_K_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_K_ik>;
	p2calc_ik[mbl_il_gap_i_mbl_il_gap_Wk_k]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_gap_Wk_ik>;
	p2calc_ik[mbl_il_gap_i_mbl_il_gap_Wk_Wi]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2ik<mbl_il_I_mbl_il_gap_Wk_ik>;


	p2calc_WiWk[hp_init] 								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_align_WiWk>;  	// hp_init_align
	p2calc_WiWk[hp_init_align_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_align_WiWk>;
	p2calc_WiWk[hp_init_align_WiWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_align_WiWk>;
	p2calc_WiWk[hp_init_gap_I_k] 						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_I_WiWk>;    // hp_init_gap_I
	p2calc_WiWk[hp_init_gap_I_WiWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_I_WiWk>;    // hp_init_gap_I
	p2calc_WiWk[hp_init_gap_K_i] 						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_K_WiWk>;    // hp_init_gap_K
	p2calc_WiWk[hp_init_gap_K_WiWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_K_WiWk>;    // hp_init_gap_K
	p2calc_WiWk[hp_init_gap_Wi_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_align_WiWk>;
	p2calc_WiWk[hp_init_gap_Wi_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_align_WiWk>;
	p2calc_WiWk[hp_init_gap_Wk_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_align_WiWk>;
	p2calc_WiWk[hp_init_gap_Wk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_align_WiWk>;
	p2calc_WiWk[hp_init_gap_iWi_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_I_WiWk>;
	p2calc_WiWk[hp_init_gap_iWi_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_I_WiWk>;
	p2calc_WiWk[hp_init_gap_kWk_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_K_WiWk>;
	p2calc_WiWk[hp_init_gap_kWk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_K_WiWk>;
	p2calc_WiWk[hp_init_gap_iWk_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_I_WiWk>;
	p2calc_WiWk[hp_init_gap_iWk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_I_WiWk>;
	p2calc_WiWk[hp_init_gap_kWi_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_K_WiWk>;
	p2calc_WiWk[hp_init_gap_kWi_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_WiWk<hp_init_gap_K_WiWk>;
	p2calc_WiWk[hp_pb_IK]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_WiWk[stem_IK] 								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_I_bWk_K>;				// bWi_I_bWk_K
	p2calc_WiWk[stem_I_stem_gap_kWk] 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_I_bWk_gap_k>;				// bWi_I_bWk_K
	p2calc_WiWk[stem_no_mbl_I_stem_gap_kWk]	 	= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_WiWk[stem_gap_iWi_stem_K] 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_gap_i_bWk_K>;				// bWi_I_bWk_K
	p2calc_WiWk[stem_no_mbl_gap_iWi_stem_K] 		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_WiWk[bi_I_bk_K]				 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2WiWk<il_I_il_K_WiWk>;
	p2calc_WiWk[bi_I_bk_gap_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2WiWk<il_I_il_gap_k_WiWk>;
	p2calc_WiWk[bi_I_bk_gap_kWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2WiWk<il_I_il_gap_k_WiWk>;
	p2calc_WiWk[bi_I_bk_gap_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2WiWk<il_I_il_K_WiWk>;
	p2calc_WiWk[bi_pb_I_bk_pb_K]  					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_WiWk[bi_gap_i_bk_K] 		 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2WiWk<il_gap_i_il_K_WiWk>;
	p2calc_WiWk[bi_gap_i_bk_gap_Wk] 		 			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2WiWk<il_gap_i_il_K_WiWk>;
	p2calc_WiWk[bi_gap_iWi_bk_K] 		 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2WiWk<il_gap_i_il_K_WiWk>;
	p2calc_WiWk[bi_gap_Wi_bk_K] 		 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2WiWk<il_I_il_K_WiWk>;
	p2calc_WiWk[bi_gap_Wi_bk_gap_k]	 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2WiWk<il_I_il_gap_k_WiWk>;
	p2calc_WiWk[bWi_I_bWk_K]  							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_I_bWk_K>;
	p2calc_WiWk[bWi_I_bWk_gap_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_I_bWk_gap_k>;
	p2calc_WiWk[bWi_I_bWk_gap_kWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_I_bWk_gap_k>;
	p2calc_WiWk[bWi_I_bWk_gap_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_I_bWk_K>;
	p2calc_WiWk[bWi_pb_I_bWk_pb_K]  					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_WiWk[bWi_gap_iWi_bWk_K] 	 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_gap_i_bWk_K>;
	p2calc_WiWk[bWi_gap_i_bWk_K] 	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_gap_i_bWk_K>;
	p2calc_WiWk[bWi_gap_i_bWk_gap_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_gap_i_bWk_K>;
	p2calc_WiWk[bWi_gap_Wi_bWk_K] 	 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_I_bWk_K>;
	p2calc_WiWk[bWi_gap_Wi_bWk_gap_k] 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2WiWk<bWi_I_bWk_gap_k>;
	p2calc_WiWk[il_I_il_K_ik]  						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_K_WiWk>;
	p2calc_WiWk[il_I_il_K_WiWk]  						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_K_WiWk>;
	p2calc_WiWk[il_I_il_gap_k_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_gap_k_WiWk>;
	p2calc_WiWk[il_I_il_gap_k_WiWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_gap_k_WiWk>;
	p2calc_WiWk[il_I_il_gap_Wk_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_K_WiWk>;
	p2calc_WiWk[il_I_il_gap_Wk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_K_WiWk>;
	p2calc_WiWk[il_I_il_gap_kWk_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_gap_k_WiWk>;
	p2calc_WiWk[il_I_il_gap_kWk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_gap_k_WiWk>;
	p2calc_WiWk[il_pb_I_il_pb_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_WiWk[il_gap_iWi_il_K_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_gap_i_il_K_WiWk>;
	p2calc_WiWk[il_gap_iWi_il_K_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_gap_i_il_K_WiWk>;
	p2calc_WiWk[il_gap_i_il_K_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_gap_i_il_K_WiWk>;
	p2calc_WiWk[il_gap_i_il_K_WiWk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_gap_i_il_K_WiWk>;
	p2calc_WiWk[il_gap_i_il_gap_Wk_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_gap_i_il_K_WiWk>;
	p2calc_WiWk[il_gap_i_il_gap_Wk_Wi]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_gap_i_il_K_WiWk>;
	p2calc_WiWk[il_gap_Wi_il_K_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_K_WiWk>;
	p2calc_WiWk[il_gap_Wi_il_K_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_K_WiWk>;
	p2calc_WiWk[il_gap_Wi_il_gap_k_i]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_gap_k_WiWk>;
	p2calc_WiWk[il_gap_Wi_il_gap_k_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2WiWk<il_I_il_gap_k_WiWk>;
	p2calc_WiWk[mblIK]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_bWi_I_mbl_bWk_K>;				 // mbl_bWi_I_mbl_bWk_K
	p2calc_WiWk[mbl_bWi_I_mbl_bWk_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_bWi_I_mbl_bWk_K>;
	p2calc_WiWk[mbl_bWi_I_mbl_bWk_gap_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_bWi_I_mbl_bWk_K>;
	p2calc_WiWk[mbl_bWi_gap_Wi_mbl_bWk_K]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_bWi_I_mbl_bWk_K>;
	p2calc_WiWk[mbl_il_I_mbl_il_K_ik]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_I_mbl_il_K_WiWk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_I_mbl_il_gap_k_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_gap_k_WiWk>;
	p2calc_WiWk[mbl_il_I_mbl_il_gap_k_WiWk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_gap_k_WiWk>;
	p2calc_WiWk[mbl_il_I_mbl_il_gap_kWk_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_gap_k_WiWk>;
	p2calc_WiWk[mbl_il_I_mbl_il_gap_kWk_Wi]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_gap_k_WiWk>;
	p2calc_WiWk[mbl_il_I_mbl_il_gap_Wk_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_I_mbl_il_gap_Wk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_pb_I_mbl_il_pb_K]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_WiWk[mbl_il_gap_Wi_mbl_il_K_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_gap_Wi_mbl_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_gap_Wi_mbl_il_gap_k_i]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_gap_k_WiWk>;
	p2calc_WiWk[mbl_il_gap_Wi_mbl_il_gap_k_Wk]	= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_I_mbl_il_gap_k_WiWk>;
	p2calc_WiWk[mbl_il_gap_iWi_mbl_il_K_k]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_gap_i_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_gap_iWi_mbl_il_K_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_gap_i_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_gap_i_mbl_il_K_k]  			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_gap_i_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_gap_i_mbl_il_K_WiWk]  		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_gap_i_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_gap_i_mbl_il_gap_Wk_k]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_gap_i_mbl_il_K_WiWk>;
	p2calc_WiWk[mbl_il_gap_i_mbl_il_gap_Wk_Wi]	= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2WiWk<mbl_il_gap_i_mbl_il_K_WiWk>;


	p2calc_i[hp_init]										= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_K_i>;		 // hp_init_align
	p2calc_i[hp_init_align_ik]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_K_i>;		 // hp_init_align
	p2calc_i[hp_init_align_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_K_i>;		 // hp_init_align
	p2calc_i[hp_init_gap_I_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_K_i>;		 // hp_init_align
	p2calc_i[hp_init_gap_I_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_K_i>;		 // hp_init_align
	p2calc_i[hp_init_gap_K_i]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_k_i<hp_init_gap_K_i>;		 // hp_init_align
	p2calc_i[hp_init_gap_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_k_i<hp_init_gap_K_i>;		 // hp_init_align
	p2calc_i[hp_init_gap_Wi_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_kWi_i>;   // hp_init_gap_Wi
	p2calc_i[hp_init_gap_Wi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_kWi_i>;   // hp_init_gap_Wi
	p2calc_i[hp_init_gap_Wk_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_kWk_i>;   // hp_init_gap_Wk
	p2calc_i[hp_init_gap_Wk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_kWk_i>;   // hp_init_gap_Wk
	p2calc_i[hp_init_gap_iWi_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_kWi_i>;  // hp_init_gap_Wi
	p2calc_i[hp_init_gap_iWi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_kWi_i>;  // hp_init_gap_Wi
	p2calc_i[hp_init_gap_kWk_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_k_i<hp_init_gap_kWk_i>;  // hp_init_gap_Wk
	p2calc_i[hp_init_gap_kWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_k_i<hp_init_gap_kWk_i>;  // hp_init_gap_Wk
	p2calc_i[hp_init_gap_iWk_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_kWk_i>;  // hp_init_gap_Wk
	p2calc_i[hp_init_gap_iWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_i<hp_init_gap_kWk_i>;  // hp_init_gap_Wk
	p2calc_i[hp_init_gap_kWi_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_k_i<hp_init_gap_kWi_i>;  // hp_init_gap_Wi
	p2calc_i[hp_init_gap_kWi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_k_i<hp_init_gap_kWi_i>;  // hp_init_gap_Wi
	p2calc_i[hp_pb_IK]    								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;  	  // hp_init_align
	p2calc_i[stem_IK]										= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2i<bi_I_bk_gap_k>;					// bi_I_bk_K
	p2calc_i[stem_I_stem_gap_kWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2i<bi_I_bk_gap_kWk>;
	p2calc_i[stem_no_mbl_I_stem_gap_kWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_i[stem_gap_iWi_stem_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2i<bi_gap_Wi_bk_gap_k>;			// bi_gap_Wi_bk_K
	p2calc_i[stem_no_mbl_gap_iWi_stem_K]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_i[bi_I_bk_K]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2i<bi_I_bk_gap_k>;
	p2calc_i[bi_I_bk_gap_k]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2i<bi_I_bk_gap_k>;  			  // bi_I_bk_K
	p2calc_i[bi_I_bk_gap_kWk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2i<bi_I_bk_gap_kWk>;			  // bi_I_bk_K
	p2calc_i[bi_I_bk_gap_Wk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2i<bi_I_bk_gap_kWk>;			  // bi_I_bk_K
	p2calc_i[bi_pb_I_bk_pb_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;			 // il_I_il_K
	p2calc_i[bi_gap_i_bk_K]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2i<bi_I_bk_gap_k>;  				 // bi_I_bk_K
	p2calc_i[bi_gap_i_bk_gap_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2i<bi_I_bk_gap_kWk>;				 // bi_I_bk_K
	p2calc_i[bi_gap_iWi_bk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2i<bi_gap_Wi_bk_gap_k>;				 // bi_I_bk_K
	p2calc_i[bi_gap_Wi_bk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2i<bi_gap_Wi_bk_gap_k>;				 // bi_I_bk_K
	p2calc_i[bi_gap_Wi_bk_gap_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2i<bi_gap_Wi_bk_gap_k>;				 // bi_I_bk_K
	p2calc_i[bWi_I_bWk_K]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2i<il_I_il_gap_k_i>;			 // il_I_il_K
	p2calc_i[bWi_I_bWk_gap_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkKgap2i<il_I_il_gap_k_i>;			 // il_I_il_K
	p2calc_i[bWi_I_bWk_gap_kWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkKgap2i<il_I_il_gap_kWk_i>; 			 // il_I_il_K
	p2calc_i[bWi_I_bWk_gap_Wk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2i<il_I_il_gap_kWk_i>; 			 // il_I_il_K
	p2calc_i[bWi_pb_I_bWk_pb_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;			 // il_I_il_K
	p2calc_i[bWi_gap_i_bWk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2i<il_I_il_gap_k_i>;			 // il_I_il_K
	p2calc_i[bWi_gap_i_bWk_gap_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2i<il_I_il_gap_kWk_i>; 			 // il_I_il_K
	p2calc_i[bWi_gap_iWi_bWk_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2i<il_gap_Wi_il_gap_k_i>; 			 // il_I_il_K
	p2calc_i[bWi_gap_Wi_bWk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2i<il_gap_Wi_il_gap_k_i>; 			 // il_I_il_K
	p2calc_i[bWi_gap_Wi_bWk_gap_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkKgap2i<il_gap_Wi_il_gap_k_i>; 			 // il_I_il_K
	p2calc_i[il_I_il_K_ik]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_I_il_gap_k_i>;					// il_I_il_K
	p2calc_i[il_I_il_K_WiWk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_I_il_gap_k_i>;					// il_I_il_K
	p2calc_i[il_I_il_gap_k_i]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2i<il_I_il_gap_k_i>;					// il_I_il
	p2calc_i[il_I_il_gap_k_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2i<il_I_il_gap_k_i>;					// il_I_il
	p2calc_i[il_I_il_gap_Wk_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_I_il_gap_kWk_i>;				// il_I_il_gap_Wk
	p2calc_i[il_I_il_gap_Wk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_I_il_gap_kWk_i>;				// il_I_il_gap_Wk
	p2calc_i[il_I_il_gap_kWk_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2i<il_I_il_gap_kWk_i>;				// il_I_il_gap_Wk
	p2calc_i[il_I_il_gap_kWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2i<il_I_il_gap_kWk_i>;				// il_I_il_gap_Wk
	p2calc_i[il_pb_I_il_pb_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;				// il_I_il_K
	p2calc_i[il_gap_iWi_il_K_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_gap_Wi_il_gap_k_i>;				// il_gap_Wi_il_K
	p2calc_i[il_gap_iWi_il_K_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_gap_Wi_il_gap_k_i>;				// il_gap_Wi_il_K
	p2calc_i[il_gap_i_il_K_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_I_il_gap_k_i>;					// il_I_il_K
	p2calc_i[il_gap_i_il_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_I_il_gap_k_i>;					// il_I_il_K
	p2calc_i[il_gap_i_il_gap_Wk_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_I_il_gap_kWk_i>;				// il_I_il_K
	p2calc_i[il_gap_i_il_gap_Wk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_I_il_gap_kWk_i>;				// il_I_il_K
	p2calc_i[il_gap_Wi_il_K_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_gap_Wi_il_gap_k_i>;				// il_gap_Wi_il_K
	p2calc_i[il_gap_Wi_il_K_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2i<il_gap_Wi_il_gap_k_i>;				// il_gap_Wi_il_K
	p2calc_i[il_gap_Wi_il_gap_k_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2i<il_gap_Wi_il_gap_k_i>;				// il_gap_Wi_il_K
	p2calc_i[il_gap_Wi_il_gap_k_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2i<il_gap_Wi_il_gap_k_i>;				// il_gap_Wi_il_K
	p2calc_i[mblIK]  										= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_k_i>;						// mbl_bi_I_mbl_bk_K
	p2calc_i[mbl_bWi_I_mbl_bWk_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_k_i>;						// mbl_bi_I_mbl_bk_K
	p2calc_i[mbl_bWi_I_mbl_bWk_gap_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_kWk_i>;						// mbl_bi_I_mbl_bk_K
	p2calc_i[mbl_bWi_gap_Wi_mbl_bWk_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_gap_Wi_mbl_il_gap_k_i>;						// mbl_bi_I_mbl_bk_K
	p2calc_i[mbl_il_I_mbl_il_K_ik]  					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_k_i>;				// mbl_il_I_mbl_il_K
	p2calc_i[mbl_il_I_mbl_il_K_WiWk]  				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_k_i>;				// mbl_il_I_mbl_il_K
	p2calc_i[mbl_il_I_mbl_il_gap_k_i]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2i<mbl_il_I_mbl_il_gap_k_i>;				// mbl_il_I_mbl_il_K
	p2calc_i[mbl_il_I_mbl_il_gap_k_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2i<mbl_il_I_mbl_il_gap_k_i>;				// mbl_il_I_mbl_il_K
	p2calc_i[mbl_il_I_mbl_il_gap_kWk_i]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2i<mbl_il_I_mbl_il_gap_kWk_i>;			// mbl_il_I_mbl_il_gap_Wk
	p2calc_i[mbl_il_I_mbl_il_gap_kWk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2i<mbl_il_I_mbl_il_gap_kWk_i>;			// mbl_il_I_mbl_il_gap_Wk
	p2calc_i[mbl_il_I_mbl_il_gap_Wk_ik]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_kWk_i>;
	p2calc_i[mbl_il_I_mbl_il_gap_Wk_Wi]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_kWk_i>;
	p2calc_i[mbl_il_pb_I_mbl_il_pb_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_i[mbl_il_gap_Wi_mbl_il_K_ik]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_gap_Wi_mbl_il_gap_k_i>;
	p2calc_i[mbl_il_gap_Wi_mbl_il_K_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_gap_Wi_mbl_il_gap_k_i>;
	p2calc_i[mbl_il_gap_Wi_mbl_il_gap_k_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2i<mbl_il_gap_Wi_mbl_il_gap_k_i>;
	p2calc_i[mbl_il_gap_Wi_mbl_il_gap_k_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2i<mbl_il_gap_Wi_mbl_il_gap_k_i>;
	p2calc_i[mbl_il_gap_iWi_mbl_il_K_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_gap_Wi_mbl_il_gap_k_i>;
	p2calc_i[mbl_il_gap_iWi_mbl_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_gap_Wi_mbl_il_gap_k_i>;
	p2calc_i[mbl_il_gap_i_mbl_il_K_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_k_i>;
	p2calc_i[mbl_il_gap_i_mbl_il_K_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_k_i>;
	p2calc_i[mbl_il_gap_i_mbl_il_gap_Wk_k]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_kWk_i>;
	p2calc_i[mbl_il_gap_i_mbl_il_gap_Wk_Wi]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2i<mbl_il_I_mbl_il_gap_kWk_i>;

	p2calc_k[hp_init]										= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_I_k>;		 // hp_init_align
	p2calc_k[hp_init_align_ik]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_I_k>;		 // hp_init_align
	p2calc_k[hp_init_align_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_I_k>;		 // hp_init_align
	p2calc_k[hp_init_gap_I_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_k<hp_init_gap_I_k>;		 // hp_init_align
	p2calc_k[hp_init_gap_I_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_k<hp_init_gap_I_k>;		 // hp_init_align
	p2calc_k[hp_init_gap_K_i]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_I_k>;		 // hp_init_align
	p2calc_k[hp_init_gap_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_I_k>;		 // hp_init_align
	p2calc_k[hp_init_gap_Wi_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_iWi_k>;   // hp_init_gap_Wi
	p2calc_k[hp_init_gap_Wi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_iWi_k>;   // hp_init_gap_Wi
	p2calc_k[hp_init_gap_Wk_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_iWk_k>;   // hp_init_gap_Wk
	p2calc_k[hp_init_gap_Wk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_iWk_k>;   // hp_init_gap_Wk
	p2calc_k[hp_init_gap_iWi_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_k<hp_init_gap_iWi_k>;  // hp_init_gap_Wi
	p2calc_k[hp_init_gap_iWi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_k<hp_init_gap_iWi_k>;  // hp_init_gap_Wi
	p2calc_k[hp_init_gap_kWk_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_iWk_k>;  // hp_init_gap_Wk
	p2calc_k[hp_init_gap_kWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_iWk_k>;  // hp_init_gap_Wk
	p2calc_k[hp_init_gap_iWk_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_k<hp_init_gap_iWk_k>;  // hp_init_gap_Wk
	p2calc_k[hp_init_gap_iWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_k<hp_init_gap_iWk_k>;  // hp_init_gap_Wk
	p2calc_k[hp_init_gap_kWi_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_iWi_k>;  // hp_init_gap_Wi
	p2calc_k[hp_init_gap_kWi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_k<hp_init_gap_iWi_k>;  // hp_init_gap_Wi
	p2calc_k[hp_pb_IK]    								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;		 // hp_init_align
	p2calc_k[stem_IK]										= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2k<bi_gap_i_bk_K>;					// bi_I_bk_K
	p2calc_k[stem_I_stem_gap_kWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2k<bi_gap_i_bk_gap_Wk>;			// bi_I_bk_gap_Wk
	p2calc_k[stem_no_mbl_I_stem_gap_kWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_k[stem_gap_iWi_stem_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2k<bi_gap_iWi_bk_K>;			// bi_gap_Wi_bk_K
	p2calc_k[stem_no_mbl_gap_iWi_stem_K]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_k[bi_I_bk_K]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2k<bi_gap_i_bk_K>;  				 // bi_I_bk_K
	p2calc_k[bi_I_bk_gap_k]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2k<bi_gap_i_bk_K>;  				 // bi_I_bk_K
	p2calc_k[bi_I_bk_gap_kWk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2k<bi_gap_i_bk_gap_Wk>;				 // bi_I_bk_K
	p2calc_k[bi_I_bk_gap_Wk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2k<bi_gap_i_bk_gap_Wk>;				 // bi_I_bk_K
	p2calc_k[bi_pb_I_bk_pb_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;			 // il_I_il_K
	p2calc_k[bi_gap_i_bk_K]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2k<bi_gap_i_bk_K>;  				 // bi_I_bk_K
	p2calc_k[bi_gap_i_bk_gap_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2k<bi_gap_i_bk_gap_Wk>;				 // bi_I_bk_K
	p2calc_k[bi_gap_iWi_bk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2k<bi_gap_iWi_bk_K>;				 // bi_I_bk_K
	p2calc_k[bi_gap_Wi_bk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2k<bi_gap_iWi_bk_K>;				 // bi_I_bk_K
	p2calc_k[bi_gap_Wi_bk_gap_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2k<bi_gap_iWi_bk_K>;				 // bi_I_bk_K
	p2calc_k[bWi_I_bWk_K]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2k<il_gap_i_il_K_k>;				// il_I_il_K
	p2calc_k[bWi_I_bWk_gap_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2k<il_gap_i_il_K_k>;			 // il_I_il_K
	p2calc_k[bWi_I_bWk_gap_kWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2k<il_gap_i_il_gap_Wk_k>; 			 // il_I_il_K
	p2calc_k[bWi_I_bWk_gap_Wk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2k<il_gap_i_il_gap_Wk_k>; 			 // il_I_il_K
	p2calc_k[bWi_pb_I_bWk_pb_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;			 // il_I_il_K
	p2calc_k[bWi_gap_i_bWk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkKgap2k<il_gap_i_il_K_k>;			 // il_I_il_K
	p2calc_k[bWi_gap_i_bWk_gap_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkKgap2k<il_gap_i_il_gap_Wk_k>; 			 // il_I_il_K
	p2calc_k[bWi_gap_iWi_bWk_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkKgap2k<il_gap_iWi_il_K_k>; 			 // il_I_il_K
	p2calc_k[bWi_gap_Wi_bWk_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2k<il_gap_iWi_il_K_k>; 			 // il_I_il_K
	p2calc_k[bWi_gap_Wi_bWk_gap_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiIbWkK2k<il_gap_iWi_il_K_k>; 			 // il_I_il_K
	p2calc_k[il_I_il_K_ik]								= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_i_il_K_k>;					// il_I_il_K
	p2calc_k[il_I_il_K_WiWk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_i_il_K_k>;					// il_I_il_K
	p2calc_k[il_I_il_gap_k_i]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_i_il_K_k>;					// il_I_il
	p2calc_k[il_I_il_gap_k_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_i_il_K_k>;					// il_I_il
	p2calc_k[il_I_il_gap_Wk_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_i_il_gap_Wk_k>;				// il_I_il_gap_Wk
	p2calc_k[il_I_il_gap_Wk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_i_il_gap_Wk_k>;				// il_I_il_gap_Wk
	p2calc_k[il_I_il_gap_kWk_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_i_il_gap_Wk_k>;				// il_I_il_gap_Wk
	p2calc_k[il_I_il_gap_kWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_i_il_gap_Wk_k>;				// il_I_il_gap_Wk
	p2calc_k[il_pb_I_il_pb_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;				// il_I_il_K
	p2calc_k[il_gap_iWi_il_K_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2k<il_gap_iWi_il_K_k>;				// il_gap_Wi_il_K
	p2calc_k[il_gap_iWi_il_K_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2k<il_gap_iWi_il_K_k>;				// il_gap_Wi_il_K
	p2calc_k[il_gap_i_il_K_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2k<il_gap_i_il_K_k>;					// il_I_il_K
	p2calc_k[il_gap_i_il_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2k<il_gap_i_il_K_k>;					// il_I_il_K
	p2calc_k[il_gap_i_il_gap_Wk_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2k<il_gap_i_il_gap_Wk_k>;				// il_I_il_K
	p2calc_k[il_gap_i_il_gap_Wk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2k<il_gap_i_il_gap_Wk_k>;				// il_I_il_K
	p2calc_k[il_gap_Wi_il_K_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_iWi_il_K_k>;				// il_gap_Wi_il_K
	p2calc_k[il_gap_Wi_il_K_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_iWi_il_K_k>;				// il_gap_Wi_il_K
	p2calc_k[il_gap_Wi_il_gap_k_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_iWi_il_K_k>;				// il_gap_Wi_il_K
	p2calc_k[il_gap_Wi_il_gap_k_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2k<il_gap_iWi_il_K_k>;				// il_gap_Wi_il_K
	p2calc_k[mblIK]  										= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_K_k>;						// mbl_bi_I_mbl_bk_K
	p2calc_k[mbl_bWi_I_mbl_bWk_K]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_K_k>;						// mbl_bi_I_mbl_bk_K
	p2calc_k[mbl_bWi_I_mbl_bWk_gap_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_gap_Wk_k>;						// mbl_bi_I_mbl_bk_K
	p2calc_k[mbl_bWi_gap_Wi_mbl_bWk_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_iWi_mbl_il_K_k>;						// mbl_bi_I_mbl_bk_K
	p2calc_k[mbl_il_I_mbl_il_K_ik]  					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_K_k>;				// mbl_il_I_mbl_il_K
	p2calc_k[mbl_il_I_mbl_il_K_WiWk]  				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_K_k>;				// mbl_il_I_mbl_il_K
	p2calc_k[mbl_il_I_mbl_il_gap_k_i]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_K_k>;				// mbl_il_I_mbl_il_K
	p2calc_k[mbl_il_I_mbl_il_gap_k_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_K_k>;				// mbl_il_I_mbl_il_K
	p2calc_k[mbl_il_I_mbl_il_gap_kWk_i]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_gap_Wk_k>;			// mbl_il_I_mbl_il_gap_Wk
	p2calc_k[mbl_il_I_mbl_il_gap_kWk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_gap_Wk_k>;			// mbl_il_I_mbl_il_gap_Wk
	p2calc_k[mbl_il_I_mbl_il_gap_Wk_ik]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_gap_Wk_k>;
	p2calc_k[mbl_il_I_mbl_il_gap_Wk_Wi]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_i_mbl_il_gap_Wk_k>;
	p2calc_k[mbl_il_pb_I_mbl_il_pb_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_k[mbl_il_gap_Wi_mbl_il_K_ik]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_iWi_mbl_il_K_k>;
	p2calc_k[mbl_il_gap_Wi_mbl_il_K_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_iWi_mbl_il_K_k>;
	p2calc_k[mbl_il_gap_Wi_mbl_il_gap_k_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_iWi_mbl_il_K_k>;
	p2calc_k[mbl_il_gap_Wi_mbl_il_gap_k_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2k<mbl_il_gap_iWi_mbl_il_K_k>;
	p2calc_k[mbl_il_gap_iWi_mbl_il_K_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2k<mbl_il_gap_iWi_mbl_il_K_k>;
	p2calc_k[mbl_il_gap_iWi_mbl_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2k<mbl_il_gap_iWi_mbl_il_K_k>;
	p2calc_k[mbl_il_gap_i_mbl_il_K_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2k<mbl_il_gap_i_mbl_il_K_k>;
	p2calc_k[mbl_il_gap_i_mbl_il_K_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2k<mbl_il_gap_i_mbl_il_K_k>;
	p2calc_k[mbl_il_gap_i_mbl_il_gap_Wk_k]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2k<mbl_il_gap_i_mbl_il_gap_Wk_k>;
	p2calc_k[mbl_il_gap_i_mbl_il_gap_Wk_Wi]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2k<mbl_il_gap_i_mbl_il_gap_Wk_k>;


	p2calc_Wi[hp_init] 									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_Wk_Wi>;  	// hp_init_align
	p2calc_Wi[hp_init_align_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_Wk_Wi>;
	p2calc_Wi[hp_init_align_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_Wk_Wi>;
	p2calc_Wi[hp_init_gap_I_k] 						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_iWk_Wi>;    // hp_init_gap_I
	p2calc_Wi[hp_init_gap_I_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_iWk_Wi>;    // hp_init_gap_I
	p2calc_Wi[hp_init_gap_K_i] 						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_kWk_Wi>;    // hp_init_gap_K
	p2calc_Wi[hp_init_gap_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_kWk_Wi>;    // hp_init_gap_K
	p2calc_Wi[hp_init_gap_Wi_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_Wk_Wi>;
	p2calc_Wi[hp_init_gap_Wi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_Wk_Wi>;
	p2calc_Wi[hp_init_gap_Wk_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_K_Wi<hp_init_gap_Wk_Wi>;
	p2calc_Wi[hp_init_gap_Wk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_K_Wi<hp_init_gap_Wk_Wi>;
	p2calc_Wi[hp_init_gap_iWi_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_iWk_Wi>;
	p2calc_Wi[hp_init_gap_iWi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_iWk_Wi>;
	p2calc_Wi[hp_init_gap_kWk_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_K_Wi<hp_init_gap_kWk_Wi>;
	p2calc_Wi[hp_init_gap_kWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_K_Wi<hp_init_gap_kWk_Wi>;
	p2calc_Wi[hp_init_gap_iWk_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_K_Wi<hp_init_gap_iWk_Wi>;
	p2calc_Wi[hp_init_gap_iWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_K_Wi<hp_init_gap_iWk_Wi>;
	p2calc_Wi[hp_init_gap_kWi_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_kWk_Wi>;
	p2calc_Wi[hp_init_gap_kWi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wi<hp_init_gap_kWk_Wi>;
	p2calc_Wi[hp_pb_IK]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wi[stem_IK] 									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wi<bWi_I_bWk_gap_Wk>;				// bWi_I_bWk_K
	p2calc_Wi[stem_I_stem_gap_kWk] 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2Wi<bWi_I_bWk_gap_kWk>;				// bWi_I_bWk_K
	p2calc_Wi[stem_no_mbl_I_stem_gap_kWk]		 	= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wi[stem_gap_iWi_stem_K] 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wi<bWi_gap_i_bWk_gap_Wk>;				// bWi_I_bWk_K
	p2calc_Wi[stem_no_mbl_gap_iWi_stem_K] 			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wi[bi_I_bk_K]				 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wi<il_I_il_gap_Wk_Wi>;
	p2calc_Wi[bi_I_bk_gap_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wi<il_I_il_gap_kWk_Wi>;
	p2calc_Wi[bi_I_bk_gap_kWk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkKgap2Wi<il_I_il_gap_kWk_Wi>;
	p2calc_Wi[bi_I_bk_gap_Wk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkKgap2Wi<il_I_il_gap_Wk_Wi>;
	p2calc_Wi[bi_pb_I_bk_pb_K]  						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wi[bi_gap_i_bk_K] 		 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wi<il_gap_i_il_gap_Wk_Wi>;
	p2calc_Wi[bi_gap_i_bk_gap_Wk] 		 			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkKgap2Wi<il_gap_i_il_gap_Wk_Wi>;
	p2calc_Wi[bi_gap_iWi_bk_K] 		 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wi<il_gap_i_il_gap_Wk_Wi>;
	p2calc_Wi[bi_gap_Wi_bk_K] 		 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wi<il_I_il_gap_Wk_Wi>;
	p2calc_Wi[bi_gap_Wi_bk_gap_k]	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wi<il_I_il_gap_kWk_Wi>;
	p2calc_Wi[bWi_I_bWk_K]  							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wi<bWi_I_bWk_gap_Wk>;
	p2calc_Wi[bWi_I_bWk_gap_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wi<bWi_I_bWk_gap_kWk>;
	p2calc_Wi[bWi_I_bWk_gap_kWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2Wi<bWi_I_bWk_gap_kWk>;
	p2calc_Wi[bWi_I_bWk_gap_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2Wi<bWi_I_bWk_gap_Wk>;
	p2calc_Wi[bWi_pb_I_bWk_pb_K]  					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wi[bWi_gap_iWi_bWk_K] 	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wi<bWi_gap_i_bWk_gap_Wk>;
	p2calc_Wi[bWi_gap_i_bWk_K] 	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wi<bWi_gap_i_bWk_gap_Wk>;
	p2calc_Wi[bWi_gap_i_bWk_gap_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2Wi<bWi_gap_i_bWk_gap_Wk>;
	p2calc_Wi[bWi_gap_Wi_bWk_K] 	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wi<bWi_I_bWk_gap_Wk>;
	p2calc_Wi[bWi_gap_Wi_bWk_gap_k] 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wi<bWi_I_bWk_gap_kWk>;
	p2calc_Wi[il_I_il_K_ik]  							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_I_il_gap_Wk_Wi>;
	p2calc_Wi[il_I_il_K_WiWk]  						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_I_il_gap_Wk_Wi>;
	p2calc_Wi[il_I_il_gap_k_i]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_I_il_gap_kWk_Wi>;
	p2calc_Wi[il_I_il_gap_k_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_I_il_gap_kWk_Wi>;
	p2calc_Wi[il_I_il_gap_Wk_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wi<il_I_il_gap_Wk_Wi>;
	p2calc_Wi[il_I_il_gap_Wk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wi<il_I_il_gap_Wk_Wi>;
	p2calc_Wi[il_I_il_gap_kWk_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wi<il_I_il_gap_kWk_Wi>;
	p2calc_Wi[il_I_il_gap_kWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wi<il_I_il_gap_kWk_Wi>;
	p2calc_Wi[il_pb_I_il_pb_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wi[il_gap_iWi_il_K_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_gap_i_il_gap_Wk_Wi>;
	p2calc_Wi[il_gap_iWi_il_K_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_gap_i_il_gap_Wk_Wi>;
	p2calc_Wi[il_gap_i_il_K_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_gap_i_il_gap_Wk_Wi>;
	p2calc_Wi[il_gap_i_il_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_gap_i_il_gap_Wk_Wi>;
	p2calc_Wi[il_gap_i_il_gap_Wk_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wi<il_gap_i_il_gap_Wk_Wi>;
	p2calc_Wi[il_gap_i_il_gap_Wk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wi<il_gap_i_il_gap_Wk_Wi>;
	p2calc_Wi[il_gap_Wi_il_K_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_I_il_gap_Wk_Wi>;
	p2calc_Wi[il_gap_Wi_il_K_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_I_il_gap_Wk_Wi>;
	p2calc_Wi[il_gap_Wi_il_gap_k_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_I_il_gap_kWk_Wi>;
	p2calc_Wi[il_gap_Wi_il_gap_k_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wi<il_I_il_gap_kWk_Wi>;
	p2calc_Wi[mblIK]										= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_bWi_I_mbl_bWk_gap_Wk>;				 // mbl_bWi_I_mbl_bWk_K
	p2calc_Wi[mbl_bWi_I_mbl_bWk_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_bWi_I_mbl_bWk_gap_Wk>;
	p2calc_Wi[mbl_bWi_I_mbl_bWk_gap_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wi<mbl_bWi_I_mbl_bWk_gap_Wk>;
	p2calc_Wi[mbl_bWi_gap_Wi_mbl_bWk_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_bWi_I_mbl_bWk_gap_Wk>;
	p2calc_Wi[mbl_il_I_mbl_il_K_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_I_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_I_mbl_il_K_WiWk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_I_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_I_mbl_il_gap_k_i]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_I_mbl_il_gap_kWk_Wi>;
	p2calc_Wi[mbl_il_I_mbl_il_gap_k_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_I_mbl_il_gap_kWk_Wi>;
	p2calc_Wi[mbl_il_I_mbl_il_gap_kWk_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wi<mbl_il_I_mbl_il_gap_kWk_Wi>;
	p2calc_Wi[mbl_il_I_mbl_il_gap_kWk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wi<mbl_il_I_mbl_il_gap_kWk_Wi>;
	p2calc_Wi[mbl_il_I_mbl_il_gap_Wk_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wi<mbl_il_I_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_I_mbl_il_gap_Wk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wi<mbl_il_I_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_pb_I_mbl_il_pb_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wi[mbl_il_gap_Wi_mbl_il_K_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_I_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_gap_Wi_mbl_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_I_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_gap_Wi_mbl_il_gap_k_i]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_I_mbl_il_gap_kWk_Wi>;
	p2calc_Wi[mbl_il_gap_Wi_mbl_il_gap_k_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_I_mbl_il_gap_kWk_Wi>;
	p2calc_Wi[mbl_il_gap_iWi_mbl_il_K_k]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_gap_i_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_gap_iWi_mbl_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_gap_i_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_gap_i_mbl_il_K_k]  			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_gap_i_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_gap_i_mbl_il_K_WiWk]  		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wi<mbl_il_gap_i_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_gap_i_mbl_il_gap_Wk_k]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wi<mbl_il_gap_i_mbl_il_gap_Wk_Wi>;
	p2calc_Wi[mbl_il_gap_i_mbl_il_gap_Wk_Wi]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wi<mbl_il_gap_i_mbl_il_gap_Wk_Wi>;

	p2calc_Wk[hp_init] 									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_Wi_Wk>;  	// hp_init_align
	p2calc_Wk[hp_init_align_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_Wi_Wk>;
	p2calc_Wk[hp_init_align_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_Wi_Wk>;
	p2calc_Wk[hp_init_gap_I_k] 						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_iWi_Wk>;    // hp_init_gap_I
	p2calc_Wk[hp_init_gap_I_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_iWi_Wk>;    // hp_init_gap_I
	p2calc_Wk[hp_init_gap_K_i] 						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_kWi_Wk>;    // hp_init_gap_K
	p2calc_Wk[hp_init_gap_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_kWi_Wk>;    // hp_init_gap_K
	p2calc_Wk[hp_init_gap_Wi_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_Wk<hp_init_gap_Wi_Wk>;
	p2calc_Wk[hp_init_gap_Wi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_Wk<hp_init_gap_Wi_Wk>;
	p2calc_Wk[hp_init_gap_Wk_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_Wi_Wk>;
	p2calc_Wk[hp_init_gap_Wk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_Wi_Wk>;
	p2calc_Wk[hp_init_gap_iWi_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_Wk<hp_init_gap_iWi_Wk>;
	p2calc_Wk[hp_init_gap_iWi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_Wk<hp_init_gap_iWi_Wk>;
	p2calc_Wk[hp_init_gap_kWk_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_kWi_Wk>;
	p2calc_Wk[hp_init_gap_kWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_kWi_Wk>;
	p2calc_Wk[hp_init_gap_iWk_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_iWi_Wk>;
	p2calc_Wk[hp_init_gap_iWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_Wk<hp_init_gap_iWi_Wk>;
	p2calc_Wk[hp_init_gap_kWi_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_Wk<hp_init_gap_kWi_Wk>;
	p2calc_Wk[hp_init_gap_kWi_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template hp_align_gap_I_Wk<hp_init_gap_kWi_Wk>;
	p2calc_Wk[hp_pb_IK]									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wk[stem_IK] 									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wk<bWi_gap_Wi_bWk_K>;				// bWi_I_bWk_K
	p2calc_Wk[stem_I_stem_gap_kWk] 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wk<bWi_gap_Wi_bWk_gap_k>;				// bWi_I_bWk_K
	p2calc_Wk[stem_no_mbl_I_stem_gap_kWk]		 	= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wk[stem_gap_iWi_stem_K] 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2Wk<bWi_gap_iWi_bWk_K>;				// bWi_I_bWk_K
	p2calc_Wk[stem_no_mbl_gap_iWi_stem_K]	 		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wk[bi_I_bk_K]				 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wk<il_gap_Wi_il_K_Wk>;
	p2calc_Wk[bi_I_bk_gap_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wk<il_gap_Wi_il_gap_k_Wk>;
	p2calc_Wk[bi_I_bk_gap_kWk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wk<il_gap_Wi_il_gap_k_Wk>;
	p2calc_Wk[bi_I_bk_gap_Wk]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wk<il_gap_Wi_il_K_Wk>;
	p2calc_Wk[bi_pb_I_bk_pb_K]  						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wk[bi_gap_i_bk_K] 		 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wk<il_gap_iWi_il_K_Wk>;
	p2calc_Wk[bi_gap_i_bk_gap_Wk] 		 			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wk<il_gap_iWi_il_K_Wk>;
	p2calc_Wk[bi_gap_iWi_bk_K] 		 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkK2Wk<il_gap_iWi_il_K_Wk>;
	p2calc_Wk[bi_gap_Wi_bk_K] 		 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkKgap2Wk<il_gap_Wi_il_K_Wk>;
	p2calc_Wk[bi_gap_Wi_bk_gap_k]	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template biIbkKgap2Wk<il_gap_Wi_il_gap_k_Wk>;
	p2calc_Wk[bWi_I_bWk_K]  							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wk<bWi_gap_Wi_bWk_K>;
	p2calc_Wk[bWi_I_bWk_gap_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wk<bWi_gap_Wi_bWk_gap_k>;
	p2calc_Wk[bWi_I_bWk_gap_kWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wk<bWi_gap_Wi_bWk_gap_k>;
	p2calc_Wk[bWi_I_bWk_gap_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wk<bWi_gap_Wi_bWk_K>;
	p2calc_Wk[bWi_pb_I_bWk_pb_K]  					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wk[bWi_gap_iWi_bWk_K] 	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wk<bWi_gap_iWi_bWk_K>;
	p2calc_Wk[bWi_gap_i_bWk_K] 	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wk<bWi_gap_iWi_bWk_K>;
	p2calc_Wk[bWi_gap_i_bWk_gap_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2Wk<bWi_gap_iWi_bWk_K>;
	p2calc_Wk[bWi_gap_Wi_bWk_K] 	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2Wk<bWi_gap_Wi_bWk_K>;
	p2calc_Wk[bWi_gap_Wi_bWk_gap_k] 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIKgap2Wk<bWi_gap_Wi_bWk_gap_k>;
	p2calc_Wk[il_I_il_K_ik]  							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_Wi_il_K_Wk>;
	p2calc_Wk[il_I_il_K_WiWk]  						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_Wi_il_K_Wk>;
	p2calc_Wk[il_I_il_gap_k_i]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_Wi_il_gap_k_Wk>;
	p2calc_Wk[il_I_il_gap_k_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_Wi_il_gap_k_Wk>;
	p2calc_Wk[il_I_il_gap_Wk_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_Wi_il_K_Wk>;
	p2calc_Wk[il_I_il_gap_Wk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_Wi_il_K_Wk>;
	p2calc_Wk[il_I_il_gap_kWk_i]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_Wi_il_gap_k_Wk>;
	p2calc_Wk[il_I_il_gap_kWk_Wi]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_Wi_il_gap_k_Wk>;
	p2calc_Wk[il_pb_I_il_pb_K]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wk[il_gap_iWi_il_K_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wk<il_gap_iWi_il_K_Wk>;
	p2calc_Wk[il_gap_iWi_il_K_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wk<il_gap_iWi_il_K_Wk>;
	p2calc_Wk[il_gap_i_il_K_k]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_iWi_il_K_Wk>;
	p2calc_Wk[il_gap_i_il_K_WiWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_iWi_il_K_Wk>;
	p2calc_Wk[il_gap_i_il_gap_Wk_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_iWi_il_K_Wk>;
	p2calc_Wk[il_gap_i_il_gap_Wk_Wi]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilK2Wk<il_gap_iWi_il_K_Wk>;
	p2calc_Wk[il_gap_Wi_il_K_ik]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wk<il_gap_Wi_il_K_Wk>;
	p2calc_Wk[il_gap_Wi_il_K_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wk<il_gap_Wi_il_K_Wk>;
	p2calc_Wk[il_gap_Wi_il_gap_k_i]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wk<il_gap_Wi_il_gap_k_Wk>;
	p2calc_Wk[il_gap_Wi_il_gap_k_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template ilIilKgap2Wk<il_gap_Wi_il_gap_k_Wk>;
	p2calc_Wk[mblIK]										= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_bWi_gap_Wi_mbl_bWk_K>;				 // mbl_bWi_I_mbl_bWk_K
	p2calc_Wk[mbl_bWi_I_mbl_bWk_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_bWi_gap_Wi_mbl_bWk_K>;
	p2calc_Wk[mbl_bWi_I_mbl_bWk_gap_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_bWi_gap_Wi_mbl_bWk_K>;
	p2calc_Wk[mbl_bWi_gap_Wi_mbl_bWk_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wk<mbl_bWi_gap_Wi_mbl_bWk_K>;
	p2calc_Wk[mbl_il_I_mbl_il_K_ik]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_Wi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_I_mbl_il_K_WiWk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_Wi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_I_mbl_il_gap_k_i]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_Wi_mbl_il_gap_k_Wk>;
	p2calc_Wk[mbl_il_I_mbl_il_gap_k_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_Wi_mbl_il_gap_k_Wk>;
	p2calc_Wk[mbl_il_I_mbl_il_gap_kWk_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_Wi_mbl_il_gap_k_Wk>;
	p2calc_Wk[mbl_il_I_mbl_il_gap_kWk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_Wi_mbl_il_gap_k_Wk>;
	p2calc_Wk[mbl_il_I_mbl_il_gap_Wk_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_Wi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_I_mbl_il_gap_Wk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_Wi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_pb_I_mbl_il_pb_K]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::noCalc;
	p2calc_Wk[mbl_il_gap_Wi_mbl_il_K_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wk<mbl_il_gap_Wi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_gap_Wi_mbl_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wk<mbl_il_gap_Wi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_gap_Wi_mbl_il_gap_k_i]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wk<mbl_il_gap_Wi_mbl_il_gap_k_Wk>;
	p2calc_Wk[mbl_il_gap_Wi_mbl_il_gap_k_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wk<mbl_il_gap_Wi_mbl_il_gap_k_Wk>;
	p2calc_Wk[mbl_il_gap_iWi_mbl_il_K_k]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wk<mbl_il_gap_iWi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_gap_iWi_mbl_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIKgap2Wk<mbl_il_gap_iWi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_gap_i_mbl_il_K_k]  			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_iWi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_gap_i_mbl_il_K_WiWk]  		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_iWi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_gap_i_mbl_il_gap_Wk_k]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_iWi_mbl_il_K_Wk>;
	p2calc_Wk[mbl_il_gap_i_mbl_il_gap_Wk_Wi]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template mIK2Wk<mbl_il_gap_iWi_mbl_il_K_Wk>;


	if (!nobranch) {
		p2calc_mbl[stem_IK] 									= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2mbl<mblIK>;
		p2calc_mbl[stem_I_stem_gap_kWk] 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2mbl<mblIK>;
		p2calc_mbl[stem_gap_iWi_stem_K] 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2mbl<mblIK>;
		p2calc_mbl[bWi_I_bWk_K]  							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiWkIK2mbl<mblIK>;
		p2calc_mbl[bWi_I_bWk_gap_k]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiWkIK2mbl<mblIK>;
		p2calc_mbl[bWi_I_bWk_gap_kWk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiWkIK2mbl<mblIK>;
		p2calc_mbl[bWi_I_bWk_gap_Wk]						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiWkIK2mbl<mblIK>;
		p2calc_mbl[bWi_gap_iWi_bWk_K] 	 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiWkIK2mbl<mblIK>;
		p2calc_mbl[bWi_gap_i_bWk_K] 	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiWkIK2mbl<mblIK>;
		p2calc_mbl[bWi_gap_i_bWk_gap_Wk]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiWkIK2mbl<mblIK>;
		p2calc_mbl[bWi_gap_Wi_bWk_K] 	 					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiWkIK2mbl<mblIK>;
		p2calc_mbl[bWi_gap_Wi_bWk_gap_k] 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template bWiWkIK2mbl<mblIK>;
		p2calc_mbl[mblIK]										= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2mbl<mblIK>;
		p2calc_mbl[mbl_bWi_I_mbl_bWk_K]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2mbl<mblIK>;
		p2calc_mbl[mbl_bWi_I_mbl_bWk_gap_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2mbl<mblIK>;
		p2calc_mbl[mbl_bWi_gap_Wi_mbl_bWk_K]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::template stemIK2mbl<mblIK>;
	}

	p2end[hp_init] 						= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_align_ik]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_align_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_I_k] 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_I_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_K_i] 				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_K_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_Wi_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_Wi_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_Wk_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_Wk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_iWi_k]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_iWi_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_kWk_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_kWk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_iWk_k]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_iWk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_kWi_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[hp_init_gap_kWi_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::hp2end;
	p2end[stem_IK]							= &foldK< FOLDK_TEMPLATE_PARAMETERS >::stem2end;
	p2end[stem_I_stem_gap_kWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::stem2end;
	p2end[stem_no_mbl_I_stem_gap_kWk]= &foldK< FOLDK_TEMPLATE_PARAMETERS >::stem2end;
	p2end[stem_gap_iWi_stem_K]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::stem2end;
	p2end[stem_no_mbl_gap_iWi_stem_K]= &foldK< FOLDK_TEMPLATE_PARAMETERS >::stem2end;
	p2end[bi_I_bk_K]				 		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_ik2end;
	p2end[bi_I_bk_gap_k]					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_ik2end;
	p2end[bi_I_bk_gap_kWk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_ik2end;
	p2end[bi_I_bk_gap_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_ik2end;
	p2end[bi_gap_i_bk_K] 		 		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_ik2end;
	p2end[bi_gap_i_bk_gap_Wk] 		 	= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_ik2end;
	p2end[bi_gap_iWi_bk_K] 		 		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_ik2end;
	p2end[bi_gap_Wi_bk_K] 		 		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_ik2end;
	p2end[bi_gap_Wi_bk_gap_k]	 		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_ik2end;
	p2end[bWi_I_bWk_K]  					= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_WiWk2end;
	p2end[bWi_I_bWk_gap_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_WiWk2end;
	p2end[bWi_I_bWk_gap_kWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_WiWk2end;
	p2end[bWi_I_bWk_gap_Wk]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_WiWk2end;
	p2end[bWi_gap_iWi_bWk_K] 	 		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_WiWk2end;
	p2end[bWi_gap_i_bWk_K] 	 			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_WiWk2end;
	p2end[bWi_gap_i_bWk_gap_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_WiWk2end;
	p2end[bWi_gap_Wi_bWk_K] 	 		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_WiWk2end;
	p2end[bWi_gap_Wi_bWk_gap_k] 		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::bulge_WiWk2end;
	p2end[il_I_il_K_ik]  				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_I_il_K_WiWk]  				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_I_il_gap_k_i]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_I_il_gap_k_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_I_il_gap_Wk_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_I_il_gap_Wk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_I_il_gap_kWk_i]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_I_il_gap_kWk_Wi]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_gap_iWi_il_K_k]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_gap_iWi_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_gap_i_il_K_k]				= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_gap_i_il_K_WiWk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_gap_i_il_gap_Wk_k]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_gap_i_il_gap_Wk_Wi]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_gap_Wi_il_K_ik]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_gap_Wi_il_K_Wk]			= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_gap_Wi_il_gap_k_i]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
	p2end[il_gap_Wi_il_gap_k_Wk]		= &foldK< FOLDK_TEMPLATE_PARAMETERS >::ilIK2end;
}

#endif /*FOLDK */
