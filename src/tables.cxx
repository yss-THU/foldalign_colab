#ifndef TABLES
#define TABLES

#include <iostream>
#include "foldalign.hxx"
#include "helper.cxx"

template< bool realigning, bool mblrealign >
class tables {

public:
	tables();
	
	bool get_right_branch(const stateType state) const {
		return branch_right[state];
	};

	bool get_right_store(const stateType state) const {
		return right_store[state];
	};

	bool isThisASeedConstraint(const stateType state) const {
		return seedConstraint[state];
	};
	
	stateType convertState2seed(const stateType state) const {
		return state2seed[state];
	};
	
	stateType convertSeed2state(const stateType state) const {
		return seed2state[state];
	}
	
private:
	bool branch_right[flow_size];
	bool right_store[flow_size];
	bool seedConstraint[flow_size];

	stateType state2seed[flow_size];
	stateType seed2state[flow_size];

};


template< bool realigning, bool mblrealign >
inline tables<realigning, mblrealign>::tables() {

	helper::init_array(state2seed, flow_size, noState);
	state2seed[stem_IK] 					= seed_bp;
	state2seed[stem_I_stem_gap_kWk] 		= seed_bpI;
	state2seed[stem_gap_iWi_stem_K] 		= seed_bpK;
/*	state2seed[bWi_I_bWk_K] 				= seed_bWi_I_bWk_K;
	state2seed[bWi_I_bWk_gap_k] 			= seed_bWi_I_bWk_gap_k;
	state2seed[bWi_I_bWk_gap_kWk] 			= seed_bWi_I_bWk_gap_kWk;
	state2seed[bWi_I_bWk_gap_Wk] 			= seed_bWi_I_bWk_gap_Wk;
	state2seed[bWi_gap_iWi_bWk_K] 			= seed_bWi_gap_iWi_bWk_K;
	state2seed[bWi_gap_i_bWk_K] 			= seed_bWi_gap_i_bWk_K;
	state2seed[bWi_gap_i_bWk_gap_Wk] 		= seed_bWi_gap_i_bWk_gap_Wk;
	state2seed[bWi_gap_Wi_bWk_K] 			= seed_bWi_gap_Wi_bWk_K;
	state2seed[bWi_gap_Wi_bWk_gap_k] 		= seed_bWi_gap_Wi_bWk_gap_k;*/
//	state2seed[mblIK] 						= seed_mblIK;
/*	state2seed[mbl_bWi_I_mbl_bWk_K] 		= seed_mbl_bWi_I_mbl_bWk_K;
	state2seed[mbl_bWi_I_mbl_bWk_gap_Wk] 	= seed_mbl_bWi_I_mbl_bWk_gap_Wk;
	state2seed[mbl_bWi_gap_Wi_mbl_bWk_K] 	= seed_mbl_bWi_gap_Wi_mbl_bWk_K;
*/
	helper::init_array(seedConstraint, flow_size, false);
	for(stateType p=0; p< flow_size; p++) {
		if (state2seed[p] != noState) {
			seedConstraint[state2seed[p]] = true;
		}
	}

//helper::printArray(seedConstraint, flow_size);
//	seedConstraint[seed_bp] = true;
//	seedConstraint[seed_bpI] = true;
//	seedConstraint[seed_bpK] = true;

	helper::init_array(seed2state, flow_size, noState);
	for(stateType p = 0; p < flow_size; p++) {
		if (state2seed[p] != noState) {
			seed2state[state2seed[p]] = p;
		}

		if (seed2state[p] == noState) {
			seed2state[p] = p;
		}
//std::cout << "State: " << int(p) << " seed2state: " << int(seed2state[p]) << " state2seed: " << int(state2seed[p]) << std::endl;
	}

	helper::init_array(branch_right, flow_size, false);
	branch_right[stem_IK] = true;
	branch_right[stem_I_stem_gap_kWk] = true;
	branch_right[stem_gap_iWi_stem_K] = true;
	branch_right[seed_bp] = true;
	branch_right[seed_bpI] = true;
	branch_right[seed_bpK] = true;

	// The alignments stored in the long term memory are either the same as those
	// which can be the right part of a multiloop (during scan), or all alignments
	// (during global or realigment). This is controled by the right_store array.
	if ( realigning ) {
		// During the "real" backtrack all states are saved.
		helper::init_array(right_store, flow_size, true);
	}
	else if ( mblrealign ) {
		// During mblrealignment and global all the states which can be involved
		// in multibranch points must be stored. Otherwise it is not possible to
		// relocate the two parts of the mbl.
		helper::copyArray(right_store, branch_right, flow_size);
		right_store[mblIK] = true;
	}
	else{
		// The scan case where only the stem states needs to be stored.
		helper::copyArray(right_store, branch_right, flow_size);
	}

//	helper::init_array(seedConstraint, flow_size, false);
//	seedConstraint[seed_bp] = true;
//	seedConstraint[seed_bpI] = true;
//	seedConstraint[seed_bpK] = true;


}
#endif /* ARRAYS */
