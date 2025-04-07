#ifndef FOLDALIGNH
#define FOLDALIGNH

// This file defines the subrevison number. The file is made by the make file
// which tries to get the subversion revision number.
#include "revnumber.hxx"

#include <string>
/******************************************************************************
*                                                                             *
*   Copyright 2004 - 2007 Jakob Hull Havgaard, hull@bioinf.ku.dk              *
*                                                                             *
*   This file is part of Foldalign                                            *
*                                                                             *
*   Foldalign is free software; you can redistribute it and/or modify         *
*   it under the terms of the GNU General Public License as published by      *
*   the Free Software Foundation; either version 2 of the License, or         *
*   (at your option) any later version.                                       *
*                                                                             *
*   Foldalign is distributed in the hope that it will be useful,              *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with Foldalign; if not, write to the Free Software                  *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA *
*                                                                             *
******************************************************************************/


//*********************************************************
//                                                        *
// NOTE: Default parameters are set in arguments.cxx file *
//                                                        *
//*********************************************************

const std::string program_name = "FOLDALIGN";
const std::string major_version = "v2.5.3";
const std::string curr_revision = (rev == "") ? "final" : rev; // rev is defined in revnumber.hxx. If it is only an hexadecimal number, try running `git fetch --tags`. Check git-describe manual to understand what this number means
const std::string version = major_version + (curr_revision == "final" ? "" : " (" + curr_revision + ")"); //append (subversion) if not final
const std::string description = "Copyright by Jakob Hull Havgaard, 2004 - 2021\nReference: Sundfeld et al. Bioinformatics 2016, 32(8):1238-1240";

// Reference_size is the number of lines in the reference array
const int reference_size = 5;
const std::string reference[reference_size] = {
	"D. Sundfeld, J.H. Havgaard, A.C. de Melo, and J. Gorodkin",
	"Foldalign 2.5: multithreaded implementation for pairwise",
	"structural RNA alignment.",
	"Bioinformatics 2016, 32(8):1238-1240,",
	"doi: 10.1093/bioinformatics/btv748"
};


/***********************************************/
// By default foldalign is compiled to run on multicore
// systems. Uncomment the define SINGLECORE line below to make
// it compile as single core. This makes the program some 
// what faster.

//#define SINGLECORE



/***********************************************/
// Define the main variable types of the program
// In some cases there are a high and a low memory choice.

	// The relatively safe high memory variable types
	// The type of alignment scores
//	typedef int     			scoreType;	// High memory choice - larger max scores
	typedef int				scoreType;  // Low memory and less safe choice

	// The state type. The states control the program flow and are also used as
	// backtrack pointers (not c++ pointers but dynamical programming pointers).
	typedef unsigned char	stateType;  		

	// Position along a sequence. This type must be able to handle the last
	// position in the longest sequence with out overflowing
	typedef long		         positionType; // More memory longer sequences
//	typedef short				positionType; // Less memory shorter sequences

	// Lenght type. Used for lambda, delta, window sizes etc.
//	typedef short	         lengthType;		   
	typedef int	         lengthType;		   

	/*********************************************/
	// Constants


	// The number of bytes in a char. (The number has been changed to fit the data)
	const long char_size = 16;
	// The memory scale (Changed to fit the data)
	const long mem_scale = 10485760; // mb

	// The main constants

	// The BIG negative score
	const scoreType 			big_neg = -20000;
	
	const scoreType warn_low  = scoreType( 0.8 * big_neg);
	const scoreType warn_high = scoreType(-0.8 * big_neg);

	// The number of different states +1 (state zero is used for the noState)
	static const stateType flow_size = 104;

	// An alignment must have state equal or greater than this before it is
	// allowed as valid structural alignment. Any state below this cut is should
	// be a hairpin state (and no hairpin state should have value >= than this)
	// If set to noState then all alignments are printed.
	static const stateType min_struc_state = 21;

	// This is the default no state which set wrongly (!= 0) will mess up the program. 
	const stateType noState                         =  0;
	// The states. They control the flow of the program.
	const stateType hp_init      							=  1;
	const stateType hp_init_align_ik						=  2;
	const stateType hp_init_align_WiWk					=  3;
	const stateType hp_init_gap_I_k						=  4;
	const stateType hp_init_gap_I_WiWk					=  5;
	const stateType hp_init_gap_K_i						=  6;
	const stateType hp_init_gap_K_WiWk					=  7;
	const stateType hp_init_gap_Wi_ik					=  8;
	const stateType hp_init_gap_Wi_Wk					=  9;
	const stateType hp_init_gap_Wk_ik					= 10;
	const stateType hp_init_gap_Wk_Wi					= 11;
	const stateType hp_init_gap_iWi_k					= 12;
	const stateType hp_init_gap_iWi_Wk					= 13;
	const stateType hp_init_gap_kWk_i					= 14;
	const stateType hp_init_gap_kWk_Wi					= 15;
	const stateType hp_init_gap_iWk_k					= 16;
	const stateType hp_init_gap_iWk_Wi					= 17;
	const stateType hp_init_gap_kWi_i					= 18;
	const stateType hp_init_gap_kWi_Wk					= 19;
	const stateType hp_pb_IK								= 20;
	const stateType stem_IK									= 21;
	const stateType stem_I_stem_gap_kWk					= 22;
	const stateType stem_no_mbl_I_stem_gap_kWk		= 23;
	const stateType stem_gap_iWi_stem_K					= 24;
	const stateType stem_no_mbl_gap_iWi_stem_K		= 25;
	const stateType bi_I_bk_K								= 26;
	const stateType bi_I_bk_gap_k							= 27;
	const stateType bi_I_bk_gap_kWk						= 28;
	const stateType bi_I_bk_gap_Wk						= 29;
	const stateType bi_pb_I_bk_pb_K						= 30;
	const stateType bi_gap_i_bk_K							= 31;
	const stateType bi_gap_i_bk_gap_Wk					= 32;
	const stateType bi_gap_iWi_bk_K						= 33;
	const stateType bi_gap_Wi_bk_K						= 34;
	const stateType bi_gap_Wi_bk_gap_k					= 35;
	const stateType bWi_I_bWk_K							= 36;
	const stateType bWi_I_bWk_gap_k						= 37;
	const stateType bWi_I_bWk_gap_kWk					= 38;
	const stateType bWi_I_bWk_gap_Wk						= 39;
	const stateType bWi_pb_I_bWk_pb_K					= 40;
	const stateType bWi_gap_i_bWk_K						= 41;
	const stateType bWi_gap_i_bWk_gap_Wk				= 42;
	const stateType bWi_gap_iWi_bWk_K					= 43;
	const stateType bWi_gap_Wi_bWk_K						= 44;
	const stateType bWi_gap_Wi_bWk_gap_k				= 45;
	const stateType il_I_il_K_ik							= 46;
	const stateType il_I_il_K_WiWk						= 47;
	const stateType il_I_il_gap_k_i						= 48;
	const stateType il_I_il_gap_k_WiWk					= 49;
	const stateType il_I_il_gap_Wk_ik					= 50;
	const stateType il_I_il_gap_Wk_Wi					= 51;
	const stateType il_I_il_gap_kWk_i					= 52;
	const stateType il_I_il_gap_kWk_Wi					= 53;
	const stateType il_pb_I_il_pb_K						= 54;
	const stateType il_gap_iWi_il_K_k					= 55;
	const stateType il_gap_iWi_il_K_Wk					= 56;
	const stateType il_gap_i_il_K_k						= 57;
	const stateType il_gap_i_il_K_WiWk					= 58;
	const stateType il_gap_i_il_gap_Wk_k				= 59;
	const stateType il_gap_i_il_gap_Wk_Wi				= 60;
	const stateType il_gap_Wi_il_K_ik					= 61;
	const stateType il_gap_Wi_il_K_Wk					= 62;
	const stateType il_gap_Wi_il_gap_k_i				= 63;
	const stateType il_gap_Wi_il_gap_k_Wk				= 64;
	const stateType mblIK									= 65;
	const stateType mbl_bWi_I_mbl_bWk_K					= 66;
	const stateType mbl_bWi_I_mbl_bWk_gap_Wk			= 67;
	const stateType mbl_bWi_gap_Wi_mbl_bWk_K			= 68;
	const stateType mbl_il_I_mbl_il_K_ik				= 69;
	const stateType mbl_il_I_mbl_il_K_WiWk				= 70;
	const stateType mbl_il_I_mbl_il_gap_k_i			= 71;
	const stateType mbl_il_I_mbl_il_gap_k_WiWk		= 72;
	const stateType mbl_il_I_mbl_il_gap_kWk_i			= 73;
	const stateType mbl_il_I_mbl_il_gap_kWk_Wi		= 74;
	const stateType mbl_il_I_mbl_il_gap_Wk_ik			= 75;
	const stateType mbl_il_I_mbl_il_gap_Wk_Wi			= 76;
	const stateType mbl_il_pb_I_mbl_il_pb_K			= 77;
	const stateType mbl_il_gap_Wi_mbl_il_K_ik			= 78;
	const stateType mbl_il_gap_Wi_mbl_il_K_Wk			= 79;
	const stateType mbl_il_gap_Wi_mbl_il_gap_k_i		= 80;
	const stateType mbl_il_gap_Wi_mbl_il_gap_k_Wk	= 81;
	const stateType mbl_il_gap_iWi_mbl_il_K_k			= 82;
	const stateType mbl_il_gap_iWi_mbl_il_K_Wk		= 83;
	const stateType mbl_il_gap_i_mbl_il_K_k			= 84;
	const stateType mbl_il_gap_i_mbl_il_K_WiWk		= 85;
	const stateType mbl_il_gap_i_mbl_il_gap_Wk_k		= 86;
	const stateType mbl_il_gap_i_mbl_il_gap_Wk_Wi	= 87;
	const stateType seed_bp                         = 88;
	const stateType seed_bpI                        = 89;
	const stateType seed_bpK                        = 90;
	const stateType seed_bWi_I_bWk_K 				= 91;
	const stateType seed_bWi_I_bWk_gap_k 			= 92;
	const stateType seed_bWi_I_bWk_gap_kWk 			= 93;
	const stateType seed_bWi_I_bWk_gap_Wk 			= 94;
	const stateType seed_bWi_gap_iWi_bWk_K 			= 95;
	const stateType seed_bWi_gap_i_bWk_K 			= 96;
	const stateType seed_bWi_gap_i_bWk_gap_Wk 		= 97;
	const stateType seed_bWi_gap_Wi_bWk_K 			= 98;
	const stateType seed_bWi_gap_Wi_bWk_gap_k 		= 99;
	const stateType seed_mblIK 						=100;
	const stateType seed_mbl_bWi_I_mbl_bWk_K 		=101;
	const stateType seed_mbl_bWi_I_mbl_bWk_gap_Wk 	=102;
	const stateType seed_mbl_bWi_gap_Wi_mbl_bWk_K 	=103;


#endif /* FOLDALIGNH */
