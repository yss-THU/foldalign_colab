#ifndef OUTPUTH
#define OUTPUTH

// These are character used to indicate base-pairs
// There are two sets. The first one is used when the base-pair is conserved
// between the sequences, the second is used for insert base-pairs.


// Normal base pairs
static const char open_bp  = '(';
static const char close_bp = ')';

// Insert base pairs
static const char open_one_seq  = '<';
static const char close_one_seq = '>';

// This is the unpaired character
static const char unpaired = '.';
	
//----------------------------------------------------
// Column format constants

// The length of front text field
static const int front_length = 22;
	
// Length of the short names
static const int length_name_field = 13;
	
// Number of nucleotides pr line in summary alignment
static const int line_length = 40;
	
// The width of a column in the column format.
static const int space = 10;


//----------------------------------------------------
// Stockholm format constants

static const int stk_feat_size = 18;

static const int stk_name_size = 18;






#endif /* OUTPUTH */
