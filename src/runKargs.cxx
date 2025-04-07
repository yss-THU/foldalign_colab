#ifndef RUNKARGS
#define RUNKARGS

#include "constraints.cxx"
#include "foldalign.hxx"
#include <pthread.h>

class runKargs {
public:
	positionType i;
	positionType k_start;
	positionType k_stop;
	positionType min_top_I;
	positionType c_i;
	constraints* cons;
	long n_pruns;
	long n_cons;
};




#endif
