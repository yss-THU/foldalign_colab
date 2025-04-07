#ifndef THREAD
#define THREAD

#include<iostream>
// This is an attempt at making the necessary functions for starting new threads

template< class aclass, class parameters >
struct threadClassFunc {
  aclass* tclass;
  void* (aclass::*func)(void* );
  parameters* para;
};

//extern "C" 
template< class aclass, class parameters >
void* startThread(void* arg) {
	
  threadClassFunc<aclass, parameters>* cf = 
    (threadClassFunc<aclass, parameters>*) arg;

  aclass& run = *cf->tclass;

  run.runK(cf->para);

  return 0;
}

#endif /* THREAD */
