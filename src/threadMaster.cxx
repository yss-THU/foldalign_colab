#ifndef THREADMASTER
#define THREADMASTER


#include "foldalign.hxx"
#include "thread.cxx"
#include "runKargs.cxx"
#include "foldThreadHandler.cxx"
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

#include <pthread.h>
#include <iostream>

#define THREADMASTER_TEMPLATE class stmcell, class ltmcell, class startCell, \
                      bool global, bool globalPrune, bool realigning, bool mblrealign, \
											bool multiThreaded


#define THREADMASTER_TEMPLATE_PARAMETERS stmcell, ltmcell, startCell, \
                                 global, globalPrune, realigning, mblrealign, multiThreaded



template< THREADMASTER_TEMPLATE >
class threadMaster {
public:
	threadMaster(const int num_threads,
				 const positionType begin_I, const positionType end_I,
				 const positionType begin_K, const positionType end_K,
				 results& resu,
				 const sequence* const seq_1, const sequence* const seq_2,
				 arguments& arg, 
				 scorematrix& s_matrix,
				 longTermMemory< ltmcell >* const ltme,
				 shortTermMemory< stmcell >* const stme,
				 output*& out,
				 const bool lastRun,
				 constraints* cons,
				 longTermMemory< startCell >* const startCoordinates,
				 stack_ssl< mbllist, int >* const mblMem);
	
	void startNewThread(positionType i, positionType k_start, positionType k_stop,
			positionType min_top_I, positionType c_i,
			constraints*& cons, long& n_pruns, long& n_cons);

	~threadMaster();

private:
   //Struct to hold arguments passed while using pthread_create
	struct poll_args_t {
		threadMaster *t; //pointer to be passed as 'this'
		int tid; //thread_id
	};

	const int n_threads;
	results& res;
	longTermMemory< ltmcell >* const ltm;
	shortTermMemory< stmcell >* const stm;
	positionType beginI;
  // This counter is increased every time a thread starts and
  // decreased when it finishes
  int active_threads_counter;

  // Locks std::cout (Mainly for printing LS lines)
  pthread_mutex_t out_mutex;
  pthread_mutex_t active_mutex;

  // For locking during STM transfer
  pthread_mutex_t transfer_mutex;
  pthread_cond_t transfer_cond;
  bool transfer;
  pthread_mutex_t start_transfer_mutex;
  pthread_cond_t start_transfer_cond;
	
	// For thread_poll mechanism
	pthread_mutex_t thread_poll_mutex;
	pthread_cond_t *thread_poll_cond;
	bool *thread_poll_job;
	pthread_cond_t thread_poll_finish_cond;
	bool thread_poll_finish; //Wish i could be atomic
	int thread_poll_finish_counter;
	struct poll_args_t *thread_poll_args;

	positionType* k_position;
	pthread_mutex_t* thread_mutex;
	pthread_cond_t* thread_cond;
	bool* thread_locked;
	foldK<THREADMASTER_TEMPLATE_PARAMETERS>** foldks;
	results* r;
	bool* active_threads;
	bool* running;
	foldThreadHandler** threadHandlers;
	pthread_t* threads;
	pthread_attr_t thread_attr;
	int curr_thread;
	runKargs* WiWkCoreArg;
	threadClassFunc< foldK<THREADMASTER_TEMPLATE_PARAMETERS>, runKargs >* parameters;
	
	int last;

	void raiseThreadsFinish();
	void waitThreadsFinish();
	void waitThreadCheckResults();
	void readyTransfer();
	void doTransfer();
	void startNextThread();
	
	int set_affinity(int tid);
    static void* workerThread_wrapper(void *arg);
    void workerThread(int tid);
	 void threadSync();
/*  inline void lock(pthread_mutex_t& mutex, pthread_mutex_t& outmutex, std::string name = "unknown") {
    //    pthread_mutex_lock(&&outmutex);
//        std::cout <<  "Main locking " << name << std::endl;
    //    pthread_mutex_unlock(&&outmutex);

    pthread_mutex_lock(&mutex);

    //    pthread_mutex_lock(&&outmutex);
//        std::cout <<  "Main got " << name << " lock" << std::endl;
    //    pthread_mutex_unlock(&&outmutex);
  }


  inline void unlock(pthread_mutex_t& mutex, pthread_mutex_t& outmutex, std::string name = "unknown") {
    pthread_mutex_unlock(&mutex);
    //    pthread_mutex_lock(&&outmutex);
//        std::cout << "Main unlocking " << name << std::endl;
    //    pthread_mutex_unlock(&&outmutex);
  }*/
};

template< THREADMASTER_TEMPLATE >
inline threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::threadMaster(
				 const int num_threads,
				 const positionType begin_I, const positionType end_I,
				 const positionType begin_K, const positionType end_K,
				 results& resu,
				 const sequence* const seq_1, const sequence* const seq_2,
				 arguments& arg, 
				 scorematrix& s_matrix,
				 longTermMemory< ltmcell >* const ltme,
				 shortTermMemory< stmcell >* const stme,
				 output*& out,
				 const bool lastRun,
				 constraints* cons,
				 longTermMemory< startCell >* const startCoordinates,
				 stack_ssl< mbllist, int >* const mblMem)
: n_threads(num_threads),
  res(resu),
  ltm(ltme),
  stm(stme),
  beginI(begin_I),
  active_threads_counter(1),
  transfer(false),
  curr_thread(0),
  last(0) {
  
	pthread_mutex_init(&out_mutex, NULL);
	pthread_mutex_init(&active_mutex, NULL);
	pthread_mutex_init(&transfer_mutex, NULL);
	pthread_cond_init(&transfer_cond, NULL);
	pthread_mutex_init(&start_transfer_mutex, NULL);
	pthread_cond_init(&start_transfer_cond, NULL);
	pthread_attr_init(&thread_attr);
	pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

	k_position = new positionType[n_threads];
	thread_mutex = new pthread_mutex_t[n_threads];
	thread_cond = new pthread_cond_t[n_threads];
	thread_locked = new bool[n_threads];
	foldks = new foldK<THREADMASTER_TEMPLATE_PARAMETERS>*[n_threads];
	r = new results[n_threads];
	active_threads = new bool[n_threads];
	running = new bool[n_threads];
	threadHandlers = new foldThreadHandler*[n_threads];
	threads = new pthread_t[n_threads];
	WiWkCoreArg = new runKargs[n_threads];
	parameters =
		new threadClassFunc< foldK<THREADMASTER_TEMPLATE_PARAMETERS>, runKargs >[n_threads];

	thread_poll_cond = new pthread_cond_t[n_threads];
	thread_poll_job = new bool[n_threads];
	pthread_mutex_init(&thread_poll_mutex, NULL);
	pthread_cond_init(&thread_poll_finish_cond, NULL);
	thread_poll_finish = false;
	thread_poll_finish_counter = 0;
	thread_poll_args = new poll_args_t[n_threads];

	for(int c=0; c < n_threads; c++) {
		k_position[c] = end_K;
		pthread_mutex_init(&thread_mutex[c], NULL);
		pthread_cond_init(&thread_cond[c], NULL);
		thread_locked[c] = false;
		active_threads[c] = false;
		running[c] = false;

		pthread_cond_init(&thread_poll_cond[c], NULL);
		thread_poll_job[c] = false;
		thread_poll_args[c].t = this;
		thread_poll_args[c].tid = c;
	}
	for(int c=0; c < n_threads; c++) {
		int pre = c > 0 ? c-1 : n_threads-1;
		int post = c < n_threads-1 ? c+1 : 0;
		threadHandlers[c] = new foldThreadHandler(c, n_threads,
			active_threads, running, out_mutex, active_mutex,
			transfer_mutex, transfer_cond, transfer,
			start_transfer_mutex, start_transfer_cond,
			k_position[pre], k_position[c], k_position[post],
			thread_mutex[c], thread_mutex[post],
			thread_cond[c], thread_cond[post],
			thread_locked[c], thread_locked[post]);
   		
		foldks[c] = new foldK<THREADMASTER_TEMPLATE_PARAMETERS>(
			begin_I, end_I, begin_K, end_K, r[c], *seq_1, *seq_2,
			arg, s_matrix, stm, ltm, out, k_position[c], threadHandlers[c],
			lastRun, cons, startCoordinates, mblMem);
  }


	// As a dummy set the zero thread to the done state. i.e.
	// i.e. A k position which is lower enough that the next thread can run to
	// the end. The mutex and conditions are set to unlocked.
	k_position[curr_thread] = -2;
	
	// Get read to work with the next thread
	if (multiThreaded)
		curr_thread++;
//std::cerr << "NTHREADS: " << n_threads << std::endl;
	for(int c=0; c < n_threads; c++) {
		int error = pthread_create(&threads[c], &thread_attr, &workerThread_wrapper, (void*) &thread_poll_args[c]);
		if (error != 0) {
			std::string error_message = "Can not open new thread. Threading error: " + error;
			throw exception(error_message, true);
		}
	}
}

template< THREADMASTER_TEMPLATE >
int threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::set_affinity(int tid) {
	cpu_set_t mask;
	CPU_ZERO(&mask);
	CPU_SET(tid, &mask);
	return pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &mask);
}

//Static function to start the worker thread
template< THREADMASTER_TEMPLATE >
void* threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::workerThread_wrapper(void *arg) {
	poll_args_t *a = (poll_args_t*) arg;
	a->t->set_affinity(a->tid);
	a->t->workerThread(a->tid);
	a->t->threadSync();
	//std::cout << "Finishing thread: " << a->tid << "\n";
	return NULL;
}

/*
 * The main worker thread
 * Keep running jobs with argument &parameters[tid] until thread_poll_finish is true
 */
template< THREADMASTER_TEMPLATE >
void threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::workerThread(int tid) {
	//std::cout << "Creating thread " << tid << std::endl;
	pthread_mutex_lock(&thread_poll_mutex);

	do
	{
		while (thread_poll_job[tid] == false && thread_poll_finish == false)
		{
			//std::cout << "Thread " << tid << " waiting for a job\n";
			pthread_cond_wait(&thread_poll_cond[tid], &thread_poll_mutex);
		}

		if (thread_poll_job[tid] == false && thread_poll_finish == true)
			break;
		//std::cout << "Thread " << tid << " starting a job\n";
		pthread_mutex_unlock(&thread_poll_mutex);

		startThread<foldK< THREADMASTER_TEMPLATE_PARAMETERS >, runKargs>((void*) &parameters[tid]);

		pthread_mutex_lock(&thread_poll_mutex);
		thread_poll_job[tid] = false;
		//std::cout << "Thread " << tid << " job done!\n";
		pthread_cond_signal(&thread_poll_cond[tid]);
	} while (thread_poll_finish == false);

	pthread_mutex_unlock(&thread_poll_mutex);
	return;
}

/*
 * Threads share data among them. Undefined behavior might occur if one finish before the other.
 */
template< THREADMASTER_TEMPLATE >
void threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::threadSync() {
	pthread_mutex_lock(&thread_poll_mutex);
	++thread_poll_finish_counter;
	if (thread_poll_finish_counter == n_threads)
		pthread_cond_broadcast(&thread_poll_finish_cond);
	else
		while (thread_poll_finish_counter < n_threads)
			pthread_cond_wait(&thread_poll_finish_cond, &thread_poll_mutex);
	pthread_mutex_unlock(&thread_poll_mutex);
}

template< THREADMASTER_TEMPLATE >
inline void threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::startNewThread(
			positionType i, positionType k_start, positionType k_stop,
			positionType min_top_I, positionType c_i,
			constraints*& cons, long& n_pruns, long& n_cons) {
//std::cerr << "CURRENT_THREAD: " << curr_thread << " " << i << std::endl;
	WiWkCoreArg[curr_thread].i = i;
	WiWkCoreArg[curr_thread].k_start = k_start;
	WiWkCoreArg[curr_thread].k_stop = k_stop;
	WiWkCoreArg[curr_thread].min_top_I = min_top_I;
	WiWkCoreArg[curr_thread].c_i = c_i;
	WiWkCoreArg[curr_thread].cons = cons;
	WiWkCoreArg[curr_thread].n_pruns = n_pruns;
	WiWkCoreArg[curr_thread].n_cons = n_cons;

	if (multiThreaded) {
		parameters[curr_thread].tclass = foldks[curr_thread];
		parameters[curr_thread].func = &foldK< THREADMASTER_TEMPLATE_PARAMETERS >::runK;
		parameters[curr_thread].para = &WiWkCoreArg[curr_thread];
				
		waitThreadCheckResults();

	 	pthread_mutex_lock(&active_mutex);
		pthread_mutex_lock(&start_transfer_mutex);

		readyTransfer();

		pthread_mutex_lock(&transfer_mutex);

		doTransfer();
	
		// Start the next thread

		if (i == beginI) {
			last = curr_thread;
		}

		k_position[curr_thread] = k_stop;
		startNextThread();
//std::cout << "Curr_thread: i-" << i << " Score:" << res.getScore() << " r[]Score:" << r[curr_thread].getScore() << std::endl;

		pthread_mutex_unlock(&transfer_mutex);
		pthread_mutex_unlock(&start_transfer_mutex);
		pthread_mutex_unlock(&active_mutex);
	
		// Update the current thread and lead thread index
		curr_thread++;
		if (curr_thread >= n_threads) {
			curr_thread = 0;
		}
	}
	else {
		doTransfer();
		void* coreArg = &WiWkCoreArg[curr_thread];
		(foldks[curr_thread]->runK)(coreArg);
//std::cout << "Curr_thread: i-" << i << " Score:" << res.getScore() << " r[]Score:" << r[curr_thread].getScore() << std::endl;

		last = curr_thread;

		if ((r[curr_thread].getLogScore() > res.getLogScore()) ||
			 ((r[curr_thread].getLogScore() == res.getLogScore()) && (r[curr_thread].getScore() > res.getScore()))) {
			res = r[curr_thread];
		}
	}
}

template< THREADMASTER_TEMPLATE >
inline void threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::startNextThread() {
			
	running[curr_thread] = true;
	active_threads[curr_thread] = false;

	//&parameters[curr_thread] must be set, before calling this
	pthread_mutex_lock(&thread_poll_mutex);
	if (thread_poll_job[curr_thread] == true)
	{
		std::cerr << "Bug at waitThreadCheckResultes(). startNextThread() have started on working thread. exiting...\n";
		exit(1);
	}
	//std::cout << "New job for thread: " << curr_thread << "\n";
	thread_poll_job[curr_thread] = true;
	pthread_cond_signal(&thread_poll_cond[curr_thread]);
	pthread_mutex_unlock(&thread_poll_mutex);
}

template< THREADMASTER_TEMPLATE >
inline void threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::readyTransfer() {
	// Lock down alle threads before the transfers

	transfer = true;

	// Wait for threads to stop. The last thread
	// sends a signal
	for(int c = 0; c < n_threads; c++) {
		if (active_threads[c] == true) {
			pthread_mutex_unlock(&active_mutex);

			pthread_cond_wait(&start_transfer_cond, &start_transfer_mutex);
			pthread_mutex_unlock(&start_transfer_mutex);
			pthread_mutex_lock(&active_mutex);
			pthread_mutex_lock(&start_transfer_mutex);
			break;
		}
	}
}

template< THREADMASTER_TEMPLATE >
inline void threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::doTransfer() {

	// Slide the stm matrix along the I-sequence.
	if (active_threads_counter >= n_threads) {
		stm->transfer();
	}
	else {active_threads_counter++;}

	if (ltm != 0) {ltm->transfer();}

	// Restart all the threads
	transfer = false;
	pthread_cond_broadcast(&transfer_cond);

}				

template< THREADMASTER_TEMPLATE >
inline void threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::waitThreadCheckResults() {

	if (active_threads_counter >= n_threads) {
		int lead_thread = curr_thread +1;
		if (lead_thread == n_threads) {lead_thread = 0;}

		// Wait for the next thread to finish
		pthread_mutex_lock(&thread_poll_mutex);
		if (thread_poll_job[lead_thread] == true)
			pthread_cond_wait(&thread_poll_cond[lead_thread], &thread_poll_mutex);
		pthread_mutex_unlock(&thread_poll_mutex);
				
		if ( ! (global || mblrealign || realigning ) ) {
			if ((r[lead_thread].getLogScore() > res.getLogScore()) ||
				 ((r[lead_thread].getLogScore() == res.getLogScore()) && (r[lead_thread].getScore() > res.getScore()))) {
				res = r[lead_thread];
			}
		}
	}
}

template< THREADMASTER_TEMPLATE >
void threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::raiseThreadsFinish() {
	//std::cout << "Telling to all threads to finish...\n";

	pthread_mutex_lock(&thread_poll_mutex);
	thread_poll_finish = true;
	for(int c=0; c < n_threads; c++) {
		pthread_cond_signal(&thread_poll_cond[c]);
	}
	pthread_mutex_unlock(&thread_poll_mutex);
}

template< THREADMASTER_TEMPLATE >
inline void threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::waitThreadsFinish() {

	// Make sure all threads are done and clean up.
	//std::cout << "Waiting for threads to finish...\n";
	for(int c=0; c < n_threads; c++) {
		raiseThreadsFinish();
		//Verify that all threads exited
		pthread_join(threads[c], NULL);

		if ( ! (global || mblrealign || realigning ) ) {
			if ((r[c].getLogScore() > res.getLogScore()) ||
				 ((r[c].getLogScore() == res.getLogScore()) && (r[c].getScore() > res.getScore()))) {
				res = r[c];
			}
		}
	}

	if (global || mblrealign || realigning ) {
		res = r[last];
	}

}

template< THREADMASTER_TEMPLATE >
inline threadMaster< THREADMASTER_TEMPLATE_PARAMETERS >::~threadMaster() {

	waitThreadsFinish();

	pthread_attr_destroy(&thread_attr);
	pthread_mutex_destroy(&out_mutex);
	pthread_mutex_destroy(&active_mutex);
	pthread_mutex_destroy(&transfer_mutex);
	pthread_cond_destroy(&transfer_cond);
	pthread_mutex_destroy(&start_transfer_mutex);
	pthread_cond_destroy(&start_transfer_cond);

	for(int c=0; c < n_threads; c++) {
		pthread_mutex_destroy(&thread_mutex[c]);
		pthread_cond_destroy(&thread_cond[c]);
		delete foldks[c];
		delete threadHandlers[c];
	}
	delete[] threads;
	delete[] threadHandlers;
	delete[] thread_mutex;
	delete[] thread_cond;

	delete[] k_position;
	delete[] thread_locked;
	delete[] WiWkCoreArg;
	delete[] parameters;

	delete[] foldks;
	delete[] r;
	delete[] active_threads;
	delete[] running;

	delete[] thread_poll_cond;
	delete[] thread_poll_job;
	delete[] thread_poll_args;
}
#endif
