#ifndef FOLDTHREADHANDLER
#define FOLDTHREADHANDLER

#include "foldalign.hxx"

#include <pthread.h>
#include <iostream>
class foldThreadHandler {
public:

	foldThreadHandler(
		const int thread_num,
		const int num_threads,
		bool*& active_thread,
		bool*& running_thread,
		pthread_mutex_t& output_mut,
		pthread_mutex_t& active_mut,
		pthread_mutex_t& transfer_mut,
		pthread_cond_t& transfer_condi,
		bool& trans_now,
		pthread_mutex_t& start_mutex,
		pthread_cond_t& start_cond,
		positionType& pre_k_pos,
		positionType& k_pos,
		positionType& post_k_pos,
		pthread_mutex_t& c_thread_mutex,
		pthread_mutex_t& po_thread_mutex,
		pthread_cond_t& c_thread_cond,
		pthread_cond_t& po_thread_cond,
		bool& c_thread_locked,
		bool& po_thread_locked
	) : thread_number(thread_num),
			n_threads(num_threads),
	    post_thread(thread_number +1 < n_threads ? thread_number+1 : 0),
			active_threads(active_thread),
			running_threads(running_thread),
			output_mutex(output_mut),
			active_mutex(active_mut),
			transfer_mutex(transfer_mut),
			transfer_cond(transfer_condi),
			transfer_now(trans_now),
			start_transfer_mutex(start_mutex),
			start_transfer_cond(start_cond),
			pre_position(pre_k_pos),
			position(k_pos),
			post_position(post_k_pos),
			curr_thread_mutex(c_thread_mutex),
			post_thread_mutex(po_thread_mutex),
			curr_thread_cond(c_thread_cond),
			post_thread_cond(po_thread_cond),
			thread_locked(c_thread_locked),
			post_thread_locked(po_thread_locked)
	{
	};

	void releasePostThread();
	void waitPreThread();
	void waitTransfer();
	void finishK();

	int getThreadNumber() const {return thread_number;};
	
	void lockAll();
	
	void lockCurrent() {
		pthread_mutex_lock(&curr_thread_mutex);
	};
	void lockPost() {
		pthread_mutex_lock(&post_thread_mutex);
	};
	void lockOutput() {
		pthread_mutex_lock(&output_mutex);
	};
	void unlockCurrent() {
		pthread_mutex_unlock(&curr_thread_mutex);
	};
	void unlockPost() {
		pthread_mutex_unlock(&post_thread_mutex);
	};
	void unlockOutput() {
		pthread_mutex_unlock(&output_mutex);
	};

private:
	foldThreadHandler();
	
	void unlockAll();
	inline void releasePost();
	inline void waitPre();
	
	const int thread_number;
	const int n_threads;
	const int post_thread;
	bool*& active_threads;
	bool*& running_threads;
	pthread_mutex_t& output_mutex;
	pthread_mutex_t& active_mutex;
	pthread_mutex_t& transfer_mutex;
	pthread_cond_t& transfer_cond;
	bool& transfer_now;
	pthread_mutex_t& start_transfer_mutex;
	pthread_cond_t& start_transfer_cond;
	positionType& pre_position;
	positionType& position;
	positionType& post_position;
	pthread_mutex_t& curr_thread_mutex;
	pthread_mutex_t& post_thread_mutex;
	pthread_cond_t& curr_thread_cond;
	pthread_cond_t& post_thread_cond;
	bool& thread_locked;
	bool& post_thread_locked;

	inline bool anyActive() {
		for(int c = 0; c < n_threads; c++) {
			if (active_threads[c] == true) {
				return true;
			}
		}
		return false;
	};

/*	
  inline void lock(pthread_mutex_t& mutex, std::string name = "unknown") {
//            pthread_mutex_lock(&&outmutex);
//	          std::cout << thread_number << " locking " << name << " " << &mutex << std::endl;
//            pthread_mutex_unlock(&&outmutex);
    pthread_mutex_lock(&mutex);
//            pthread_mutex_lock(&&outmutex);
//    std::cout << " " << thread_number << " got " << name << " lock" << std::endl;
//            pthread_mutex_unlock(&&outmutex);
  };


  inline void unlock(pthread_mutex_t& mutex, std::string name = "unknown") {
//            pthread_mutex_lock(&&outmutex);
//    std::cout << thread_number << " unlocking " << name << " " << &mutex << std::endl;
//            pthread_mutex_unlock(&&outmutex);
    pthread_mutex_unlock(&mutex);
//           pthread_mutex_lock(&&outmutex);
//    std::cout << " " << thread_number << " unlocked " << name << std::endl;
//            pthread_mutex_unlock(&&outmutex);
  };
*/	
};

inline void foldThreadHandler::lockAll() {
	pthread_mutex_lock(&active_mutex);
	pthread_mutex_lock(&post_thread_mutex);
	pthread_mutex_lock(&curr_thread_mutex);
	pthread_mutex_lock(&start_transfer_mutex);
	pthread_mutex_lock(&transfer_mutex);
}

inline void foldThreadHandler::unlockAll() {
	pthread_mutex_unlock(&active_mutex);
	pthread_mutex_unlock(&post_thread_mutex);
	pthread_mutex_unlock(&curr_thread_mutex);
	pthread_mutex_unlock(&start_transfer_mutex);
	pthread_mutex_unlock(&transfer_mutex);
}

inline void foldThreadHandler::releasePostThread() {

	active_threads[thread_number] = true;

	while (transfer_now) {
		active_threads[thread_number] = false;

		if (!anyActive()) {

			pthread_cond_signal(&start_transfer_cond);
		}
     
		pthread_mutex_unlock(&post_thread_mutex);
		pthread_mutex_unlock(&curr_thread_mutex);
		pthread_mutex_unlock(&start_transfer_mutex);
		pthread_mutex_unlock(&active_mutex);

		pthread_cond_wait(&transfer_cond, &transfer_mutex);
		pthread_mutex_unlock(&transfer_mutex);

		lockAll();
		
		active_threads[thread_number] = true;
	}

	releasePost();

	if (position <= pre_position+1) {
		active_threads[thread_number] = false;
		waitPre();
	}
	else {
		unlockAll();
	}
}

inline void foldThreadHandler::releasePost() {

	if (post_thread_locked && post_position > position+1) {
		post_thread_locked = false;
		if (running_threads[post_thread]) {
			active_threads[post_thread] = true;
		}
		pthread_cond_signal(&post_thread_cond);
	}

}

inline void foldThreadHandler::waitPre() {

	thread_locked = true;

	pthread_mutex_unlock(&post_thread_mutex);
	pthread_mutex_unlock(&start_transfer_mutex);
	pthread_mutex_unlock(&transfer_mutex);
	pthread_mutex_unlock(&active_mutex);

	pthread_cond_wait(&curr_thread_cond, &curr_thread_mutex);

	pthread_mutex_unlock(&curr_thread_mutex);
//	lockAll();

}

inline void foldThreadHandler::finishK() {

	active_threads[thread_number] = false;
	running_threads[thread_number] = false;

	releasePost();

	if (transfer_now && !anyActive()) {
		pthread_cond_signal(&start_transfer_cond);
	}

	unlockAll();
}

/*

inline void foldThreadHandler::releasePostThread() {
////	lock(post_thread_mutex, "post");
	pthread_mutex_lock(&post_thread_mutex);
	// Release the thread with position i-1 (post)
	if (post_thread_locked && post_position > position+1) {
std::cout << thread_number << " Post pos " << post_position << " > " << position << "+1" << std::endl;
		post_thread_locked = false;
////		lock(start_transfer_mutex, "start_transfer");
		pthread_mutex_lock(&start_transfer_mutex);
		const int post_thread = thread_number +1 < n_threads ? thread_number+1 : 0;
		if (running_threads[post_thread]) {
////			lock(active_mutex, "anyActive");
			pthread_mutex_lock(&active_mutex);
			active_threads[post_thread] = true;
//			unlock(active_mutex, "anyActive");
			pthread_mutex_unlock(&active_mutex);
		}
		pthread_cond_signal(&post_thread_cond);
//		unlock(start_transfer_mutex, "start_transfer");
		pthread_mutex_unlock(&start_transfer_mutex);
std::cout << thread_number << " Post pos released" << std::endl;
	}
//	unlock(post_thread_mutex, "post");
	pthread_mutex_unlock(&post_thread_mutex);
}


inline void foldThreadHandler::waitPreThread() {
	// Wait for the thread at i +1 to move ahead on the K sequence
////	lock(curr_thread_mutex, "curr");
	pthread_mutex_lock(&curr_thread_mutex);
	if (position <= pre_position+1) {
std::cout << thread_number << " wait pre " << position << " <= " << pre_position << "+1 TN:" << transfer_now << std::endl;
////		lock(start_transfer_mutex, "start_transfer");
		pthread_mutex_lock(&start_transfer_mutex);
////		lock(active_mutex, "anyActive");
		pthread_mutex_lock(&active_mutex);
		active_threads[thread_number] = false;
//		unlock(active_mutex, "anyActive");
		pthread_mutex_unlock(&active_mutex);

		if (transfer_now && !anyActive()) {
//			unlock(curr_thread_mutex, "curr");
			pthread_mutex_unlock(&curr_thread_mutex);
////			lock(transfer_mutex, "transfer");
			pthread_mutex_lock(&transfer_mutex);
			pthread_cond_signal(&start_transfer_cond);
//			unlock(start_transfer_mutex, "start_transfer");
			pthread_mutex_unlock(&start_transfer_mutex);

			pthread_cond_wait(&transfer_cond, &transfer_mutex);
//			unlock(transfer_mutex, "transfer");
			pthread_mutex_unlock(&transfer_mutex);
////			lock(curr_thread_mutex, "curr");
			pthread_mutex_lock(&curr_thread_mutex);
			if (position <= pre_position+1) {
				thread_locked = true;
				pthread_cond_wait(&curr_thread_cond, &curr_thread_mutex);
			}
		}
		else {
//			unlock(start_transfer_mutex, "start_transfer");
			pthread_mutex_unlock(&start_transfer_mutex);
			thread_locked = true;
			pthread_cond_wait(&curr_thread_cond, &curr_thread_mutex);
		}
std::cout << thread_number << " wait for pre done " << position << std::endl;

	}
//	unlock(curr_thread_mutex, "curr");
	pthread_mutex_unlock(&curr_thread_mutex);
}



/inline void foldThreadHandler::waitPreThread() {
	// Wait for the thread at i +1 to move ahead on the K sequence
//	lock(curr_thread_mutex, "curr");
pthread_mutex_lock(&curr_thread_mutex);
	if (position <= pre_position+1) {

//		unlock(curr_thread_mutex, "curr");
pthread_mutex_unlock(&curr_thread_mutex);
		waitTransfer();
//		lock(curr_thread_mutex, "curr");
pthread_mutex_lock(&curr_thread_mutex);

//		lock(active_mutex, "anyActive");
pthread_mutex_lock(&active_mutex);
		active_threads[thread_number] = false;
//		unlock(active_mutex, "anyActive");
pthread_mutex_unlock(&active_mutex);
		thread_locked = true;
		pthread_cond_wait(&curr_thread_cond, &curr_thread_mutex);

	}
//	unlock(curr_thread_mutex, "curr");
pthread_mutex_unlock(&curr_thread_mutex);
}


inline void foldThreadHandler::waitTransfer() {
	// Wait for the transfer if necessary
//	lock(start_transfer_mutex, "start_transfer");
	pthread_mutex_lock(&start_transfer_mutex);
	while (transfer_now) {
std::cout << thread_number << " wait for transfer" << std::endl;
//		lock(transfer_mutex, "transfer");
		pthread_mutex_lock(&transfer_mutex);
//		lock(active_mutex, "anyActive");
		pthread_mutex_lock(&active_mutex);
		active_threads[thread_number] = false;
//		unlock(active_mutex, "anyActive");
		pthread_mutex_unlock(&active_mutex);

		if (!anyActive()) {pthread_cond_signal(&start_transfer_cond);}
     
//		unlock(start_transfer_mutex, "start_transfer");
		pthread_mutex_unlock(&start_transfer_mutex);
		pthread_cond_wait(&transfer_cond, &transfer_mutex);
//		lock(active_mutex, "anyActive");
		pthread_mutex_lock(&active_mutex);
		active_threads[thread_number] = true;
//		unlock(active_mutex, "anyActive");
		pthread_mutex_unlock(&active_mutex);
//		unlock(transfer_mutex, "transfer");
		pthread_mutex_unlock(&transfer_mutex);

//		lock(start_transfer_mutex, "start_transfer");
		pthread_mutex_lock(&start_transfer_mutex);
std::cout << thread_number << " transfer done" << std::endl;
	}
//	unlock(start_transfer_mutex, "start_transfer");
	pthread_mutex_unlock(&start_transfer_mutex);
}


inline void foldThreadHandler::finishK() {
//  lock(start_transfer_mutex, "start_transfer restart end");
	pthread_mutex_lock(&start_transfer_mutex);
//  lock(transfer_mutex, "transfer restart end");
	pthread_mutex_lock(&transfer_mutex);

//	lock(active_mutex, "anyActive");
	pthread_mutex_lock(&active_mutex);
	active_threads[thread_number] = false;
////	unlock(active_mutex, "anyActive");
	pthread_mutex_unlock(&active_mutex);
	if (transfer_now && !anyActive()) {
		pthread_cond_signal(&start_transfer_cond);
	}
	running_threads[thread_number] = false;

////  unlock(transfer_mutex, "transfer end");
	pthread_mutex_unlock(&transfer_mutex);
////  unlock(start_transfer_mutex, "start transfer end");
	pthread_mutex_unlock(&start_transfer_mutex);
}
*/
#endif /* FOLDTHREADHANDLER */

