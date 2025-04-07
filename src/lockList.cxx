#ifndef LOCKLIST
#define LOCKLIST

#include "foldalign.hxx"
#include "exception.cxx"
#include "helper.cxx"

#include <pthread.h>
#include <iostream>
#include <sstream>

class lockList {

public:
	inline lockList(const int size);
	
	inline void lock(positionType i, positionType k);
	inline void unlock(positionType i, positionType k);
	
	inline ~lockList();

private:

	const int listSize;
			
	pthread_mutex_t mem_lock;

	struct lockpos {
		positionType i;
		positionType k;
		pthread_mutex_t mutex;
		pthread_cond_t cond;
		int n_wait;
	};

	lockpos** locks;

	int* active;
	int n_active;

	inline bool useExistingLockIfExists(positionType i, positionType k);
	inline void makeNewLock(positionType i, positionType k);
	inline std::string buildUnlockErrorMessage(positionType i, positionType k);
};

inline lockList::lockList(const int size)
 : listSize(size) {
	 pthread_mutex_init(&mem_lock, NULL);

	locks = new lockpos*[size];
	active = new int[size];
	n_active = 0;
	
	for(int c = 0; c < size; c++) {
		active[c] = c;
		locks[c] = new lockpos;
		pthread_mutex_init(&locks[c]->mutex, NULL);
		pthread_cond_init(&locks[c]->cond, NULL);
		locks[c]->n_wait = 0;
	}

}

inline lockList::~lockList() {


	pthread_mutex_destroy(&mem_lock);

	delete[] active;

	for(int c=0; c < listSize; c++) {
		pthread_mutex_destroy(&locks[c]->mutex);
		pthread_cond_destroy(&locks[c]->cond);
		delete locks[c];
	}
	delete[] locks;

}

inline void lockList::lock(positionType i, positionType k) {

	pthread_mutex_lock(&mem_lock);

	if (! useExistingLockIfExists(i,k)) {
		makeNewLock(i,k);
	}
	pthread_mutex_unlock(&mem_lock);
}

inline bool lockList::useExistingLockIfExists(positionType i, positionType k) {

	bool found = false;
	for(int cc = 0; cc < n_active; cc++) {
		const int c = active[cc];
		if (i == locks[c]->i && k == locks[c]->k) {
			locks[c]->n_wait++;
			
			pthread_mutex_lock(&locks[c]->mutex);
			pthread_mutex_unlock(&mem_lock);
			pthread_cond_wait(&locks[c]->cond, &locks[c]->mutex);
			pthread_mutex_unlock(&locks[c]->mutex);
			pthread_mutex_lock(&mem_lock);

			found = true;
			break;
		}
	}

	return found;
}

inline void lockList::makeNewLock(positionType i, positionType k) {

	const int c = active[n_active];
	n_active++;
	if (n_active > listSize) {
		std::cerr << "n_active overflow in lockList. Size: " << listSize << " n_active: " << n_active << std::endl;
	}
		
	locks[c]->i = i;
	locks[c]->k = k;
	locks[c]->n_wait = 1;

}


inline void lockList::unlock(positionType i, positionType k) {

	pthread_mutex_lock(&mem_lock);

	for(int cc=0; cc < n_active; cc++) {
		const int c = active[cc];
		if (locks[c]->i == i && locks[c]->k == k) {
			locks[c]->n_wait--;
			if (locks[c]->n_wait > 0) {
				
				pthread_mutex_lock(&locks[c]->mutex);
				pthread_cond_signal(&locks[c]->cond);
				pthread_mutex_unlock(&locks[c]->mutex);
			}
			else {
				if (cc < n_active -1) {
					helper::swap(active[cc], active[n_active-1]);
				}
				n_active--;
			}

			pthread_mutex_unlock(&mem_lock);
			return;
		}
	}


	std::string error = buildUnlockErrorMessage(i, k);
	
	pthread_mutex_unlock(&mem_lock);
	throw exception(error, true);

}

inline std::string lockList::buildUnlockErrorMessage(positionType i,
							positionType k) {

	std::stringstream error;
	error << "Program error! lockList unlock. i=" << i << " k=" << k << "\n";
	for(int c = 0; c < listSize; c++) {
		error << "lock " << c << " " << locks[c]->i << " ";
		error << locks[c]->k << " " << locks[c]->n_wait << "\n";
	}

	return error.str();
}
#endif /* LOCKLIST */
