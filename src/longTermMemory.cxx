#ifndef LONGTEMMEMORY
#define LONGTEMMEMORY
#include "longCell.cxx"
#include "helper.cxx"
#include "two_link.cxx"
#include "foldalign.hxx"
#include "exception.cxx"
#include "lockList.cxx"

#include <string>

#define TEMP_DEF_LTM template< class cell >
#define TEMP_SPES_LTM cell

TEMP_DEF_LTM
class longTermMemory {
public:

	inline longTermMemory(const positionType k_dim,
			const lengthType i_dim,
			positionType i_pos,
			positionType k_pos,
			const bool& global,
			const int num_threads = 1,
			const positionType k_off = 0,
			const positionType i_fix = -1);
	
	inline cell* getNextPos(positionType i, positionType k,
				lengthType& Wi, lengthType& Wk) {
	// This function is used when i & k are known but Wi and Wk are unknown
		// and the needed cell is just the next valid cell
		(this->*p2calcPos)(i, k);
		if (matrix[i][k] == 0) {
			Wi = 0;
			Wk = 0; 
			return 0;
		}
		return matrix[i][k]->getNext(Wi,Wk);
	};

	inline cell* getPos(positionType i, positionType k,
					lengthType& Wi, lengthType& Wk, const int thread = 0) {
		// This function is used when all coordinates are known (i,k,Wi,Wk)
		// and a specific cell is needed.
		(this->*p2calcPos)(i, k);
		if (matrix[i][k] == 0) {
			return 0;
		}
		return matrix[i][k]->getPos(Wi, Wk);
	};

	inline cell* putPos(positionType i, positionType k,
					const lengthType Wi, const lengthType Wk, const int thread = 0) {
		plock(mem_lock, thread, "rc memlock");
		(this->*p2calcPos)(i, k);
		if (matrix[i][k] == 0) {
			matrix[i][k] = new two_link< cell >();
			if (matrix[i][k] == 0) {
				std::string error = "LTM error. Most likely cause: Out of memory";
				punlock(mem_lock, thread, "rc memlock");
				throw exception(error, false);
			}
		}
		punlock(mem_lock, thread, "rc memlock");
		return matrix[i][k]->putNext(Wi, Wk);
	};

	inline two_link< cell >* getLockAndResetSubMatrix(positionType i, positionType k, const int thread = 0);

	inline void deletePos(positionType i, positionType k) {
		(this->*p2calcPos)(i, k);
		if (matrix[i][k] != 0) {
			matrix[i][k]->removeRest();
		}
	};


	inline void switchWi(positionType i, positionType k) {
		(this->*p2calcPos)(i, k);
		if (matrix[i][k] != 0) {
			matrix[i][k]->lastWk();
		}
	};
	
	inline void transfer() {
		(this->*p2trans)();
	};
	
	inline void setPositions(const positionType i_pos,
				 const positionType k_pos, const int thread = 0) {
		// This function is not checked by the testLongTermMemory.cxx class
		i_position = i_pos;
		k_position = k_pos;
	};
	
	inline positionType getIStart() const {
		return i_start;
	}
	
	inline long getSize() const {
		return matSize + two_link< cell >::getSize();
	}

	inline void resetCurrent(positionType i, positionType k, const int thread =0);

	inline void unlock(positionType i, positionType k, const int thread=0);
	
	inline ~longTermMemory();

private:
	inline void plock(pthread_mutex_t& mutex, const int thread = -1, std::string name = "unknown") {
//						pthread_mutex_lock(&outmutex);
//				std::cerr << thread << " locking " << name << std::flush;
//						pthread_mutex_unlock(&outmutex);
		pthread_mutex_lock(&mutex);
//						pthread_mutex_lock(&outmutex);
//				std::cerr << thread << " got " << name << " lock" << std::endl;
//						pthread_mutex_unlock(&outmutex);
	}

	inline void punlock(pthread_mutex_t& mutex, const int thread = -1, std::string name = "unknown") {
//						pthread_mutex_lock(&outmutex);
//				std::cerr << thread << " unlocking " << name << std::flush;
//						pthread_mutex_unlock(&outmutex);
		pthread_mutex_unlock(&mutex);
//						pthread_mutex_lock(&outmutex);
//				std::cerr << thread <<" unlocked " << name << std::endl;
//						pthread_mutex_unlock(&outmutex);
	}

	// No default constructor
	inline longTermMemory();
	
	// Changes the i, k coordinates into the internal coordinates
	// Scan coordinates	
	inline void calcPos_local(positionType& i, positionType& k) const
	{
		//positionType oi = i; positionType ok = k;
		i -= i_position;
		k = k_position - k;
		//	if ( i >= i_dimension || i < 0) {std::cerr << "LTM I spot trouble " << oi << " " << ok << " " << i << " " << k << " " <<	k_dimension << " " << k_offset << std::endl;}
		//	if ( k >= k_dimension || k < 0) {std::cerr << "K spot trouble " << k << std::endl;}

	};
	
	// Changes the i, k coordinates into the internal coordinates	
	// global and realignment coordinates
	inline void calcPos_global(positionType& i, positionType& k) const
	{
		//		positionType oi = i; positionType ok = k;
		k = k - i + k_offset;
		i -= i_position;
		//		if ( i >= i_dimension || i < 0) {std::cerr << "LTM I spot trouble " << oi << " " << ok << " " << i << " " << k << " " <<	k_dimension << " " << k_offset << std::endl;}
		//		if ( k >= k_dimension || k < 0) {std::cerr << "LTM K spot trouble " << oi << " " << ok << " " << i << " " << k << " " <<	k_dimension << " " << k_offset << std::endl;}
	};

	// The dimension along the k-chunk
	const positionType k_dimension;
	
	// The dimension along the I-sequence.
	const positionType i_dimension;
	
	// The matrixs position on the I-sequence
	positionType i_position;
	
	// The matrixs position on the K-sequence
	positionType k_position;

	// The data matrix
	two_link< cell >*** matrix;
	
	// The offset between the k and i start positions in global alignment
	const positionType k_offset;

	// Global align? true or false
	const bool global_align;
	
	// The number of threads
	const int n_threads;
	
	// Fixed i start position
	const positionType i_start;
	
	// Different verions of the calcPos function is needed for local and global
	// alignment. This pointer points to the function in use.
	void (longTermMemory::* const p2calcPos)(positionType& i, positionType& k) const;

	// A transfer function is needed for local but not for global alignment.
	// The pointer keeps track of the transfer function used.
	void (longTermMemory::* const p2trans)();

	inline void transfer_local();
	inline void transfer_global() {};
	
	const long matSize;
	
	pthread_mutex_t mem_lock;
	
	lockList* locks;
};

TEMP_DEF_LTM
inline longTermMemory< TEMP_SPES_LTM >::longTermMemory(
			const positionType k_dim,
			const lengthType i_dim,
			positionType i_pos,
			positionType k_pos,
			const bool& global,
			const int num_threads,
			const positionType k_off,
			const positionType i_fix)
: k_dimension(k_dim),
	i_dimension(positionType(i_dim)),
	i_position(i_pos),
	k_position(k_pos), 
	k_offset(k_off),
	global_align(global),
	n_threads(num_threads+1),
	i_start(i_fix),
	p2calcPos( global_align ? &longTermMemory::calcPos_global : &longTermMemory::calcPos_local),
	p2trans( global_align ? &longTermMemory::transfer_global : &longTermMemory::transfer_local),
	matSize(char_size * i_dimension * k_dimension * sizeof(cell*) / mem_scale)
{

//	std::cout << "LTM START: " << i_dimension << " " << k_dimension << " " << " " << global << " " << this << std::endl;
	pthread_mutex_init(&mem_lock, NULL);

	locks = new lockList(2*n_threads);
	
	// I do not suffer fools gladly, and fools with badges never.
	// I want no interference from this hulking boob. Is that clear? 

	std::string error = "Could not allocate longTermMemory. Most likely cause: Out of memory";
	matrix = new two_link< cell >**[i_dimension];
	if (matrix == 0) {throw exception(error, false);}
	for(positionType i = 0; i < i_dimension; i++) {
		matrix[i] = new two_link< cell >*[k_dimension];
		if (matrix[i] == 0) {throw exception(error, false);}
		for( positionType k = 0; k < k_dimension; k++ ) {matrix[i][k] = 0;}
	}
//std::cerr << "LTM DIM. I: " << i_dimension << " K: " << k_dimension << std::endl;
}

TEMP_DEF_LTM
inline longTermMemory< TEMP_SPES_LTM >::~longTermMemory() {

//std::cout << "I am the LTM DESTORYER! " << this << std::endl;

	pthread_mutex_destroy(&mem_lock);

	delete locks;

	for(positionType i = 0; i < i_dimension; i++) {
		for( positionType k = 0; k < k_dimension; k++ ) {
			if (matrix[i][k] != 0) {
				delete matrix[i][k];
				matrix[i][k] = 0;
			}
		}
		delete[] matrix[i];
	}
	delete[] matrix;
	//			std::cerr << " and i am done" << std::endl;
}

TEMP_DEF_LTM
inline void longTermMemory< TEMP_SPES_LTM >::transfer_local() {
	plock(mem_lock,	-1, "ul memlock");

	// There is no transfer required during global alignment
	// Empty the oldest i array
	positionType i = i_dimension -1;
	for(positionType k = 0; k < k_dimension; k++) {
		if (matrix[i][k] != 0) {
			delete matrix[i][k];
			matrix[i][k] = 0;
		}
//		matrix[i][k]->resetCurrent();
//		matrix[i][k]->removeRest();
	}

	// Shift the matrix
	two_link< cell >** tmp;
	tmp = matrix[i];
	for( ; i > 0; i--) {
		matrix[i] = matrix[i-1];
	}
	matrix[0] = tmp;

	i_position--;

	punlock(mem_lock,	-2, "ul memlock");
}

TEMP_DEF_LTM
inline void longTermMemory< TEMP_SPES_LTM >::resetCurrent(positionType i,
								positionType k,
								const int thread) {

	locks->lock(i, k);
	plock(mem_lock, thread, "rc memlock");

	(this->*p2calcPos)(i, k);
	if (matrix[i][k] != 0) {matrix[i][k]->resetCurrent();}

	punlock(mem_lock, thread, "rc memlock");
}
	
TEMP_DEF_LTM
inline two_link< cell >* longTermMemory< TEMP_SPES_LTM >::getLockAndResetSubMatrix(positionType i,
								positionType k,
								const int thread) {

	locks->lock(i, k);
	plock(mem_lock, thread, "rc memlock");

	(this->*p2calcPos)(i,k);
	if (matrix[i][k] == 0) {
		matrix[i][k] = new two_link< cell >();
		if (matrix[i][k] == 0) {
			std::string error = "LTM error. Most likely cause: Out of memory";
			punlock(mem_lock, thread, "rc memlock");
			throw exception(error, false);
		}
	}
	matrix[i][k]->resetCurrent();
	punlock(mem_lock, thread, "rc memlock");
	return matrix[i][k];
}

 
TEMP_DEF_LTM
inline void longTermMemory< TEMP_SPES_LTM >::unlock(positionType i, 
								positionType k, const int thread) {

	locks->unlock(i, k);
}


#endif /* SHORTTEMMEMORY */
