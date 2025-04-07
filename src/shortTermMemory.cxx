#ifndef SHORTTEMMEMORY
#define SHORTTEMMEMORY

#include "arguments.cxx"
#include "helper.cxx"
#include "scorematrix.cxx"
#include "exception.cxx"
#include "foldalign.hxx"
#include "stmWiWkSubMatrix.cxx"

#include <pthread.h>
//#include <stack>

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

#define TEMP_DEF_SHORTTERMMEMORY template< class cell >
#define TEMP_SPES_SHORTTERMMEMORY cell


TEMP_DEF_SHORTTERMMEMORY
class shortTermMemory {
public:

  inline shortTermMemory(const positionType& i_dim, 
			 const positionType& k_dim,
			 const lengthType& Wi_dim,
			 arguments& argu,
			 scorematrix& sm,
			 const int num_threads = 1,
			 const positionType k_off = 0, 
			 const lengthType Wi_minimum = 0);
	
  inline cell* getPos(positionType i, positionType k,
		      lengthType Wi, lengthType Wk, const int thread = 0) {

//    pthread_mutex_lock(&memlock);
		(this->*p2calcpos)(i, k);
		return matrix[i][k]->getPos(Wi, Wk);
  };
	
  inline cell* putPos(positionType i, positionType k, 
		      lengthType Wi, lengthType Wk, const int thread = 0) {
		(this->*p2calcpos)(i, k);
		return matrix[i][k]->putPos(Wi, Wk);
	};

  inline lengthType getWimax(positionType i, positionType k,
			     const int thread = 0) {
		(this->*p2calcpos)(i, k);
		return matrix[i][k]->getWimax();
  };

  inline void transfer();
	inline void clear_old_cells(int cells_to_clean);

	inline stmWiWkSubMatrix< cell >* getSubMatrix(positionType i, positionType k) {

		// Before this function can be called, make sure that the lock
		// function has been called
		(this->*p2calcpos)(i, k);
		stmWiWkSubMatrix< cell >* returnValue = matrix[i][k];
		
		return returnValue;
	};
	
	inline positionType getIdimension() const {return i_dimension;}
	inline positionType getKdimension() const {return k_dimension;}
	inline lengthType getWk_zero() const {return Wk_zero;};
	inline lengthType getWi_minimum() const {return Wi_min;};
	inline lengthType getK_offset() const {return k_offset;};
	
  inline void setPositions(const positionType i_pos, 
			   const positionType k_pos) {
    i_position = i_pos;
    k_position = k_pos;
  };

  inline void set_I_Position(const positionType pos) {i_position = pos;};

  inline void set_K_Position(const positionType pos) {k_position = pos;};

  inline long getSize() const {
    return 0;
    //    return char_size*(cell_count * cell_size + mat3d_size +
    //		      matWk_cell_count * ptr_size)/mem_scale;
  };
	
	inline void lock() {pthread_mutex_lock(&memlock);};
	inline void unlock() {pthread_mutex_unlock(&memlock);};
	

  inline ~shortTermMemory();

private:

  const positionType i_dimension;

  // The dimension along the k-chunk
  const positionType k_dimension;
	
  // The window size along sequence I dimension
  const lengthType Wi_dimension;

  // The window size along sequence K dimension
  // This is not const since it is depended on the Wi. In the extremes of Wi
  // it is not const.
  lengthType Wk_dimension;
  lengthType Wk_zero;

  arguments& arg;
	
  const int n_threads;

  // The matrixs position on the I-sequence
  positionType i_position;
	
  positionType k_position;

  // The data matrix
  stmWiWkSubMatrix< cell >*** matrix;
	
  // The length of the sequences
  const positionType seq_length1;
  const positionType seq_length2;

  lengthType wimax_init;

  inline void calcPos_local(positionType& i, positionType& k) const {
//    positionType oi = i; positionType ok = k;
    i = i_position - i;
    k = k_position - k;
//    if (i < 0 || i >= i_dimension || k < 0 || k >= k_dimension) {
//	std::cerr << "STM trouble: " << i << " " << i_dimension << " " << oi << " " << k << " " << k_dimension << " " << ok << std::endl;
//    }
  };

  inline void calcPos_global(positionType& i, positionType& k) const {
//    positionType oi = i; positionType ok = k;
    k = k - i + k_offset;
    i = i_position - i;
//    if (i < 0 || i >= i_dimension || k < 0 || k >= k_dimension) {
//		std::cerr << "STM trouble: " << i << " " << i_dimension << " " << oi << " " << k << " " << k_dimension << " " << ok << std::endl;
//      }
   };
	
  const positionType k_offset;

  const lengthType Wi_min;
	
  void (shortTermMemory::* const p2calcpos)(positionType&, positionType&) const;

  const long cell_size;
  const long ptr_size;
  const long mat3d_size;
  //	long cell_count;
  //	long matWk_cell_count;
	
	pthread_mutex_t memlock;

};

TEMP_DEF_SHORTTERMMEMORY
inline shortTermMemory< TEMP_SPES_SHORTTERMMEMORY >::shortTermMemory
   (const positionType& i_dim,
    const positionType& k_dim,
    const lengthType& Wi_dim,
    arguments& argu, 
    scorematrix& sm,
    const int num_threads,
    const positionType k_off,
    const lengthType Wi_minimum) : 
     i_dimension(i_dim),
     k_dimension(k_dim),
     Wi_dimension(Wi_dim),
     Wk_dimension(lengthType(2*argu.ltOpt("-max_diff")+3)),
     Wk_zero(lengthType(argu.ltOpt("-max_diff")+1)), 
     arg(argu),
     n_threads(num_threads),
     i_position(arg.ptOpt("-j")),
     k_position(arg.ptOpt("-l")),
     seq_length1(arg.ptOpt("lenSeq1")), 
     seq_length2(arg.ptOpt("lenSeq2")), 
     wimax_init(lengthType(1)),
     k_offset(k_off),
     Wi_min(Wi_minimum),
     p2calcpos( arg.boolOpt("-global") || arg.boolOpt("realigning") || 
		arg.boolOpt("mblrealign") ? &shortTermMemory::calcPos_global : 
		&shortTermMemory::calcPos_local),
     cell_size(sizeof(cell)),
     ptr_size(sizeof(cell*)),
     mat3d_size( i_dimension * k_dimension * ((Wi_dimension +1)* ptr_size + sizeof(lengthType)))
{
	pthread_mutex_init(&memlock, NULL);
//std::cout << "STM START " << i_dimension << " " << k_dimension << " " << Wi_dimension << " " << Wk_dimension << " " << this << std::endl;
	if (wimax_init >= Wi_dimension) {wimax_init = Wi_dimension -1;}

	if (Wk_dimension > seq_length1 + seq_length2 +3) {
		Wk_dimension = seq_length1 + seq_length2 +3;
		Wk_zero = lengthType(seq_length1 +1);
	}

	// The out of memory error handling could probally be improved.
	std::string error = "Could not allocate the shortTermMemory memory stack. Most likely cause: Out of memory.";
	try {
		// The first dimension of the matrix is only two. The current start and
		// the next

		error = "Could not allocate shortTerm memory. Most likely cause: Out of memory";

		matrix = new stmWiWkSubMatrix< cell >**[i_dimension];
		if (matrix == 0) {throw exception(error, false);}

		for(positionType i = 0; i < i_dimension; i++) {

			// The second dimension is the chunk_size
			matrix[i] = new stmWiWkSubMatrix< cell >*[k_dimension];

			if (matrix[i] == 0) {throw exception(error, false);}

			for( positionType k = 0; k < k_dimension; k++ ) {

				// There has to be room for lambda positions in the third dimension
				// plus the 0 value.
				matrix[i][k] = new stmWiWkSubMatrix<cell>(Wi_dimension, Wi_minimum,
													wimax_init, Wk_dimension, Wk_zero);
				if (matrix[i][k] == 0) {throw exception(error, false);}
			}
		}
	}
  catch ( exception ) {throw;}
  catch ( ... ) {
    throw exception(error, false);
  }

  //I have seen some slip-shod backwater burgs, but this place takes the cake.
}

TEMP_DEF_SHORTTERMMEMORY
inline void shortTermMemory< TEMP_SPES_SHORTTERMMEMORY >::transfer() {

	pthread_mutex_lock(&memlock);
	//std::cout << "stm: transfer: Must clean: " << i_position << " (matrix 0)" << std::endl;
	// Switch the two I-positions
	stmWiWkSubMatrix<cell>** tmp_matrix = matrix[0];
	for(positionType i=0; i < i_dimension -1; i++) {
		matrix[i] = matrix[i+1];
	}
	matrix[i_dimension-1] = tmp_matrix;

	// All cells matrix[i_dimension-1][k] have been cleared on the runK function. (clear_old_cells function)
	i_position--;
	pthread_mutex_unlock(&memlock);
}

TEMP_DEF_SHORTTERMMEMORY
inline void shortTermMemory< TEMP_SPES_SHORTTERMMEMORY >::clear_old_cells(int cells_to_clean) {
	//std::cout << "stm: clear_old_cells: Cleaning: " << cells_to_clean << " (matrix " << i_position - cells_to_clean << ")\n";

	// Delete all the old i==1 cells
	pthread_mutex_lock(&memlock);
	const positionType i = i_position - cells_to_clean;
	stmWiWkSubMatrix<cell>** tmp_mtx = matrix[i];
	pthread_mutex_unlock(&memlock);

	for (positionType k=0; k < k_dimension; k++) {
		tmp_mtx[k]->reset();
	}
}

TEMP_DEF_SHORTTERMMEMORY
inline shortTermMemory< TEMP_SPES_SHORTTERMMEMORY >::~shortTermMemory() {
//std::cout << "STM GONE " << this << std::endl;
	pthread_mutex_destroy(&memlock);

	for(positionType i = 0; i < i_dimension; i++) {
		for( positionType k = 0; k < k_dimension; k++ ) {
			delete matrix[i][k];
		}
		delete[] matrix[i];
	}
	delete[] matrix;
}
#endif /* SHORTTEMMEMORY */
