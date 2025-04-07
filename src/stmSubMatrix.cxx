#ifndef STMSUBMATRIX
#define STMSUBMATRIX

#include "arguments.cxx"
#include "helper.cxx"
#include "scorematrix.cxx"
#include "exception.cxx"
#include "foldalign.hxx"
#include "shortTermMemory.cxx"
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

#define TEMP_DEF_STMSUBMATRIX template< class cell >
#define TEMP_SPES_STMSUBMATRIX cell


TEMP_DEF_STMSUBMATRIX
class stmSubMatrix {
public:

  inline stmSubMatrix(
			shortTermMemory< cell >* const shortTM,
			positionType i_pos,
			positionType k_pos)
			: stm(shortTM),
  		  i_position(i_pos),
			  k_position(k_pos),
			  Wk_zero(stm->getWk_zero()), 
			  Wi_min(stm->getWi_minimum()),
			  k_offset(stm->getK_offset()) {

		initSubMatrix();
	//I have seen some slip-shod backwater burgs, but this place takes the cake.
	};
	
	inline void newK(const positionType k) {
		k_position = k;
		initSubMatrix();
	}

  inline cell* getPos(positionType i, positionType k,
		      lengthType Wi, lengthType Wk, const int thread = 0) {

		calcpos(i, k);
    return matrix[i][k]->getPos(Wi,Wk);
  }
	
  inline cell* putPos(positionType i, positionType k, 
		      lengthType Wi, lengthType Wk, const int thread = 0) {
		calcpos(i, k);
		return matrix[i][k]->putPos(Wi, Wk);
	};


  inline lengthType getWimax(positionType i, positionType k,
			     const int thread = 0) {
		calcpos(i, k);
		return matrix[i][k]->getWimax();
  };

  inline long getSize() const {
    return 0;
    //    return char_size*(cell_count * cell_size + mat3d_size +
    //		      matWk_cell_count * ptr_size)/mem_scale;
  }
	
  inline ~stmSubMatrix() {
		matrix[0][0] = 0;
		matrix[1][0] = 0;
		matrix[0][1] = 0;
		matrix[1][1] = 0;
	};

private:
	stmSubMatrix();
  // The data matrix
  stmWiWkSubMatrix< cell >* matrix[2][2];

	shortTermMemory< cell >* stm;

  // The matrixs position on the I-sequence
  const positionType i_position;
	
  positionType k_position;

  const lengthType Wk_zero;
	
  const lengthType Wi_min;

  const positionType k_offset;
	
  inline void calcpos(positionType& i, positionType& k) const {
//    positionType oi = i; positionType ok = k;
    i = i_position - i;
    k = k_position - k;
//    if (i < 0 || i >= 2 || k < 0 || k >= 2) {
//	std::cerr << "STM trouble: " << i << " " << 2 << " " << oi << " " << k << " " << 2 << " " << ok << std::endl;
//    }
  };
	
	inline void initSubMatrix() {
		stm->lock();

//std::cerr << i_position << " " << k_position << std::endl;
		matrix[0][0] = stm->getSubMatrix(i_position, k_position);
		if (i_position > 0) {
			matrix[1][0] = stm->getSubMatrix(i_position-1, k_position);
		}
		else {
			matrix[1][0] = 0;
		}
		if (k_position > 0) {
			matrix[0][1] = stm->getSubMatrix(i_position, k_position-1);
			if (i_position > 0) {
				matrix[1][1] = stm->getSubMatrix(i_position-1, k_position-1);
			}
			else {
				matrix[1][1] = 0;
			}
		}
		else {
			matrix[0][1] = 0;
			matrix[1][1] = 0;
		}
//std::cerr << i_position << " " << k_position << " " <<  matrix[0][0] << " " << matrix[1][0] << " " << matrix[0][1] << " " << matrix[1][1] << std::endl;
		stm->unlock();

	}

};
#endif /* STMSUBMATRIX */
