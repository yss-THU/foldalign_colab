#ifndef STMWIWKSUBMATRIX
#define STMWIWKSUBMATRIX

#include "arguments.cxx"
#include "helper.cxx"
#include "scorematrix.cxx"
#include "exception.cxx"
#include "foldalign.hxx"

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

#define TEMP_DEF_STMWIWKSUBMATRIX template< class cell >
#define TEMP_SPES_STMWIWKSUBMATRIX cell

TEMP_DEF_STMWIWKSUBMATRIX
class stmWiWkSubMatrix {
public:

	inline stmWiWkSubMatrix(const lengthType& Wi_dim,
				const lengthType Wi_minimum,
				lengthType wimax_init,
				const lengthType Wk_dim,
				const lengthType WkZero);
	
	inline cell* getPos(lengthType Wi, lengthType Wk, const int thread = 0) {

		calcpos2(Wi, Wk);
		if (matrix[Wi] == 0) {
			return 0;
		}
		return matrix[Wi][Wk];
	}
	
	inline cell* putPos(lengthType Wi, lengthType Wk, const int thread = 0);

	inline lengthType getWimax() const {
		return Wimax;
	};

	inline long getSize() const {
		return 0;
    //    return char_size*(cell_count * cell_size + mat3d_size +
    //		      matWk_cell_count * ptr_size)/mem_scale;
	}
	
	inline void reset() {
		for(lengthType Wi = 0; Wi < Wi_dimension; Wi++ ) {
			if (matrix[Wi] != 0) {
				for(lengthType Wk = 0; Wk < Wk_dimension; Wk++) {
					if (matrix[Wi][Wk] != 0) {
						delete matrix[Wi][Wk];
						matrix[Wi][Wk] = 0;
					}
				}
				delete[] matrix[Wi];
				matrix[Wi] = 0;
			}
		}
		Wimax = Wimax_init;
	};

	
	inline ~stmWiWkSubMatrix();

private:

	// The window size along sequence I dimension
	const lengthType Wi_dimension;
	const lengthType Wi_min;
	lengthType Wimax;
	lengthType Wimax_init;

	// The window size along sequence K dimension
	// This is not const since it is depended on the Wi. In the extremes of Wi
	// it is not const.
	const lengthType Wk_dimension;
	const lengthType Wk_zero;

	// The data matrix
	cell*** matrix;

	inline void calcpos2(lengthType& Wi, lengthType& Wk) const {
//		positionType oWi = Wi; positionType oWk = Wk;
		Wk = Wk - Wi + Wk_zero;
		Wi -= Wi_min;
//		if (Wi < 0 || Wi >= Wi_dimension || Wk < 0 || Wk >= Wk_dimension) {
//			std::cerr << "STM WiWk trouble: " << Wi << " " << Wi_dimension << " " << oWi << " " << Wk << " " << Wk_dimension << " " << oWk << std::endl;
//		}
	};

	inline void newCell (const lengthType Wi, const lengthType Wk) {
   	matrix[Wi][Wk] = new cell();
	}

};

TEMP_DEF_STMWIWKSUBMATRIX
inline stmWiWkSubMatrix< TEMP_SPES_STMWIWKSUBMATRIX >::stmWiWkSubMatrix
			(const lengthType& Wi_dim,
			const lengthType Wi_minimum,
			lengthType wimax_init,
			const lengthType Wk_dim,
			const lengthType WkZero) 
 : Wi_dimension(Wi_dim - Wi_minimum),
   Wi_min(Wi_minimum),
	Wimax(wimax_init),
	Wimax_init(wimax_init),
   Wk_dimension(Wk_dim),
   Wk_zero(WkZero)
{

	// The out of memory error handling could probally be improved.
	std::string error = "Could not allocate stmWiWkSubMatrix. Most likely cause: Out of memory";
	try {
		// The first dimension of the matrix is only two. The current start and
		// the next


		matrix = new cell**[Wi_dimension];
		if (matrix == 0) {throw exception(error, false);}

		for(lengthType Wi = 0; Wi < Wi_dimension; Wi++ ) {
			matrix[Wi] = 0;
		}
	}
	catch ( exception ) {throw;}
	catch ( ... ) {
		throw exception(error, false);
	}
	//I have seen some slip-shod backwater burgs, but this place takes the cake.
}

TEMP_DEF_STMWIWKSUBMATRIX
inline cell* stmWiWkSubMatrix< TEMP_SPES_STMWIWKSUBMATRIX >::putPos(
		lengthType Wi, lengthType Wk, const int thread) {

	if (Wi > Wimax) {Wimax = Wi;}
	
	calcpos2(Wi, Wk);

	if (matrix[Wi] == 0) {

		matrix[Wi] = new cell*[Wk_dimension];

		for ( lengthType tWk = 0; tWk < Wk_dimension; tWk++ ) {
			// Initially every slot is empty
			matrix[Wi][tWk] = 0;
		}
		matrix[Wi][Wk] = new cell();
	}
	else {
		if (matrix[Wi][Wk] == 0) {
			matrix[Wi][Wk] = new cell();
		}
	}

	return matrix[Wi][Wk];
}

TEMP_DEF_STMWIWKSUBMATRIX
inline stmWiWkSubMatrix< TEMP_SPES_STMWIWKSUBMATRIX >::~stmWiWkSubMatrix() {

	for( lengthType Wi = 0; Wi < Wi_dimension; Wi++ ) {
		if (matrix[Wi] != 0) {
			for ( lengthType Wk = 0; Wk < Wk_dimension; Wk++ ) {
				if ( matrix[Wi][Wk] != 0 ) {
					delete matrix[Wi][Wk];
				}
			}
			delete[] matrix[Wi];
		}
	}
	delete[] matrix;
}
#endif /* STMWIWKSUBMATRIX */
