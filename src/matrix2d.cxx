#ifndef MATRIX2D
#define MATRIX2D


#include "foldalign.hxx"
#include <iostream>

template< class data>
class matrix2d {
public:
	matrix2d(const int size_);

	matrix2d& operator=(const matrix2d& m){
		if (this != &m) {
			clean();
			copy(m);
		}
		return *this;
	};

	matrix2d(const matrix2d& m) {
		copy(m);
	};

	~matrix2d() {
		clean();
	};

	void init(data value);

	data get(const int i, const int j) const {
		return matrix[i][j];
	}

	void set(const int i, const int j, data value) {
		matrix[i][j] = value;
	}

	int getSize() const {
		return size;
	}

	void print_matrix() {
		for(int i=0; i< size; i++) {
			for (int j=0; j<size; j++) {
				std::cout << "  " << matrix[i][j];
			}
			std::cout << std::endl;
		}
	}
private:

	matrix2d() {
		std::cerr << "No default constructor for matrix2d.cxx" << std::endl;
	}
	void makeMatrix();

	void copy(const matrix2d& m);
	void clean();

	int size;
	data** matrix;
};

template<class data>
matrix2d<data>::matrix2d(const int size_) : size(size_) {
	makeMatrix();
}

template<class data>
void matrix2d<data>::init(data value) {

	for(int i=0; i<size; i++) {
		for(int j=0; j < size; j++) {
			matrix[i][j] = value;
		}
	}
}


template<class data>
void matrix2d<data>::copy(const matrix2d& m) {

	size = m.getSize();

	matrix = new data*[size];
	for(int i=0; i<size; i++) {
		matrix[i] = new data[size];
		for(int j=0; j < size; j++) {
			matrix[i][j] = m.get(i,j);
		}
	}
}

template<class data>
void matrix2d<data>::makeMatrix() {
	matrix = new data*[size];
	for(int i=0; i<size; i++) {
		matrix[i] = new data[size];
	}
}

template<class data>
void matrix2d<data>::clean() {

	for(int i=0; i<size; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
}
#endif
