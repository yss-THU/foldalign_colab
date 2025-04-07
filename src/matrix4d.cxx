#ifndef MATRIX4D
#define MATRIX4D

#include "foldalign.hxx"
#include <iostream>

template< class data>
class matrix4d {
public:
	matrix4d(const int size_);

	matrix4d& operator=(const matrix4d& m){
		if (this != &m) {
			clean();
			copy(m);
		}
		return *this;
	};

	matrix4d(const matrix4d& m) {
		copy(m);
	};

	~matrix4d() {
		clean();
	};

	void init(data value);

	data get(const int i, const int j, const int k, const int l) const {
		return matrix[i][j][k][l];
	}

	void set(const int i, const int j, const int k, const int l, data value) {
		matrix[i][j][k][l] = value;
	}

	int getSize() const {
		return size;
	}

	void print_matrix() {
		for(int i=0; i<size; i++) {
			for(int j=0; j < size; j++) {
				for(int k=0; k<size; k++) {
					for(int l=0; l<size; l++) {
						std::cout << "  " << matrix[i][j][k][l];
					}
				}
				std::cout << std::endl;
			}
		}
		std::cout << std::endl;
		std::cout << std::endl;

	}

private:

	matrix4d() {
		std::cerr << "No default constructor for matrix4d.cxx" << std::endl;
	}
	void makeMatrix();

	void copy(const matrix4d& m);
	void clean();

	int size;
	data**** matrix;
};

template<class data>
matrix4d<data>::matrix4d(const int size_) : size(size_) {
	makeMatrix();
}

template<class data>
void matrix4d<data>::init(data value) {
	for(int i=0; i<size; i++) {
		for(int j=0; j < size; j++) {
			for(int k=0; k<size; k++) {
				for(int l=0; l<size; l++) {
					matrix[i][j][k][l] = value;
				}
			}
		}
	}
}


template<class data>
void matrix4d<data>::copy(const matrix4d& m) {

	size = m.getSize();

	matrix = new data***[size];
	for(int i=0; i<size; i++) {
		matrix[i] = new data**[size];
		for(int j=0; j < size; j++) {
			matrix[i][j] = new data*[size];
			for(int k=0; k<size; k++) {
				matrix[i][j][k] = new data[size];
				for(int l=0; l<size; l++) {
					matrix[i][j][k][l] = m.get(i,j,k,l);
				}
			}
		}
	}
}

template<class data>
void matrix4d<data>::makeMatrix() {
	matrix = new data***[size];
	for(int i=0; i<size; i++) {
		matrix[i] = new data**[size];
		for(int j=0; j < size; j++) {
			matrix[i][j] = new data*[size];
			for(int k=0; k < size; k++) {
				matrix[i][j][k] = new data[size];
			}
		}
	}

}

template<class data>
void matrix4d<data>::clean() {

	for(int i=0; i<size; i++) {
		for(int j=0; j < size; j++) {
			for(int k=0; k<size; k++) {
				delete[] matrix[i][j][k];
			}
			delete[] matrix[i][j];
		}
		delete[] matrix[i];
	}
	delete[] matrix;
}
#endif
