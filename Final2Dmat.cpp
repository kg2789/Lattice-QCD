#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

template <typename T>
class Final2Dmat{
    int row;
    int column;
    int size;
    T** mat;
    public:
    Final2Dmat(int _row, int _column){
        this->row = _row;
        this->column = _column;
        mat = new T*[column];
        for ( int i = 0; i < column; i++ ){
	        mat[i] = new T[row];
        }
    }
    void set( const T& t, const int& x, const int& y){
	    mat[y][x] = t;
    }

    T get( const int& x, const int& y ){
	    return mat[y][x];
    }
};

