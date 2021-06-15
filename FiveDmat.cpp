#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <stdexcept>

template <typename T>
class FiveDmat{
    int dim1;
    int dim2;
    int dim3;
    int dim4;
    int dim5;
    T***** FDmat;
    public:
    FiveDmat(int _dim1, int _dim2, int _dim3, int _dim4,int _dim5){
        this->dim1 = _dim1;
        this->dim2 = _dim2;
        this->dim3 = _dim3;
        this->dim4 = _dim4;
        this->dim5 = _dim5;
        FDmat = new T****[dim1];
        for(int i =0; i <dim1; i++){
            FDmat[i] = new T***[dim2];
            for (int j = 0; j < dim2; j++){
                FDmat[i][j] = new T**[dim3];
                for (int k = 0; k < dim3; k++){
                    FDmat[i][j][k] = new T*[dim4];
                    for (int l = 0; l < dim3; l++){
                        FDmat[i][j][k][l] = new T[dim5];
                    }
                }
            }
        }
    }
    void set(T value, int i, int j, int k, int l, int m){
        FDmat[i][j][k][l][m] = value;
    }

    T get(int pos1, int pos2, int pos3, int pos4, int pos5){
        return FDmat[pos1][pos2][pos3][pos4][pos5];
    }
};
