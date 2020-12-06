#ifndef MATRIX_H
#define MATRIX_H

#include<iostream>

struct Mat_T{
    float *Mat;
    int ref = 0;
};

Mat_T* Adder(Mat_T*, Mat_T*, int);
Mat_T* Minus(Mat_T* , Mat_T*, int);
Mat_T* Multiply(Mat_T*, Mat_T*, int, int, int, int);

class Matrix{
private:
    int row;
    int column;
    Mat_T *Mat;
public:
    Matrix();
    Matrix(int, int);
    Matrix(const Matrix &);
//    Matrix(Matrix * const );
    Matrix& operator =(const Matrix&);
    Matrix operator +(const Matrix &) const;
    Matrix operator -(const Matrix &) const;
    Matrix operator *(const Matrix &) const;
    void Initial_Matrix();
    void clone(const Matrix &);
    void self_copy();
    void resetMat(int value, int col, int row);
    int getRow();
    int getCol();
    float *getMat();
    friend Matrix operator *(float scale, const Matrix& other);
    friend Matrix operator *(const Matrix& other, float scale);
    friend std::ostream & operator << (std::ostream&, const Matrix&);
    friend Matrix & operator >>(Matrix &m1, Matrix &m2);
    ~Matrix();
};

#endif // MATRIX_H
