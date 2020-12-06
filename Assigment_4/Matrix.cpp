#include <iostream>
#include "Matrix.h"
#include "/home/bill/桌面/openblas/include/cblas.h"
using namespace std;

Matrix::Matrix(int a, int b){
    column = b;
    row = a;
    Mat = new Mat_T;
    Mat->Mat = new float[a * b]();
}

Matrix::Matrix(){
    column = 0;
    row = 0;
    Mat = new Mat_T;
    Mat->Mat = new float[0]();
}

Matrix::Matrix(const Matrix &a){
    column = a.column;
    row = a.row;
    Mat = a.Mat;
    a.Mat->ref ++;
}

Matrix& Matrix::operator=(const Matrix &M){
    if(Mat->ref <= 0){
        delete Mat;
    }
    else{
        Mat->ref--;
    }
    column = M.column;
    row = M.row;
    Mat = M.Mat;
    M.Mat->ref ++;
    return *this;
}


Matrix Matrix::operator+(const Matrix &other) const{

    if(column != other.column || row != other.row){
        cout << "The Matrix are not match!" << endl;
        cout << "Using default constructors ..." << endl;
        return Matrix();
    }
    Matrix *b = new Matrix(column, row);
    b->Mat->ref--;
    b->Mat = Adder(Mat, other.Mat, column * row);

    return *b;
}

Matrix Matrix::operator-(const Matrix &other) const{

    if(column != other.column || row != other.row){
        cout << "The Matrix are not match!" << endl;
        cout << "Using default constructors ..." << endl;
        return Matrix();
    }

    Matrix *b = new Matrix(column, row);
    b->Mat->ref--;
    b->Mat = Minus(Mat, other.Mat, column * row);

    return *b;
}

Matrix Matrix::operator*(const Matrix &other) const{

    if(column != other.row){
        cout << "The Matrix are not match!" << endl;
        cout << "Using default constructors ..." << endl;
        return Matrix();
    }

    Matrix *b = new Matrix(column, row);
    b->Mat->ref--;
    b->Mat = Multiply(Mat, other.Mat, column * row, row, other.column, column);

    return *b;
}

Mat_T* Multiply(Mat_T* a, Mat_T* b, int c, int row, int column, int k){
    Mat_T *out = new Mat_T;
    out->Mat = new float[c];
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row , column, k, 1.0, a->Mat, k, b->Mat, column, 0.0, out->Mat, column);
    return out;
}

Mat_T* Minus(Mat_T* a, Mat_T* b, int length){
        Mat_T* to_ret = new Mat_T;
        to_ret->Mat = new float[length];

        for(int i = 0; i < length; i++){
            to_ret->Mat[i] = a->Mat[i] - b->Mat[i];
        }
        return to_ret;
}

Mat_T* Adder(Mat_T* a, Mat_T* b, int length){
        Mat_T* to_ret = new Mat_T;
        to_ret->Mat = new float[length];

        for(int i = 0; i < length; i++){
            to_ret->Mat[i] = a->Mat[i] + b->Mat[i];
        }
        return to_ret;
}

void Matrix::Initial_Matrix(){
    cout << "Initial Start" << endl;
    int length = column * row;
    for(int i = 0; i < length; i++){
       cin >> Mat->Mat[i];
    }
    cout << "Initial end" << endl;
}

std::ostream & operator<< (std::ostream& os, const Matrix& other){
    int length = other.column * other.row;
    for(int i = 0; i < length; i++){
        if(i % other.column == 0){
            os << endl;
        }
        os << other.Mat->Mat[i] << " ";
    }
    os << endl;
    return os;
}

Matrix & operator>>(Matrix &m1, Matrix &m2){
    m1.clone(m2);
    return m1;
}


Matrix::~Matrix(){
    if(Mat->ref <= 0){
        delete Mat;
        cout << "delete completely" << endl;
    }
    else{
        Mat->ref--;
        cout << "No delete completely" << endl;
    }
}

Matrix operator*(float scale, const Matrix& other){
    Matrix *a = new Matrix(other.row, other.column);
    for(int i = 0; i < other.column * other.row; i++){
        a->Mat->Mat[i] = other.Mat->Mat[i] * scale;
    }
    a->Mat->ref--;
    return *a;
}

Matrix operator*(const Matrix& other, float scale){
    Matrix *a = new Matrix(other.row, other.column);
    for(int i = 0; i < other.column * other.row; i++){
        a->Mat->Mat[i] = other.Mat->Mat[i] * scale;
    }
    a->Mat->ref--;
    return *a;
}

void Matrix::clone(const Matrix &a){
    this->Mat->ref--;
    if(this->Mat->ref < 0){
        cout << "Old Mat delete!" << endl;
        delete this->Mat;
    }

    column = a.column;
    row = a.row;
    Mat = new Mat_T;
    Mat->Mat = new float[column * row];
    for(int i = 0; i < a.column * a.row; i++){
        Mat->Mat[i] = a.Mat->Mat[i];
    }
    cout << "clone completely!" << endl;
}

void Matrix::self_copy(){
    float *tmp = new float[column * row];
    for(int i = 0; i < column * row; i++){
        tmp[i] = Mat->Mat[i];
    }
    Mat->Mat = tmp;
    cout << "clone completely!" << endl;
}

void Matrix::resetMat(int value, int rowc, int col){
//    cout << Mat->ref << endl;
    if(Mat->ref == 0){
        Mat->Mat[(rowc -1) * column + col - 1] = value; //考虑到人的习惯
        cout << "set successfully!" << endl;
    }
    else{
        Mat->ref--;
        self_copy();
        Mat->Mat[(rowc - 1) * column + col - 1] = value;
        cout << "set successfully!" << endl;
    }
}

int Matrix::getRow(){
    return row;
}

int Matrix::getCol(){
    return column;
}
float * Matrix::getMat(){
    return Mat->Mat;
}





















