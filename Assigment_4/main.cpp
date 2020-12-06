#include <iostream>
#include "Matrix.h"
using namespace std;

void test(Matrix a){
    cout << a;
}

int main(){
    cout << "------------Test to initial Matrix----------" << endl;
    cout << "Enter row column" << endl;
    int row, column;
    cin >> row >> column;
    Matrix a(row, column);
    a.Initial_Matrix();
    cout << "------using cout << Mat to show Matrix------" << endl;
    cout << a;
    cout << "-----------Test 2 * Mat or Mat * 2-----------" << endl;
    Matrix c = a * 2.0;
    cout << c;
    cout << 2 * a;
    cout << "--------Test Mat * Mat(using openblas)-------" << endl;
    cout << a * c;
    cout << "---Test function(Mat delivery) && destructor---" << endl;
    test(c);
    cout << "-----------------Test resetMat--------------" << endl;
    c.resetMat(7,row/2,column/2);
    cout << c;
    cout << "---------------Test Mat + Mat---------------" << endl;
    cout << (a + c);
    cout << "--------------Test Mat = Mat----------------" << endl;
    Matrix d = c;
    cout << d;
    cout << "--------------Test Mat >> Mat---------------" << endl;
    d >> c;
    cout << c;
    cout << "-----------Test getRow()/getCol()----------" << endl;
    cout << c.getCol() << " " << c.getRow() << endl;
    cout << "-------------Test Mat(Mat &)---------------" << endl;
    Matrix x(a);
    cout << x;
    cout << "----------------Test End------------------" << endl;
}
