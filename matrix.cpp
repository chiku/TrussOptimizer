#include <iostream>
#include <cmath>

#include "cmatrix/src/matrix.h"

using namespace std;

#ifndef __LOCAL_MATRIX_H__
#define __LOCAL_MATRIX_H__

typedef CMatrix::Matrix<double> Matrix;

#endif

/*
int main()
{
	Matrix A;
	A.setSize(3, 3);
	A(0, 0) = 1; A(0, 1) = 3; A(0, 2) = 4;
	A(1, 0) = 3; A(1, 1) = 2; A(1, 2) = 1;
	A(2, 0) = 7; A(2, 1) = 8; A(2, 2) = 8;
	cout <<A;
	Matrix B = A;
	
	B.invert();
	cout <<endl <<B;
	Matrix C;
	C = A*B;
	cout <<endl <<C;
	cout <<endl <<B.rows() <<B.columns();
}
*/
