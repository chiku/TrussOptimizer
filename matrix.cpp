#include <iostream>
#include <cmath>

using namespace std;

#ifndef __MATRIX_H__
#define __MATRIX_H__

const int MAX2 = 30;

class Matrix;

// matrix addition
Matrix operator +(const Matrix& A, const Matrix& B);
Matrix operator -(const Matrix& A, const Matrix& B);
Matrix operator *(const Matrix& A, const Matrix& B);

ostream& operator <<(ostream& s, Matrix& A);
istream& operator >>(istream& s, Matrix& A);

void zeroMatrix(Matrix& A);

/* IMPLEMENTING MATRIX */
class Matrix
{
	double mat[MAX2][MAX2];
	long int row, col;

public:
	long int rows() { return row; };
	long int columns() { return col; };

    friend Matrix operator + (const Matrix& A, const Matrix& B);
    friend Matrix operator - (const Matrix& A, const Matrix& B);
    friend Matrix operator * (const Matrix& A, const Matrix& B);

    friend ostream& operator <<(ostream& s, Matrix& A);
    friend istream& operator >>(istream& s, Matrix& A);

    friend void zeroMatrix(Matrix& A);

    Matrix() { row = 0; col = 0; }
    Matrix(long int r, long int c) { row = r; col = c; }
    Matrix operator + () const;
    Matrix operator - () const;
    Matrix scalarMult(double c);

    double& operator()(long int r, long int c);
    double operator()(long int r, long int c) const;

	void mutateToInclude(double value, long int r, long int c);
	void setSize(long int r, long int c);

	Matrix invert();
};

inline double Matrix::operator() (long int r, long int c) const
{
	return mat[r][c];
}

inline double& Matrix::operator() (long int r, long int c)
{
	return mat[r][c];
}

inline void Matrix::mutateToInclude(double value, long int r, long int c)
{
	mat[r][c] = value;
	setSize(max(row, r+1), max(col, c+1));
}

inline void Matrix::setSize(long int r, long int c)
{
	row = r;
	col = c;
}

// Make all elements of the matrix 0
inline void zeroMatrix(Matrix& A)
{
	for (int i=0; i<A.row; i++)
		for (int j=0; j<A.col; j++)
			A.mat[i][j] = 0;
}


// matrix addition
Matrix operator +(const Matrix& A, const Matrix& B)
{
	Matrix C;
	C.setSize(A.row, A.col);
	if (A.row!=B.row || A.col!=B.col)
	{
		cout <<"Matries incompatible for addition. (" <<A.row <<", " <<A.col <<") and ("<< B.row <<", " <<B.col <<")\n";
		return C;
	}
	for (int i=0; i<A.row; i++)
		for (int j=0; j<A.col; j++)
			C.mat[i][j] = A.mat[i][j]+B.mat[i][j];
	return A;
}

// matrix subtraction
Matrix operator -(const Matrix& A, const Matrix& B)
{
	Matrix C;
	C.setSize(A.row, A.col);
	if (A.row!=B.row || A.col!=B.col)
	{
		cout <<"Matries incompatible for subtraction. (" <<A.row <<", " <<A.col <<") and ("<< B.row <<", " <<B.col <<")\n";
		return C;
	}
	for (int i=0; i<A.row; i++)
		for (int j=0; j<A.col; j++)
			C.mat[i][j]=A.mat[i][j]-B.mat[i][j];
	return C;
}

// matrix multiplication
Matrix operator *(const Matrix& A, const Matrix& B)
{
	Matrix C;
	C.setSize(A.row, B.col);
	zeroMatrix(C);
	if (A.col!=B.row)
	{
		cout <<"Matries incompatible for multiplication. (" <<A.row <<", " <<A.col <<") and ("<< B.row <<", " <<B.col <<")\n";
		return C;
	}
	int i, j, k;
	for (i=0; i<A.row; i++)
			for (j=0; j<B.col; j++)
				for (k=0; k<A.col; k++)
				C.mat[i][j] += A.mat[i][k] * B.mat[k][j];

	return C;
}

// unary +
Matrix Matrix::operator +() const
{
	return (*this);
}

// unary -
Matrix Matrix::operator -() const
{
	Matrix B;
	B.setSize(row, col);
	for (int i=0; i<row; i++)
		for (int j=0; j<col; j++)
			B.mat[i][j] = -mat[i][j];
	return B;
}

// scalar multiplicaiton
Matrix Matrix::scalarMult(double x)
{
	Matrix B;
	B.setSize(row, col);
	for (int i=0; i<row; i++)
		for (int j=0; j<col; j++)
			B.mat[i][j] = x*mat[i][j];
	return B;
}

// matrix inversion guass jordan elimination using full pivoting
Matrix Matrix::invert()
{
	int n = row;
	if (row != col)
	{
		cout <<"Matrix not square.(" <<row <<", " <<col <<")\n";
		return (*this);
	}
	
	int i, icol, irow, j, k, l, ll;
	float big, dum, pivinv;
	int* indxc = new int[n];
	int* indxr = new int[n]; 
	int* ipiv = new int[n]; /* The integer arrays ipiv, indxr, and indxc are
	used for bookkeeping on the pivoting. */

	for (i=0; i<n; i++) 
		ipiv[i]=-1;
		
	for (i=0; i<n; i++) 
	{ 
		// This is the main loop over the columns to be reduced.
		big = 0.0; 
		for (j=0; j<n; j++) // This is the outer loop of the search for a pivot element.
			if (ipiv[j] != 0) 
				for (k=0; k<n; k++)
				{
					if (ipiv[k] == -1)
					{
						if (fabs(mat[j][k]) >= big)
						{
							big = fabs(mat[j][k]);
							irow = j;
							icol = k;
						}
					} 
					else if (ipiv[k] > 1) 
					{
						cout <<"Singular Matrix.\n";
						return (*this);
					}
				}
		(ipiv[icol])++;
		/* We now have the pivot element, so we interchange rows, if needed, to put the pivot
		element on the diagonal. The columns are not physically interchanged, only relabeled:
		indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
		indxr[i] is the row in which that pivot element was originally located. If indxr[i] 6=
		indxc[i] there is an implied column interchange. With this form of bookkeeping, the
		solution b's will end up in the correct order, and the inverse matrix will be scrambled
		by columns. */
		if (irow != icol)
		{
			for (l=0; l<n; l++)
				swap(mat[irow][l],mat[icol][l]);
		}
	
		indxr[i] = irow; // We are now ready to divide the pivot row by the
		indxc[i] = icol; // pivot element, located at irow and icol.
		if (mat[icol][icol] == 0.0) 
		{
			cout <<"Singular Matrix.\n";
			return (*this);
		}
		pivinv = 1.0/mat[icol][icol];
		mat[icol][icol]=1.0;
		for (l=0; l<n; l++)
			mat[icol][l] *= pivinv;

		for (ll=0; ll<n; ll++) // Next, we reduce the rows...
			if (ll != icol) //  ...except for the pivot one, of course.
			{	
				dum = mat[ll][icol];
				mat[ll][icol] = 0.0;
				for (l=0; l<n;l++) 
					mat[ll][l] -= mat[icol][l]*dum;
			}
	}
	/* This is the end of the main loop over columns of the reduction. It only remains to unscramble
	the solution in view of the column interchanges. We do this by interchanging pairs of
	columns in the reverse order that the permutation was built up. */
	for (l=n-1; l>=0; l--)
	{
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				swap(mat[k][indxr[l]],mat[k][indxc[l]]);
	} // And we are done.
	delete[] indxc;
	delete[] indxr;
	delete[] ipiv;

	return (*this);
}


// accept matrix from console
istream& operator >>(istream& s, Matrix& A)
{
	for (int i=0; i<A.row; i++)
		for (int j=0; j<A.col; j++)
			s >>A.mat[i][j];
	return s;
}

// prints matrix to console
ostream& operator <<(ostream& s, Matrix& A)
{
	for (int i=0; i<A.row; i++)
	{
		s <<endl;
		for (int j=0; j<A.col; j++)
			s <<A.mat[i][j] <<'\t';
	}
	return s;
}

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
