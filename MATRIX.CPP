#include <iostream>
#include <cmath>
#include <conio.h>

using namespace std;

const double TINY = 1e-10;

/* IMPLEMENTING MATRIX */

struct Matrix
{
	double mat[13][13];
	int row, col;
};

// changes order of matrix
void changeOrder(Matrix& A, int r, int c)
{
	A.row=r;
	A.col=c;
}



// matrix addition
Matrix operator +(Matrix a, Matrix b)
{
	Matrix c;
	changeOrder(c, a.row, a.col);
	if (a.row!=b.row || a.col!=b.col)
	{
		cout <<"Matries incompatible for addition.\n";
		return c;
	}
	for (int i=0; i<a.row; i++)
		for (int j=0; j<a.col; j++)
			c.mat[i][j]=a.mat[i][j]+b.mat[i][j];
	return c;
}

// matrix subtraction
Matrix operator -(Matrix A, Matrix B)
{
	Matrix C;
	changeOrder(C, A.row, A.col);
	if (A.row!=B.row || A.col!=B.col)
	{
		cout <<"Matries incompatible for subtraction.\n";
		return C;
	}
	for (int i=0; i<A.row; i++)
		for (int j=0; j<A.col; j++)
			C.mat[i][j]=A.mat[i][j]-B.mat[i][j];
	return C;
}

// matrix multiplication
Matrix operator *(Matrix A, Matrix B)
{
	Matrix C;
	changeOrder(C, A.row, B.col);
	if (A.col!=B.row)
	{
		cout <<"Matries incompatible for multiplication.\n";
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
Matrix operator +(Matrix A)
{
	return A;
}

// unary -
Matrix operator -(Matrix A)
{
	Matrix B;
	changeOrder(B, A.row, A.col);
	for (int i=0; i<A.row; i++)
		for (int j=0; j<A.col; j++)
			B.mat[i][j] = -A.mat[i][j];
	return B;
}

// scalar multiplicaiton
Matrix scalarMult(Matrix A, double x)
{
	Matrix B;
	changeOrder(B, A.row, A.col);
	for (int i=0; i<A.row; i++)
		for (int j=0; j<A.col; j++)
			B.mat[i][j] = x*A.mat[i][j];
	return B;
}

#define NMAX 100
/* Function to perform the partial-pivoting Gaussian elimination.
	 a[][] is the original matrix in the input and transformed
	 matrix plus the pivoting element ratios below the diagonal
	 in the output.	indx[] records the pivoting order.
 */
void pivot (Matrix& a, int indx[])
{
	int i, j, k, itmp;
	double c1, pi, pi1, pj;
	double c[NMAX];
	int n = a.row;

	if (n > NMAX)
	{
		cout <<"The matrix dimension is too large.\n";
		exit(1);
	}

	/* Initialize the index */	

	for (i = 0; i < n; ++i)
		indx[i] = i;

	/* Find the rescaling factors, one from each row */
 
	for (i = 0; i < n; ++i)
	{
		c1 = 0;
		for (j = 0; j < n; ++j)
			if (fabs(a.mat[i][j]) > c1) 
				c1 = fabs(a.mat[i][j]);
		c[i] = c1;
	}

	/* Search the pivoting (largest) element from each column */ 

	for (j = 0; j < n-1; ++j)
	{
		pi1 = 0;
		for (i = j; i < n; ++i)
		{
			pi = fabs(a.mat[indx[i]][j]) / c[indx[i]];
			if (pi > pi1)
			{
				pi1 = pi;
				k = i;
			}
		}

		/* Interchange the rows via indx[] to record pivoting order */

		itmp = (int)indx[j];
		indx[j] = indx[k];
		indx[k] = itmp;
		for (i = j+1; i < n; ++i)
		{
			pj = a.mat[indx[i]][j]/a.mat[indx[j]][j];

			/* Record pivoting ratios below the diagonal */

			a.mat[indx[i]][j] = pj;

			/* Modify other elements accordingly */

			for (k = j+1; k < n; ++k)
				a.mat[indx[i]][k] = a.mat[indx[i]][k]-pj*a.mat[indx[j]][k];
		}
	}
}

/* Function to invert matrix a[][] with the inverse stored
	 in x[][] in the output. */
void inv(Matrix& a, Matrix& x, int indx[])
{
	int n = a.row;
	int i,j,k;
	double b[NMAX][NMAX];
	x.row = x.col = a.row;

	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			b[i][j] = 0;

	for (i = 0; i < n; ++i)
		b[i][i] = 1;

	pivot(a, indx);

	for (i = 0; i < n-1; ++i)
		for (j = i+1; j < n; ++j)
			for (k = 0; k < n; ++k)
				b[indx[j]][k] = b[indx[j]][k]-a.mat[indx[j]][i]*b[indx[i]][k];


	for (i = 0; i < n; ++i)
		x.mat[n-1][i] = b[indx[n-1]][i]/a.mat[indx[n-1]][n-1];
		for (j = n-2; j >= 0; j = j-1)
		{
			x.mat[j][i] = b[indx[j]][i];
			for (k = j+1; k < n; ++k)
				x.mat[j][i] = x.mat[j][i]-a.mat[indx[j]][k]*x.mat[k][i];
			x.mat[j][i] = x.mat[j][i]/a.mat[indx[j]][j];
		}
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

int main()
{
	Matrix A;
	changeOrder(A, 3, 3);
	A.mat[0][0] = 1; A.mat[0][1] = 3; A.mat[0][2] = 4;
	A.mat[1][0] = 3; A.mat[1][1] = 2; A.mat[1][2] = 1;
	A.mat[2][0] = 7; A.mat[2][1] = 8; A.mat[2][2] = 8;
	cout <<A;
	Matrix B = A;
	int indx[100];
	inv(A, B, indx);
	cout <<endl <<B;
	cout <<endl <<B.row <<B.col;
	getch();
}
	
