#include <iostream.h>
/* IMPLEMENTING MATRIX CLASS */

template <class T>
class Matrix
{
	private:
   	T mat[20][20];
      int row, col;

      void changeOrder(int, int);

   public:
   	Matrix matAdd(Matrix, Matrix);
      Matrix matSub(Matrix, Matrix);
      Matrix matMul(Matrix, Matrix);

      T element(int, int);
      T determinant(Matrix);
      void getMatrix();
      void printMatrix();

   public:
   	Matrix();
      Matrix(int, int);
};

// no-arguement constructor
template <class T>
Matrix<T>::Matrix()
{
	row=col=1;
	mat[0][0]=0;
}

// two-arguement constructor
template <class T>
Matrix<T>::Matrix(int r, int c)
{
   if (r<=0 || c<=0)
   {
   	cerr <<"Order of matrix cannot be non-positive.\n";
      return;
   }
   row=r;
   col=c;
   for (int i=0; i<r; i++)
   	for (int j=0; j<r; j++)
      	mat[i][j]=0;
}


// matrix addition
template <class T>
Matrix Matrix<T>::matAdd(Matrix a, Matrix b)
{
	Matrix c;
   if (a.row!=b.row || a.col!=b.col)
   {
   	cerr <<"Matries incompatible for addition.\n";
      return c;
   }
   c.changeOrder(a.row, a.col);
   for (int i=0; i<a.row; i++)
   	for (int j=0; j<a.col; j++)
      	c.element(i, j)=a.element(i, j)+b.element(i, j);
   return c;
}

// matrix subtraction
template <class T>
Matrix Matrix<T>::matSub(Matrix a, Matrix b)
{
	Matrix c;
   if (a.row!=b.row || a.col!=b.col)
   {
   	cerr <<"Matries incompatible for subtraction.\n";
      return c;
   }
   c.changeOrder(a.row, a.col)
   for (int i=0; i<a.row; i++)
   	for (int j=0; j<a.col; j++)
      	c.element(i, j)=a.element(i, j)-b.element(i, j);
   return c;
}

// matrix multiplication
template <class T>
Matrix Matrix<T>::matMul(Matrix a, Matrix b)
{
	Matrix c;
   if (a.col!=b.row)
   {
   	cerr <<"Matries incpmpatible for multiplication.\n";
      return c;
   }
   c.changeOrder(a.row, b.col);
   for (int i=0; i<a.row; i++)
   	for (int j=0; j<b.col; b++)
        	for (int k=0; k<a.col; k++)
         	c[i][j] += a[i][k] + b[i][j];
   return c;
}

/*
// determines corresponding determinant
template <class T>
T Matrix<T>::determinant(T A)
{
	Matrix aux;
	T d=0;

	if (n==1)
		return A[0][0];

	if (n==2)
		return (A[0][0]*A[1][1] - A[1][0]*A[0][1]);

	for (int a=0; a<n; a++)
	{
		int i, j;
		for (i=0; i<n; i++)
		{
			for (j=0; j<a; j++)
				aux[i][j] = A[i+1][j];
		}

		for (i=0; i<n; i++)
		{
			for (j=a; j<n; j++)
				aux[i][j] = A[i+1][j+1];
		}

		if (a%2)
			d -= A[0][a]*determinant(aux, n-1);
		else
			d += A[0][a]*determinant(aux, n-1);
	}
	return d;
}    */

// accept matrix from console
template <class T>
void Matrix<T>::getMatrix()
{
	for (int i=0; i<row; i++)
      for (int j=0; j<col; j++)
      	cin >>mat[i][j];
}

// prints matrix to console
template <class T>
void Matrix<T>::printMatrix()
{
	for (int i=0; i<row; i++)
   {
      cout <<endl;
      for (int j=0; j<col; j++)
      	cout <<mat[i][j] <<'\t';
   }
}


void main()
{
   Matrix<double> A(2, 2), B(2, 2);
   A.getMatrix();
   B.getMatrix();
   A.printMatrix();
   B.printMatrix();
   Matrix<double> C(2, 2);
   C=C.matAdd(A, B);
   C.printMatrix();
}



/*
void part(double A[MAT_MAX_SIZE][MAT_MAX_SIZE], double B[MAT_MAX_SIZE][MAT_MAX_SIZE], int n, int x, int y)
{
	int i,j;
	for (i=0; i<x; i++)
	{
	for (j=0; j<y; j++)
		B[i][j] = A[i][j];
	}
	for (i=0; i<x; i++)
	{
		for (j=y; j<n; j++)
			B[i][j] = A[i][j+1];
	}
	for (i=x; i<n; i++)
	{
		for (j=0; j<y; j++)
			B[i][j] = A[i+1][j];
	}
	for (i=x; i<n; i++)
	{
		for (j=y; j<n; j++)
			B[i][j] = A[i+1][j+1];
	}
}


void invert(double A[MAT_MAX_SIZE][MAT_MAX_SIZE], double B[MAT_MAX_SIZE][MAT_MAX_SIZE], int n)
{
	double detA=determinant(A, n);
	double temp[MAT_MAX_SIZE][MAT_MAX_SIZE];
	int i, j;
	if (detA == 0)
	{
		cout <<"The matrix has no inverse!\n";
		return;
	}
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			double aux[MAT_MAX_SIZE][MAT_MAX_SIZE];
			double t;
			part(A, aux, n, i, j);
			t=determinant(aux, n-1);
			temp[i][j] = t/detA;
			if ((i+j)%2)
				temp[i][j] = -temp[i][j];
		}
	}

	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
			B[i][j] = temp[j][i];
	}
}
*/