#include <iostream>
#include <strstream>

/* IMPLEMENTING MATRIX CLASS */

class Matrix
{
	private:
   	double mat[20][20];
      int row, col;

      void changeOrder(int, int);
      Matrix part(int, int);

   public:
     	Matrix();
      Matrix(int, int);
      Matrix(int, int, char[]);

   	double element(int, int);
      double determinant();

      friend Matrix operator +(Matrix, Matrix);
      friend Matrix operator -(Matrix, Matrix);
      friend Matrix operator *(Matrix, Matrix);

      friend Matrix operator +(Matrix);
      friend Matrix operator -(Matrix);

      Matrix scalarMult(double);
      friend Matrix operator ~(Matrix);

      friend ostream& operator <<(ostream& , Matrix&);
      friend istream& operator >>(istream& , Matrix&);
};

// no-arguement constructor
Matrix::Matrix()
{
	row=col=1;
	mat[0][0]=0;
}

// two-arguement constructor
Matrix::Matrix(int r, int c)
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

// three-arguement constructor
Matrix::Matrix(int r, int c, char str[])
{
   if (r<=0 || c<=0)
   {
   	cerr <<"Order of matrix cannot be non-positive.\n";
      return;
   }
   row=r;
   col=c;

   istrstream s(str);
   for (int i=0; i<r; i++)
   	for (int j=0; j<c; j++)
      	s >>mat[i][j];
}


// returns element
double Matrix::element(int r, int c)
{
	return mat[r][c];
}

// changes order of matrix
void Matrix::changeOrder(int r, int c)
{
	row=r;
   col=c;
}

// determines corresponding determinant
double Matrix::determinant()
{
	Matrix aux;
	double d=0;

   int n;
   if (row!=col)
   {
   	cerr <<"\nOnly square matrices can have corresponding determinants.";
      return d;
   }
   n=row;

   if (n==1)
		return mat[0][0];

	if (n==2)
		return (mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1]);

   aux.changeOrder(n-1, n-1);
	for (int a=0; a<n; a++)
	{
		int i, j;
		for (i=0; i<n; i++)
			for (j=0; j<a; j++)
				aux.mat[i][j] = mat[i+1][j];

		for (i=0; i<n; i++)
			for (j=a; j<n; j++)
				aux.mat[i][j] = mat[i+1][j+1];

		if (a%2)
			d -= mat[0][a]*aux.determinant();
		else
			d += mat[0][a]*aux.determinant();
	}
	return d;
}

// used in finding the inverse of matrix
Matrix Matrix::part(int m, int n)
{
   Matrix B;
   B.changeOrder(row-1, col-1);
	int i,j;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			B.mat[i][j] = element(i, j);
	for (i=0; i<m; i++)
		for (j=n; j<col; j++)
			B.mat[i][j] = element(i, j+1);

	for (i=m; i<row; i++)
		for (j=0; j<n; j++)
			B.mat[i][j] = element(i+1, j);
	for (i=m; i<row; i++)
		for (j=n; j<col; j++)
			B.mat[i][j] = element(i+1, j+1);
   return B;
}


// matrix addition
Matrix operator +(Matrix a, Matrix b)
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
      	c.mat[i][j]=a.element(i, j)+b.element(i, j);
   return c;
}

// matrix subtraction
Matrix operator -(Matrix A, Matrix B)
{
	Matrix C;
   if (A.row!=B.row || A.col!=B.col)
   {
   	cerr <<"Matries incompatible for subtraction.\n";
      return C;
   }
   C.changeOrder(A.row, A.col);
   for (int i=0; i<A.row; i++)
   	for (int j=0; j<A.col; j++)
      	C.mat[i][j]=A.element(i, j)-B.element(i, j);
   return C;
}

// matrix multiplication
Matrix operator *(Matrix A, Matrix B)
{
	Matrix C;
   if (A.col!=B.row)
   {
   	cerr <<"Matries incpmpatible for multiplication.\n";
      return C;
   }
   C.changeOrder(A.row, B.col);
   int i, j, k;
   for (i=0; i<A.row; i++)
	  	for (j=0; j<B.col; j++)
        	for (k=0; k<A.col; k++)
         	C.mat[i][j] += A.element(i, k) * B.element(i, j);

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
	for (int i=0; i<A.row; i++)
   	for (int j=0; j<A.col; j++)
      	A.mat[i][j] = -A.mat[i][j];
	return A;
}

// scalar multiplicaiton
Matrix Matrix::scalarMult(double x)
{
   Matrix B;
   B.changeOrder(row, col);
	for (int i=0; i<row; i++)
   	for (int j=0; j<col; j++)
      	B.mat[i][j] = x*mat[i][j];
	return B;
}



// finds the inverse of the matrix
Matrix operator ~(Matrix A)
{
	double detA=A.determinant();
	Matrix temp;
	int i, j;
	if (detA == 0 || A.row!=A.col)
	{
		cout <<"The matrix has no inverse!\n";
		return temp;
	}
   temp.changeOrder(A.row, A.col);
   int n=A.row;

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
		{
			double t;
         Matrix aux(n, n);
			aux=A.part(i, j);
			t=aux.determinant();
			temp.mat[j][i] = t/detA;
			if ((i+j)%2)
				temp.mat[j][i] = -temp.mat[j][i];
	}
   return temp;
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


void main()
{
   int a;
   /*Matrix A(3, 3);
   cin >>A;
   cout <<A;
   cout <<endl;
   cout <<"Determinant = " <<A.determinant() <<endl;
   A=-A.scalarMult(3);*/
   Matrix C(3, 3, "2 3 2 5.7 -2 -5.4 5 3 3 4");
   cout <<endl <<C;

   cin >>a;
}

