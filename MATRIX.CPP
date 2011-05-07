#include <iostream>
#include <strstream>

/* IMPLEMENTING MATRIX */

struct Matrix
{
  	double mat[10][10];
   int row, col;
};


// changes order of matrix
void changeOrder(Matrix& A, int r, int c)
{
	A.row=r;
   A.col=c;
}

// determines corresponding determinant
double determinant(Matrix A)
{
	Matrix aux;
	double d=0;

   int n;
   if (A.row!=A.col)
   {
   	cerr <<"\nOnly square matrices can have corresponding determinants.";
      return d;
   }
   n=A.row;

   if (n==1)
		return A.mat[0][0];

	if (n==2)
		return (A.mat[0][0]*A.mat[1][1] - A.mat[1][0]*A.mat[0][1]);

   changeOrder(aux, n-1, n-1);
	for (int a=0; a<n; a++)
	{
		int i, j;
		for (i=0; i<n; i++)
			for (j=0; j<a; j++)
				aux.mat[i][j] = A.mat[i+1][j];

		for (i=0; i<n; i++)
			for (j=a; j<n; j++)
				aux.mat[i][j] = A.mat[i+1][j+1];

		if (a%2)
			d -= A.mat[0][a]*determinant(aux);
		else
			d += A.mat[0][a]*determinant(aux);
	}
	return d;
}

// used in finding the inverse of matrix
Matrix part(Matrix A, int m, int n)
{
   Matrix B;
   changeOrder(B, A.row-1, A.col-1);
	int i,j;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			B.mat[i][j] = A.mat[i][j];
	for (i=0; i<m; i++)
		for (j=n; j<A.col; j++)
			B.mat[i][j] = A.mat[i][j+1];

	for (i=m; i<A.row; i++)
		for (j=0; j<n; j++)
			B.mat[i][j] = A.mat[i+1][j];
	for (i=m; i<A.row; i++)
		for (j=n; j<A.col; j++)
			B.mat[i][j] = A.mat[i+1][j+1];
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
   changeOrder(c, a.row, a.col);
   for (int i=0; i<a.row; i++)
   	for (int j=0; j<a.col; j++)
      	c.mat[i][j]=a.mat[i][j]+b.mat[i][j];
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
   changeOrder(C, A.row, A.col);
   for (int i=0; i<A.row; i++)
   	for (int j=0; j<A.col; j++)
      	C.mat[i][j]=A.mat[i][j]-B.mat[i][j];
   return C;
}

// matrix multiplication
Matrix operator *(Matrix A, Matrix B)
{
	Matrix C;
   if (A.col!=B.row)
   {
   	cerr <<"Matries incompatible for multiplication.\n";
      return C;
   }
   changeOrder(C, A.row, B.col);
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



// finds the inverse of the matrix
Matrix operator ~(Matrix A)
{
	double detA=determinant(A);
	Matrix temp;
	int i, j;
	if (detA == 0 || A.row!=A.col)
	{
		cout <<"The matrix has no inverse!\n";
		return temp;
	}
   changeOrder(temp, A.row, A.col);
   int n=A.row;

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
		{
			double t;
         Matrix aux;
         changeOrder(aux, n, n);
			aux=part(A, i, j);
			t=determinant(aux);
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

/*
void main()
{
   int a;
   Matrix C, D;
   changeOrder(C, 1, 4);
   changeOrder(D, 4, 1);
   C.mat[0][0]=1; C.mat[0][1]=2;
   C.mat[0][2]=2; C.mat[0][3]=3;

	D.mat[0][0]=3; D.mat[1][0]=2;
   D.mat[2][0]=4; D.mat[3][0]=3;

   Matrix E=C*D;
   cout <<endl <<E;
   cout <<endl <<C;
   cout <<endl <<D;

   cin >>a;
}
  */
