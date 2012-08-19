// Fem analysis of a 2D truss

#include <iostream>
#include <fstream>
#include <cmath>
#include "matrix.cpp"

#ifndef __TRUSS_H__
#define __TRUSS_H__

using namespace std;

const int MAX = 30;

struct Node
{
	double x;
	double y;	// (x, y) are the coordinates
};


// implementing 2D Truss class
class Truss
{
	private:
		int total_nodes, total_members;
		Node *N;
		Matrix force, displacement;
		Matrix kglobal, kglobalcond, *klocal;
		Matrix k11, k12, k21, k22;
		Matrix Fknown, Funknown, uknown, uunknown;
		double E;

		char (*connectivity)[MAX];
		char *knowledgeF, *knowledgeu;

   public:
		double *area;
		double *length;
		Matrix *locforce;
		Matrix uglobal;
				
		void getData();	// accepts all input data
		Truss();
		~Truss();
		void findKLocal();
		void findKGlobal();
		void condense();
		void solve();
		
		int members() { return total_members; };
		int nodes() { return total_nodes; };
		void modifyArea(double new_area[]);
		void printMatrices();
};


Truss::Truss()
{
	N = new Node[MAX];
	klocal = new Matrix[MAX];
	connectivity = new char[MAX][MAX];
	knowledgeF = new char[MAX];
	knowledgeu = new char[MAX];
	locforce = new Matrix[MAX];
	area = new double[MAX];
	length = new double[MAX];
}


// accepts the input data
void Truss::getData()
{
	int i, j;
	ifstream F("truss.dat");

	total_members = 0; // the actual value is found in findKLocal()
	// the number of nodes
	F >>total_nodes;

	changeOrder(kglobal, 2*total_nodes, 2*total_nodes);
	changeOrder(force, 2*total_nodes, 1);
	changeOrder(displacement, 2*total_nodes, 1);

	for (i=0; i<2*total_nodes; i++)
		for (j=9; j<2*total_nodes; j++)
			kglobal(i, j) = 0;

	for (i=0; i<total_nodes; i++)
	{
		// coordinates of the nodes
		F >>N[i].x >>N[i].y;
		connectivity[i][i] = 'n';  // avoids node being connected to itself

		// input forces
		// Is the force in x-direction known(y/n)?
		F >>knowledgeF[i*2];
		if (knowledgeF[i*2]=='y')
			F >>force(2*i, 0);

		// Is the force in y-direction known(y/n)?
		F >>knowledgeF[i*2+1];
		if (knowledgeF[i*2+1]=='y')
			F >>force(2*i+1, 0);

		// input displacement
		// Is the displacement in x-direction known(y/n)?
		F >>knowledgeu[i*2];
		if (knowledgeu[i*2]=='y')
			F >>displacement(2*i, 0);

		// Is the displacement in y-direction known(y/n)?
		F >>knowledgeu[i*2+1];
		if (knowledgeu[i*2+1]=='y')
			F >>displacement(2*i+1, 0);
	}

	int a=0;
	for (i=0; i<total_nodes; i++)
		for (j=i+1; j<total_nodes; j++)
		{
			// Is node i connected to node j (y/n)?;
			F >>connectivity[i][j];
			connectivity[j][i] = connectivity[i][j];
			// The area of the node (dummy value for optimization)
			if (connectivity[i][j] == 'y')
				F >>area[a++];
		}
	// The value of E
	F >>E;
	F.close();

	findKLocal(); // now the number of members present is known
}


Truss::~Truss()
{
	delete[] N;
	delete[] klocal;
	delete[] connectivity;
	delete[] knowledgeF;
	delete[] knowledgeu;
	delete[] locforce;
	delete[] area;
	delete[] length;
}


// finds all the local k matrices
void Truss::findKLocal()
{
	total_members=0;
	int i, j;
	double Cx, Cy, len;

	for(i=0; i<total_nodes; i++)
		for(j=i+1; j<total_nodes; j++)
			if(connectivity[i][j] == 'y' || connectivity[i][j] == 'Y')
			{
				length[total_members] = len = sqrt ( (N[i].x-N[j].x)*(N[i].x-N[j].x) + (N[i].y-N[j].y)*(N[i].y-N[j].y) );
				Cx = (N[j].x-N[i].x) / len;
				Cy = (N[j].y-N[i].y) / len;

				changeOrder(klocal[total_members], 4, 4);
				klocal[total_members](0, 0) = klocal[total_members](2, 2) = Cx*Cx/len;
				klocal[total_members](1, 1) = klocal[total_members](3, 3) = Cy*Cy/len;
				klocal[total_members](0, 2) = klocal[total_members](2, 0) = -Cx*Cx/len;
				klocal[total_members](3, 1) = klocal[total_members](1, 3) = -Cy*Cy/len;
				klocal[total_members](0, 1) = klocal[total_members](1, 0) =
						klocal[total_members](3, 2) = klocal[total_members](2, 3) = Cx*Cy/len;
				klocal[total_members](3, 0) = klocal[total_members](0, 3)=
						klocal[total_members](1, 2) = klocal[total_members](2, 1)= -Cx*Cy/len;
				klocal[total_members] = klocal[total_members].scalarMult(E*area[total_members]);
				total_members++;
			}
}

// finds the uncondensed k global matrix
void Truss::findKGlobal()
{
	int i, j, mem_no=0;
 	zeroMatrix(kglobal);
	for (i=0; i<total_nodes-1; i++)
		for (j=i; j<total_nodes; j++)
	 	if ( (connectivity[i][j] == 'Y' || connectivity[i][j] == 'y') && i!=j)
		 {
				kglobal(i*2  , i*2  ) += klocal[mem_no](0, 0);
				kglobal(i*2  , i*2+1) += klocal[mem_no](0, 1);
				kglobal(i*2  , j*2  ) += klocal[mem_no](0, 2);
				kglobal(i*2  , j*2+1) += klocal[mem_no](0, 3);

				kglobal(i*2+1, i*2  ) += klocal[mem_no](1, 0);
				kglobal(i*2+1, i*2+1) += klocal[mem_no](1, 1);
				kglobal(i*2+1, j*2  ) += klocal[mem_no](1, 2);
				kglobal(i*2+1, j*2+1) += klocal[mem_no](1, 3);

				kglobal(j*2  , i*2  ) += klocal[mem_no](2, 0);
				kglobal(j*2  , i*2+1) += klocal[mem_no](2, 1);
				kglobal(j*2  , j*2  ) += klocal[mem_no](2, 2);
				kglobal(j*2  , j*2+1) += klocal[mem_no](2, 3);

				kglobal(j*2+1, i*2  ) += klocal[mem_no](3, 0);
				kglobal(j*2+1, i*2+1) += klocal[mem_no](3, 1);
				kglobal(j*2+1, j*2  ) += klocal[mem_no](3, 2);
				kglobal(j*2+1, j*2+1) += klocal[mem_no](3, 3);
				
				mem_no++;
		 }
}

/* static condensation scheme to find k11, k12, k21 and k22
	as well as Fknown, Funknown, uknownand uknown*/
void Truss::condense()
{
	int i, j;
	int count_fknown=0, count_funknown=0, count_uknown=0, count_uunknown=0;
	Matrix temp;	 // holds the force condensation
	changeOrder(temp, total_nodes*2, total_nodes*2);
	changeOrder(kglobalcond, total_nodes*2, total_nodes*2);

	// condensing for forces
	// known forces at top
	for (i=0; i<2*total_nodes; i++)
		if (knowledgeF[i]=='y')
		{
			Fknown.mutateToInclude(force(i, 0), count_fknown, 0);
			for (j=0; j<2*total_nodes; j++)
				temp.mutateToInclude(kglobal(i, j), count_fknown, j);
			count_fknown++;
		}

	// unknown forces at bottom
	for (i=0; i<2*total_nodes; i++)
		if (knowledgeF[i]=='n')
		{
			Funknown.mutateToInclude(force(i, 0), count_funknown, 0);
			for (j=0; j<2*total_nodes; j++)
				temp.mutateToInclude(kglobal(i, j), count_fknown+count_funknown, j);
			count_funknown++;
		}

	//condensing for displacements
	// known displacements at the top
	for (i=0; i<2*total_nodes; i++)
		if (knowledgeu[i]=='y')
		{
			uknown.mutateToInclude(displacement(i, 0), count_uknown, 0);
			for (j=0; j<2*total_nodes; j++)
				kglobalcond.mutateToInclude(temp(j, i), j, count_uknown);
			count_uknown++;
		}
		
	// unknown displacements at the bottom
	for (i=0; i<2*total_nodes; i++)
		if (knowledgeu[i]=='n')
		{
			uunknown.mutateToInclude(displacement(i, 0), count_uunknown, 0);
			for (j=0; j<2*total_nodes; j++)
				kglobalcond.mutateToInclude(temp(j, i), j, count_uknown+count_uunknown);
			count_uunknown++;
		}

	changeOrder(Fknown, count_fknown, 1);
	changeOrder(Funknown, count_funknown, 1);
	changeOrder(uknown, count_uknown, 1);
	changeOrder(uunknown, count_uunknown, 1);

	changeOrder(k11, count_fknown, count_uknown);
	changeOrder(k12, count_fknown, count_uunknown);
	changeOrder(k21, count_funknown, count_uknown);
	changeOrder(k22, count_funknown, count_uunknown);

	// create the matrices K11, K12, K21 and K22
	// f known, u known
	for (i=0; i<count_fknown; i++)
		for (j=0; j<count_uknown; j++)
			k11(i, j) = kglobalcond(i, j);
	// f known, u unknown
	for (i=0; i<count_fknown; i++)
		for (j=0; j<count_uunknown; j++)
			k12(i, j) = kglobalcond(i, count_uknown+j);
	// f unknown, u known
	for (i=0; i<count_funknown; i++)
		for (j=0; j<count_uknown; j++)
			k21(i, j) = kglobalcond(count_fknown+i, j);
	// f unknown, u unknown
	for (i=0; i<count_funknown; i++)
		for (j=0; j<count_uunknown; j++)
			k22(i, j) = kglobalcond(count_fknown+i, count_uknown+j);
}


void Truss::solve()
{
	Matrix k12inv; k12inv = k12; k12inv.invert();
	uunknown = k12inv * (Fknown - k11*uknown);
	Funknown = k21*uknown + k22*uunknown;

	int i, j;
	Matrix ulocal;
	changeOrder(uglobal, 2*total_nodes, 1);
	changeOrder(ulocal, 4, 1);

	// find global displacement matrix
	int a=0, b=0;
	for (i=0; i<2*total_nodes; i++)
	{
		if (knowledgeu[i] == 'y')
			uglobal(i, 0) = uknown(a++, 0);
		else
			uglobal(i, 0) = uunknown(b++, 0);
	}

	a=0;
	for (i=0; i<total_nodes-1; i++)
	 // find all local displacement matrices
		for (j=i+1; j<total_nodes; j++)
	 		if (connectivity[i][j] == 'y')
			{
			 	ulocal(0, 0) = uglobal(2*i  , 0);
			 	ulocal(1, 0) = uglobal(2*i+1, 0);
			 	ulocal(2, 0) = uglobal(2*j  , 0);
		 		ulocal(3, 0) = uglobal(2*j+1, 0);

				locforce[a] = klocal[a]*ulocal;
				a++;
		 	}
}


// Modify the area of the truss elements
void Truss::modifyArea(double new_area[])
{
	for (int i=0; i<total_members; i++)
		area[i] = new_area[i];
}

// prints all the matrices
void Truss::printMatrices()
{
	int i;
	cout <<"\n\nThe local k matrices are";
	for (i=0; i<total_members; i++)
		cout <<"\n\n\tklocal" <<i << klocal[i];

	cout <<"\n\n\nThe uncondensed global k matrix is\n" <<kglobal;
	cout <<"\n\n\nThe condensed global k matrix is\n" <<kglobalcond;

	cout <<"\n\n\n\t\tK11=" << k11 <<"\n\n\t\tK12=" << k12
			<<"\n\n\t\tK21=" << k21 <<"\n\n\t\tK22=" << k22;

	cout <<"\n\nFknown=" <<Fknown;
	cout <<"\n\nUknown=" <<uknown;

	cout <<"\n\n\nUunknown=" <<uunknown;
	cout <<"\n\nFunknown=" <<Funknown;

	cout <<"\n\n\nThe forces in the members are";
	for (i=0; i<total_members; i++)
		cout <<"\n\nMember #" <<i <<locforce[i];
		
	cout <<"\n\nThe global displacement matrix is" <<uglobal;
}

#endif

/*
int main()
{
	Truss T;
	T.getData();
	T.findKLocal();
	T.findKGlobal();
	T.condense();
	T.solve();
	T.printMatrices();
}
*/
