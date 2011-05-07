// Fem analysis of a 2D truss

#include <iostream>
#include <strstream>
#include <conio.h>
#include <math.h>
#include "matrix.cpp"

const int MAX = 10;

struct Node
{
	double x;
	double y;	// (x, y) are the coordinates
	void getData();
};

// accepts the coordinated of the node
void Node::getData()
{
	cout <<"Enter the coordinates of the node(x, y): ";
	cin >>x >>y;

}


// implementing 2D Truss class
class Truss
{
	private:
		int total_nodes, total_members;
      Node *N;
      Matrix force, displacement;
      double *area;
      Matrix kglobal, kglobalcond, *klocal;
		Matrix k11, k12, k21, k22;
      Matrix Fknown, Funknown, uknown, uunknown;
      Matrix *locforce;
      double E;

		char (*connectivity)[MAX];
      char *knowledgeF, *knowledgeu;

   public:
		Truss();	// constructor accepts all input data
      ~Truss();
		void findKLocal();
      void findKGlobal();
      void condense();
      void solve();

      void printMatrices();
};

// constructor accepts the input data
Truss::Truss()
{
   int i, j;
   N = new Node[MAX];
   klocal = new Matrix[MAX];
   connectivity = new char[MAX][MAX];
   knowledgeF = new char[MAX];
   knowledgeu = new char[MAX];
   locforce = new Matrix[MAX];
   area = new double[MAX];

   total_members = 0; // the actual value is found in findKLocal()
	cout <<"Enter the number of nodes: ";
   cin >>total_nodes;

   changeOrder(kglobal, 2*total_nodes, 2*total_nodes);
   changeOrder(force, 2*total_nodes, 1);
   changeOrder(displacement, 2*total_nodes, 1);

   for (i=0; i<2*total_nodes; i++)
   	for (j=9; j<2*total_nodes; j++)
      	kglobal.mat[i][j] = 0;

   for (i=0; i<total_nodes; i++)
   {
   	cout <<"\n\tNODE " <<i <<endl;
      N[i].getData();
      connectivity[i][i] = 'n';  // avoids node being connected to itself

      // input forces
      cout <<"Is the force in x-direction known(y/n)?";
      cin >>knowledgeF[i*2];
      if (knowledgeF[i*2]=='y')
      	cin >>force.mat[2*i][0];

      cout <<"Is the force in y-direction known(y/n)?";
      cin >>knowledgeF[i*2+1];
      if (knowledgeF[i*2+1]=='y')
      	cin >>force.mat[2*i+1][0];

      // input displacement
      cout <<"Is the displacement in x-direction known(y/n)?";
      cin >>knowledgeu[i*2];
      if (knowledgeu[i*2]=='y')
      	cin >>displacement.mat[2*i][0];

      cout <<"Is the displacement in y-direction known(y/n)?";
      cin >>knowledgeu[i*2+1];
      if (knowledgeu[i*2+1]=='y')
      	cin >>displacement.mat[2*i+1][0];
   }

	int a=0;
   for (i=0; i<total_nodes; i++)
   	for (j=i+1; j<total_nodes; j++)
      {
      	cout <<"Is node " <<i <<" connected to node " <<j <<" (y/n)?";
         cin >>connectivity[i][j];
         connectivity[j][i] = connectivity[i][j];
         cout <<"Enter the area of the node: ";
         cin >>area[a++];
      }
   cout <<"Enter the value of E: ";
   cin >>E;
}


Truss::~Truss()
{
	delete[] klocal;
   delete[] N;
   delete[] klocal;
   delete[]connectivity;
   delete[]knowledgeF;
   delete[]knowledgeu;
   delete[] locforce;
   delete[] area;
}


// finds all the local k matrices
void Truss::findKLocal()
{
	total_members=0;
   int i, j;
   double Cx, Cy, len;

   for(i=0; i<total_nodes; i++)
   	for(j=i+1; j<total_nodes; j++)
      {
      	if(connectivity[i][j] == 'y' || connectivity[i][j] == 'Y')
			{
				len = sqrt( pow((N[i].x-N[j].x),2.0) + pow((N[i].y-N[j].y),2.0) );
            Cx = (N[j].x-N[i].x) / len;
				Cy = (N[j].y-N[i].y) / len;

            changeOrder(klocal[total_members], 4, 4);
				klocal[total_members].mat[0][0]=klocal[total_members].mat[2][2]=Cx*Cx/len;
				klocal[total_members].mat[1][1]=klocal[total_members].mat[3][3]=Cy*Cy/len;
            klocal[total_members].mat[0][2]=klocal[total_members].mat[2][0]=-Cx*Cx/len;
				klocal[total_members].mat[3][1]=klocal[total_members].mat[1][3]=-Cy*Cy/len;
            klocal[total_members].mat[0][1]=klocal[total_members].mat[1][0]=
            				klocal[total_members].mat[3][2]=klocal[total_members].mat[2][3]=Cx*Cy/len;
				klocal[total_members].mat[3][0]=klocal[total_members].mat[0][3]=
            				klocal[total_members].mat[1][2]=klocal[total_members].mat[2][1]=-Cx*Cy/len;
            klocal[total_members] = scalarMult(klocal[total_members], E*area[total_members]);
            total_members++;
         }
      }
}

// finds the uncondensed k global matrix
void Truss::findKGlobal()
{
	int i, j, mem_no=0;
   for (i=0; i<total_nodes-1; i++)
   	for (j=i; j<total_nodes; j++)
      	if (connectivity[i][j] == 'Y' || connectivity[i][j] == 'y' && i!=j)
         {
				kglobal.mat[i*2][i*2]     += klocal[mem_no].mat[0][0];
				kglobal.mat[i*2][i*2+1]   += klocal[mem_no].mat[0][1];
				kglobal.mat[i*2][j*2]     += klocal[mem_no].mat[0][2];
				kglobal.mat[i*2][j*2+1]   += klocal[mem_no].mat[0][3];

				kglobal.mat[i*2+1][i*2]   += klocal[mem_no].mat[1][0];
				kglobal.mat[i*2+1][i*2+1] += klocal[mem_no].mat[1][1];
				kglobal.mat[i*2+1][j*2]   += klocal[mem_no].mat[1][2];
				kglobal.mat[i*2+1][j*2+1] += klocal[mem_no].mat[1][3];

				kglobal.mat[j*2][i*2]     += klocal[mem_no].mat[2][0];
				kglobal.mat[j*2][i*2+1]   += klocal[mem_no].mat[2][1];
				kglobal.mat[j*2][j*2]     += klocal[mem_no].mat[2][2];
				kglobal.mat[j*2][j*2+1]   += klocal[mem_no].mat[2][3];

				kglobal.mat[j*2+1][i*2]   += klocal[mem_no].mat[3][0];
				kglobal.mat[j*2+1][i*2+1] += klocal[mem_no].mat[3][1];
				kglobal.mat[j*2+1][j*2]   += klocal[mem_no].mat[3][2];
				kglobal.mat[j*2+1][j*2+1] += klocal[mem_no].mat[3][3];
            mem_no++;
         }
}

/* static condensation scheme to find k11, k12, k21 and k22
	as well as Fknown, Funknown, uknownand uknown*/
void Truss::condense()
{
	int i, j;
   int count_fknown=0, count_funknown=0, count_uknown=0, count_uunknown=0;
   Matrix temp;    // holds the force condensation
   changeOrder(temp, total_nodes*2, total_nodes*2);
   changeOrder(kglobalcond, total_nodes*2, total_nodes*2);

   // condensing for forces
   // known forces at top
   for (i=0; i<2*total_nodes; i++)
     	if (knowledgeF[i]=='y')
		{
      	Fknown.mat[count_fknown][0]=force.mat[i][0];
         for (j=0; j<2*total_nodes; j++)
         	temp.mat[count_fknown][j]=kglobal.mat[i][j];
         count_fknown++;
	   }

   // unknown forces at bottom
   for (i=0; i<2*total_nodes; i++)
     	if (knowledgeF[i]=='n')
		{
      	Funknown.mat[count_funknown][0]=force.mat[i][0];
         for (j=0; j<2*total_nodes; j++)
         	temp.mat[count_fknown+count_funknown][j]=kglobal.mat[i][j];
         count_funknown++;
      }

   //condensing for displacements
   // known displacements at the top
   for (i=0; i<2*total_nodes; i++)
     	if (knowledgeu[i]=='y')
		{
      	uknown.mat[count_uknown][0]=displacement.mat[i][0];
         for (j=0; j<2*total_nodes; j++)
         	kglobalcond.mat[j][count_uknown]=temp.mat[j][i];
         count_uknown++;
      }
   // unknown displacements at the bottom
   for (i=0; i<2*total_nodes; i++)
     	if (knowledgeu[i]=='n')
		{
      	uknown.mat[count_uunknown][0]=displacement.mat[i][0];
         for (j=0; j<2*total_nodes; j++)
         	kglobalcond.mat[j][count_uknown+count_uunknown]=temp.mat[j][i];
         count_uunknown++;
      }

   changeOrder(Fknown, count_fknown, 1);
   changeOrder(Funknown, count_funknown, 1);
   changeOrder(uknown, count_uknown, 1);
   changeOrder(uknown, count_uunknown, 1);

   // create the matrices K11, K12, K21 and K22
   // f known, u known
   for (i=0; i<count_fknown; i++)
		for (j=0; j<count_uknown; j++)
        	k11.mat[i][j]=kglobalcond.mat[i][j];
   // f known, u unknown
   for (i=0; i<count_fknown; i++)
		for (j=0; j<count_uunknown; j++)
        	k12.mat[i][j]=kglobalcond.mat[i][count_uknown+j];
   // f unknown, u known
   for (i=0; i<count_funknown; i++)
		for (j=0; j<count_uknown; j++)
        	k21.mat[i][j]=kglobalcond.mat[count_fknown+i][j];
   // f unknown, u unknown
   for (i=0; i<count_funknown; i++)
		for (j=0; j<count_uunknown; j++)
        	k22.mat[i][j]=kglobalcond.mat[count_fknown+i][count_uknown+j];
   changeOrder(k11, count_fknown, count_uknown);
   changeOrder(k12, count_fknown, count_uunknown);
   changeOrder(k21, count_funknown, count_uknown);
   changeOrder(k22, count_funknown, count_uunknown);
}

void Truss::solve()
{
	uunknown = (~k12) * (Fknown - k11*uknown);
   Funknown = k21*uknown + k22*uunknown;

	int i, j;
   Matrix ulocal, uglobal;
   changeOrder(uglobal, 2*total_nodes, 1);
   changeOrder(ulocal, 4, 1);

   // find global displacement matrix
   int a=0, b=0;
   for (i=0; i<2*total_nodes; i++)
   {
      if (knowledgeu[i] == 'y')
   		uglobal.mat[i][0] = uknown.mat[a++][0];
      else
      	uglobal.mat[i][0] = uunknown.mat[b++][0];
   }

   a=0;
   for (i=0; i<total_nodes-1; i++)
      // find all local displacement matrices
   	for (j=i+1; j<total_nodes; j++)
      {
      	if (connectivity[i][j] == 'y')
         {
         	ulocal.mat[0][0] = uglobal.mat[2*i][0];
         	ulocal.mat[1][0] = uglobal.mat[2*i+1][0];
         	ulocal.mat[2][0] = uglobal.mat[2*j][0];
         	ulocal.mat[3][0] = uglobal.mat[2*j+1][0];

            locforce[a] = klocal[a]*ulocal;
            a++;
            changeOrder(locforce[a], 4, 1);
         }
      }
}

// prints all the matrices
void Truss::printMatrices()
{
	int i;
   cout <<"\n\nThe local k matrices are";
   for (i=0; i<total_members; i++)
   	cout <<"\n\n\tk" <<i << klocal[i];

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

}

void main()
{
   cout.precision(3);
	Truss T;
   T.findKLocal();
   T.findKGlobal();
   T.condense();
   T.solve();
   T.printMatrices();
}



