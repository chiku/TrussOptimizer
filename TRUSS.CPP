// Fem analysis of a 2D truss

#include <iostream>
#include <strstream>
#include <conio.h>
#include <math.h>
#include "matrix.cpp"

const int MAX = 20;
#define EA 1E10


struct Node
{
	double x;
	double y;	// (x, y) are the coordinates

	double u;
	double v;	// (u, v) are the displacements
   char dispx;	// y means that displacement in x direction at that point is known
	char dispy;	// y means that displacement in y direction at that point is known

	double Fx;
	double Fy;	// (Fx, Fy) are the forces
   char forx;  // y means that force in x direction at that point is known
   char fory;	// y means that force in y direction at that point is known

	Node();		// constructor
	void getData();
};

Node::Node()
{
	x=y=u=v=Fx=Fy=0;
	forx=fory='n';
	dispx=dispy='n';
}

void Node::getData()
{
	cout <<"Enter the coordinates of the node(x, y): ";
	cin >>x >>y;

	cout <<"Is the force in x direction at  the node known(y/n)?";
	cin >>forx;
	if (forx == 'y' || forx == 'Y')
	{
		cout <<"Enter the x force of the node Fx: ";
		cin >>Fx;
	}
	cout <<"Is the force in y direction at  the node known(y/n)?";
   cin >>fory;
   if (fory == 'y' || fory == 'Y')
   {
	cout <<"Enter the y force of the node Fy: ";
	  cin >>Fy;
   }

	cout <<"Is the displacement in x direction at  the node known(y/n)?";
	cin >>dispx;
	if (dispx == 'y' || dispx == 'Y')
	{
		cout <<"Enter the x displacement of the node u: ";
		cin >>u;
	}
	cout <<"Is the displacement in y direction at  the node known(y/n)?";
   cin >>dispy;
   if (dispy == 'y' || dispy == 'Y')
   {
   	cout <<"Enter the y displacement of the node v: ";
   	cin >>v;
   }
}

class Truss
{
	private:
	  Node N[MAX];
     int total_nodes;
	  int total_members;

		Matrix kglobal;
		Matrix klocal[MAX];
		char connectivity[MAX][MAX];

   public:
		Truss();	// constructor accepts all input data
		void findKLocal();
      void findKGlobal();

      void printKLocal();
      void printKGlobal();
};

Truss::Truss()
{
   int i, j;
   total_members = 0; // the actual value is found in findKLocal()
	cout <<"Enter the number of nodes: ";
   cin >>total_nodes;

   changeOrder(kglobal, 2*total_nodes, 2*total_nodes);
   for (i=0; i<2*total_nodes; i++)
   	for (j=9; j<2*total_nodes; j++)
      	kglobal.mat[i][j] = 0;

   for (i=0; i<total_nodes; i++)
   {
   	cout <<"\n\tNODE " <<i <<endl;
      N[i].getData();
      connectivity[i][i] = 'n';  // a node is always connected to itself
   }

   for (i=0; i<total_nodes; i++)
   	for (j=i+1; j<total_nodes; j++)
      {
      	cout <<"Is node " <<i <<" connected to node " <<j <<" (y/n)?";
         cin >>connectivity[i][j];
         connectivity[j][i] = connectivity[i][j];
      }
}


// finds all the local k matrices
void Truss::findKLocal()
{
	total_members=0;
   int i, j;
   double Cx, Cy, len;

   for(i=0; i<total_nodes; i++)
	{
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
            total_members++;
         }
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
   //kglobal = scalarMult(kglobal, 0.5);
}



// prints all the k local matrices
void Truss::printKLocal()
{
	int i;
   cout <<"The local k matrices are\n";
   for (i=0; i<total_members; i++)
   {
   	cout <<"\n\n\tk" <<i <<"/EA=";
      cout <<klocal[i];
   }
}

// prints all the k local matrices
void Truss::printKGlobal()
{
	int i;
   cout <<"\n\nThe global k matrix is\n";
   cout <<"\n\n\tk" <<i <<"/EA=";
   cout <<kglobal;
}


void main()
{

   cout.precision(3);
	Truss T;
   T.findKLocal();
   T.printKLocal();
   T.findKGlobal();
   T.printKGlobal();
}


// finds k global
/* The members
Truss::findKGlobal()
{






void main(void)
{
	clrscr();
	int n,i,j,k=1,chk,memno=0,gd[4],p,q;
	// n = number of nodes
	cout<<"\nEnter the no. of nodes: ";
	cin>>n;
	static double kglobal[27][27],klocal[27][4][4],connectmat[15][15];;
	double dircosx,dircosy,len;
	points nodes[15];

	for(i=0;i<n;i++)
	{
		cout <<"\nNode " <<(i+1) <<endl;
		nodes[i]=getdata();
		if(i==0 || i==n-1)
		{
			nodes[i].displacex=0;
			nodes[i].displacey=0;
		}
		else
		{
			nodes[i].displacex=k++;
			nodes[i].displacey=k++;
		}
	}

	cout<<"\nEnter 1 if the following nodes are connected else 0.\n";

	for(i=0;i<n-1;i++)
	{
		for(j=i+1;j<n;j++)
		{
		cout<<"\nIs node "<<i+1<<" connected to node "<<j+1<<"?";
		cin>>chk;
			if(chk==1)
			{
				len=pow((nodes[i].nodex-nodes[j].nodex),2.0)+pow((nodes[i].nodey-nodes[j].nodey),2.0);
				connectmat[i][j]=connectmat[j][i]=sqrt(len);
			}
		}
	}

	// Finding the local stiffness matrices
	cout <<"\n\nThe local stiffness matrices are\n\n";
	memno=0;
	for(i=0;i<n-1;i++)
	{
		for(j=i+1;j<n;j++)
		{
			if(connectmat[i][j]>0)
			{
				dircosx=(nodes[j].nodex-nodes[i].nodex)/connectmat[i][j];
				dircosy=(nodes[j].nodey-nodes[i].nodey)/connectmat[i][j];
				klocal[memno][0][0]=klocal[memno][2][2]=pow(dircosx,2)/connectmat[i][j];
				klocal[memno][1][1]=klocal[memno][3][3]=pow(dircosy,2)/connectmat[i][j];
				klocal[memno][0][1]=klocal[memno][1][0]=klocal[memno][3][2]=klocal[memno][2][3]=dircosx*dircosy/connectmat[i][j];
				klocal[memno][3][0]=klocal[memno][0][3]=klocal[memno][1][2]=klocal[memno][2][1]=-dircosx*dircosy/connectmat[i][j];
				klocal[memno][0][2]=klocal[memno][2][0]=-pow(dircosx,2)/connectmat[i][j];
				klocal[memno][3][1]=klocal[memno][1][3]=-pow(dircosy,2)/connectmat[i][j];

				for(p=0;p<4;p++)
				{
					cout<<"\n";
					for(q=0;q<4;q++)
					{
						cout<<'\t'<<klocal[memno][p][q];
					}
				}
				cout<<"\n\n";
				getch();

				gd[0]=nodes[i].displacex;gd[1]=nodes[i].displacey;
				gd[2]=nodes[j].displacex;gd[3]=nodes[j].displacey;

				for(p=0;p<4;p++)
					cout<<gd[p]<<"   \n";

				for(p=0;p<4;p++)
				{
					for(q=0;q<4;q++)
					{
						if(gd[p]!=0 && gd[q]!=0)
							kglobal[gd[p]-1][gd[q]-1]+=klocal[memno][p][q];
					}
				}

			memno++;
			}
		}
	}

	// Displaying global matrix
	clrscr();
	cout <<" \n\nGlobal stiffness matrix is \n";
	for(i=0;i<k-1;i++)
	{
	cout<<"\n";
		for(j=0;j<k-1;j++)
		{
		kglobal[i][j]*=EA;
		cout<<'\t'<<kglobal[i][j];
		}
	}
	static double kglobalinv[27][27];

	cout <<endl<<endl;
	invert(kglobal, kglobalinv, k-1);
	cout.precision(4);

	// The forces
	for(i=1;i<n-1;i++)
	{
	nodes[i].displacex=0;
	nodes[i].displacey=0;
		for(p=0;p<(k-1)/2;p++)
		{
			nodes[i].displacex+=kglobalinv[2*i-2][2*p]*nodes[p+1].forcex+kglobalinv[2*i-2][2*p+1]*nodes[p+1].forcey;
			nodes[i].displacey+=kglobalinv[2*i-1][2*p]*nodes[p+1].forcex+kglobalinv[2*i-1][2*p+1]*nodes[p+1].forcey;
		}
	}

	cout<<"\n\nThe displacements are\n";

	for(p=1;p<n-1;p++)
	{
		cout<<"\n"<<nodes[p].displacex;
		cout<<"\n"<<nodes[p].displacey;
	}

	getch();
}              */