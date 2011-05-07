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
	void getData();
};

// accepts the coordinated of the node
void Node::getData()
{
	cout <<"Enter the coordinates of the node(x, y): ";
	cin >>x >>y;

}

class Truss
{
	private:
		Node N[MAX];
      Matrix force, displacement;
		int total_nodes;
		int total_members;

		Matrix kglobal;
		Matrix klocal[MAX];

		char connectivity[MAX][MAX];
      char knowledgeF[MAX], knowledgeu[MAX];
      Matrix k11, k12, k21, k22;
      Matrix Fknown, Funknown;
      Matrix uknown, uunknown;

   public:
		Truss();	// constructor accepts all input data
		void findKLocal();
      void findKGlobal();
      void condense();

      void printKLocal();
      void printKGlobal();
      void printCondensed();
};

// implementing 2D Truss class
Truss::Truss()
{
   int i, j;
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
}

/* static condensation scheme to find k11, k12, k21 and k22
	as well as Fknown, Funknown, uknownand uknown*/

void Truss::condense()
{
	int i, j;
   int k11i=0, k11j=0, k12i=0, k12j=0,  k21i=0, k21j=0,  k22i=0, k22j=0;
   for (i=0; i<2*total_nodes; i++)
   {
   	for (j=0; j<2*total_nodes; j++)
      {
        	if (knowledgeF[i]=='y')
         {
           	if (knowledgeu[j]=='y')
				{
               // force known + displacement known -> K11
            	k11.mat[k11i][k11j] = kglobal.mat[i][j];
               /*Fknown.mat[k11i][0] = force.mat[i][0];
               uknown.mat[k11j][0] = displacement.mat[j][0];*/
               k11j++;
            }
            else
            {
            	// force known + displacement unknown -> K12
            	k12.mat[k12i][k12j] = kglobal.mat[i][j];
               //Fknown.mat[k12i][0] = force.mat[i][0];
               k12j++;
            }
            k11i++;
            k12i++;
         }
         else
         {
            if ( knowledgeu[j]=='y')
				{
            	// force unknown + displacement known -> K21
            	k21.mat[k21i][k21j] = kglobal.mat[i][j];
               //uknown.mat[k21j][0] = displacement.mat[j][0];
               k21j++;
            }
            else
            {
            	// force unknown + displacement known -> K22
            	k22.mat[k22i][k22j] = kglobal.mat[i][j];
               k22j++;
            }
            k21i++;
            k22i++;
         }
      }
   }
   changeOrder(k11, k11i, k11j);
   changeOrder(k21, k21i, k21j);
   changeOrder(k12, k12i, k12j);
   changeOrder(k22, k22i, k22j);

   changeOrder(Fknown, k11j, 1);
   changeOrder(Funknown, k21j, 1);

   changeOrder(uknown, k11i, 1);
   changeOrder(uunknown, k21i, 1);
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
inline void Truss::printKGlobal()
{
	cout <<"\n\nThe uncondensed global k matrix is\n";
   cout <<"\n\n\tk/EA=";
   cout <<kglobal;
}


// prints
inline void Truss::printCondensed()
{
	cout <<"\n\n\n\tThe condensed matrices are\n";
   cout <<"\n\n\t\t\tK11/EA=" << k11;
   cout <<"\n\n\t\t\tK12/EA=" << k12;
   cout <<"\n\n\t\t\tK21/EA=" << k21;
   cout <<"\n\n\t\t\tK22/EA=" << k22;
   /*
   cout <<"\n\n\\t\t\tFknown=" <<Fknown;
   cout <<"\n\n\t\t\tUknown=" <<uknown;*/
}

void main()
{
   cout.precision(3);
	Truss T;
   T.findKLocal();
   T.findKGlobal();
   T.condense();

   T.printKLocal();
   T.printKGlobal();
   T.printCondensed();


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