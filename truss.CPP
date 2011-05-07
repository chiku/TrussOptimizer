// Fem analysis of a simple truss

#include <iostream.h>
#include <conio.h>
#include <math.h>
#define EA 1E10

struct points
{
	double nodex;
	double nodey;
	double displacex;
	double displacey;
	double forcex;
	double forcey;
};

points getdata(void)
{
	points dummy;
	cout<<"Enter the coordinates of the node(x y): ";
	cin>>dummy.nodex>>dummy.nodey;
	cout<<"Enter the forces at the node(Fx Fy): ";
	cin>>dummy.forcex>>dummy.forcey;
	return dummy;
}

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
}