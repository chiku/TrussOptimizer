// Optimization of truss using DE

#include <iostream>
#include <fstream>
#include <conio.h>
#include <cmath>
#include "truss.cpp"

#ifndef __DE_H__
#define __DE_H__

using namespace std;

inline int randomBetween(int a, int b) // [a, b)
{
	return ( rand() * (b-a) /(RAND_MAX) + a);
}

inline double randomBetween(double a, double b) // [a, b)
{
	return ( rand() * (b-a) /(0.0+RAND_MAX) + a);
}


class DE
{
	private:
		Truss *T;
		double *fitness;
		double avg_fitness, best_fitness, total_fitness;
		int best_fitness_loc;
		
		double MIN_AREA, MAX_AREA;  // geometric constraints
		double MAX_STRESS, MAX_DISP; // behaviour constraints
		int POPULATION;  // DE parameter
		double F_MIN, F_MAX, CR_MIN, CR_MAX;	 // DE parameters
		double PENALTY;  // constant penalty

	protected:
		void findFitness();
		void findBestFitness(); // the fittest solution
		void findFitnessTrial(); // fitness of the trial solution

	public:
		DE();
		~DE();
		void getData();
		void evolution();
		void printResult()
		{
			int i;
			cout <<"\n\nThe stress in members are\n";
			for (i=0; i<T[0].members(); i++)
			{
				double t =  T[best_fitness_loc].locforce[i].mat[0][0] * T[best_fitness_loc].locforce[i].mat[0][0] + 
							T[best_fitness_loc].locforce[i].mat[1][0] * T[best_fitness_loc].locforce[i].mat[1][0];
				cout <<"\t" <<i <<">  " <<sqrt(t)/T[best_fitness_loc].area[i] <<endl;
			}		
			cout <<"\n\nThe nodal displacements are\n";
			for (i=0; i<T[0].uglobal.row; i++)
				cout <<"\t" <<i <<"> " <<T[best_fitness_loc].uglobal.mat[i][0] <<endl;
			cout <<"\n\nThe member area are\n";
			for (i=0; i<T[0].members(); i++)
				cout <<"Member no "<< i <<" " <<T[best_fitness_loc].area[i] <<endl;
		}
};

// Constructor
DE::DE()
{
	getData();
	
	T = new Truss[POPULATION+1];
	fitness = new double[POPULATION+1];
	for (int i=0; i<POPULATION+1; i++)
	{
		T[i].getData();
		// Random initilization of area
		for (int j=0; j<T[i].members(); j++)
			T[i].area[j] = randomBetween(MIN_AREA, MAX_AREA);
	}
}


DE::~DE()
{
	delete[] T;
	delete[] fitness;
}


// Get the DE parameters
void DE::getData()
{
	ifstream File("de.dat");
	File >>MIN_AREA >>MAX_AREA;
	File >>MAX_STRESS >>MAX_DISP;
	File >>PENALTY;
	File >>POPULATION >>F_MIN >>F_MAX >>CR_MIN >>CR_MAX;
	File.close();
}

// Find the fitnesses inclusive of penalty
void DE::findFitness()
{
	int i, j;
	for (i=0; i<POPULATION; i++)
	{
		// Solve and obtain the forces in the members
		T[i].findKLocal(); 
		T[i].findKGlobal(); 
		T[i].condense(); 
		T[i].solve();

		fitness[i] = 0;
		for (j=0; j<T[0].members(); j++)
		{
			fitness[i] += T[i].area[j];
			// stress penalty
			double t = sqrt( T[i].locforce[j].mat[0][0]*T[i].locforce[j].mat[0][0] + 
				T[i].locforce[j].mat[1][0]*T[i].locforce[j].mat[1][0] );
			if ( fabs(t / T[i].area[j]) >= MAX_STRESS)
				fitness[i] += PENALTY;
		}

		// displacement penalty		
		for (j=0; j<T[0].uglobal.row; j++)
			if ( fabs(T[i].uglobal.mat[j][0]) >= MAX_DISP )
				fitness[i] += PENALTY;
	}
}


// Find the fittest solution and the average fitness
void DE::findBestFitness()
{
	best_fitness = fitness[0]; 
	best_fitness_loc = 0;
	total_fitness = fitness[0]; 
	for (int i=1; i<POPULATION; i++)
	{
		if (fitness[i] < fitness[best_fitness_loc])
		{
			best_fitness_loc = i;
			best_fitness = fitness[i];
		}
		total_fitness += fitness[i];
	}
	avg_fitness = total_fitness/POPULATION;
}


// fitness of the trial solution
void DE::findFitnessTrial()
{
	int j;
	// Solve and obtain the forces in the members
	T[POPULATION].findKLocal();
	T[POPULATION].findKGlobal();
	T[POPULATION].condense();
	T[POPULATION].solve();
	
	fitness[POPULATION] = 0;
	for (j=0; j<T[POPULATION].members(); j++)
	{
		fitness[POPULATION] += T[POPULATION].area[j];
		// stress penalty
		double t = sqrt( T[POPULATION].locforce[j].mat[0][0]*T[POPULATION].locforce[j].mat[0][0] + 
			T[POPULATION].locforce[j].mat[1][0]*T[POPULATION].locforce[j].mat[1][0] );
		if ( fabs(t / T[POPULATION].area[j]) >= MAX_STRESS)
			fitness[POPULATION] += PENALTY;
	}
	// displacement penalty		
	for (j=0; j<T[POPULATION].uglobal.row; j++)
		if ( fabs(T[POPULATION].uglobal.mat[j][0]) >= MAX_DISP)
			fitness[POPULATION] += PENALTY;
}



// Main algorithm for evolution
void DE::evolution()
{
	int i, j, k;
	
	Truss temp;
	temp.getData();
	findFitness();
	findBestFitness();

	long int GENERATION = 0;
	do
	{
		GENERATION++;

		cout <<"Generation: " <<GENERATION <<"\tBest fit.: " <<best_fitness
			<<"\tAvg. fit.: " <<avg_fitness <<endl;
		int r1 = randomBetween(0, POPULATION);
		int r2 = randomBetween(0, POPULATION);
		int r3 = randomBetween(0, POPULATION);
		int r4 = randomBetween(0, POPULATION);
		int r5 = randomBetween(0, POPULATION);

		// Differential mutation		
		for (j=0; j<=T[0].members(); j++)
		{
			double F = randomBetween(F_MIN, F_MAX);
			temp.area[j] = T[best_fitness_loc].area[j] + F * (T[r2].area[j] - T[r3].area[j]);
			if (temp.area[j] <= MIN_AREA) temp.area[j] = MIN_AREA;
			if (temp.area[j] >= MAX_AREA) temp.area[j] = MAX_AREA;
		}

		// Recombination
		for (j=0; j<=T[0].members(); j++)
		{
			double CR = randomBetween(CR_MIN, CR_MAX);
			if ( randomBetween(0.0, 1.0) < CR )
				T[POPULATION].area[j] = temp.area[j];
			else
				T[POPULATION].area[j] = T[r4].area[j];
		}

		// Selection
		findFitnessTrial();
		if (fitness[POPULATION] < fitness[r5])
		{
			for (j=0; j<T[0].members(); j++)
				T[r5].area[j] = T[POPULATION].area[j];
			
			// update fitnesses
			total_fitness = total_fitness - fitness[r5] + fitness[POPULATION];
			avg_fitness = total_fitness / POPULATION;
			fitness[r5] = fitness[POPULATION];
			
			// check if we found the fittest solution
			if (fitness[r5] < best_fitness)
			{
				best_fitness = fitness[r5];
				best_fitness_loc = r5;
			}
		}
	}while (avg_fitness - best_fitness >= 0.0001 & GENERATION < 5000);
}

#endif

/*
int main()
{
	cout.precision(5);
	for (int i=0; i<10000; i++)
		cout <<randomBetween(1.0, 10.0) <<endl;	
	getch();
}
*/


