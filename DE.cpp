// Optimization of truss using DE

#include <iostream>
#include <fstream>
#include <conio.h>
#include <cmath>
#include "truss.cpp"

#ifndef __DE_H__
#define __DE_H__




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
		double avg_fitness, best_fitness;
		int best_fitness_loc;
		char file[]; // name of the data file
		
		double MIN_AREA, MAX_AREA, MAX_STRESS; // constraints
		int POPULATION;  // DE parameter
		double CR, F;	 // DE parameters
		double PENALTY;  // constant penalty

	protected:
		void findFitness();
		void findFitnessLast();		

	public:
		DE();
		void getData();
		void evolution();
};

// Constructor
DE::DE()
{
	T = new Truss[POPULATION+1];
	fitness = new double[POPULATION+1];
	strcpy(file, "truss.dat");
	for (int i=0; i<POPULATION+1; i++)
	{
		T[i].getData();
	
		// Random initilization of area
		double *a = new double[POPULATION+1];
		for (i=0; i<T[i].members(); i++)
			a[i] = randomBetween(MAX_AREA, MIN_AREA);
		T[i].modifyArea(a);
		delete[] a;
	}
	findFitness();
}


// Get the DE parameters
void DE::getData()
{
	ifstream File("de.dat");
	File >>MIN_AREA >>MAX_AREA >>MAX_STRESS;
	File >>PENALTY;
	File.close();
}

// Find the fitnesses inclusive of penalty
// Also finds the best and average fitness
void DE::findFitness()
{
	int i, j;
	for (i=0; i<POPULATION; i++)
	{
		// Solve and obtain the forces in the members
		T[i].findKLocal(); T[i].findKGlobal(); T[i].condense(); T[i].solve();
	
		fitness[i] = T[i].area[0];
		for (j=1; j<=T[i].members(); j++)
		{
			fitness[i] += T[i].area[j];
			double t = T[i].locforce[j].mat[0][0]*T[i].locforce[j].mat[0][0] + 
				T[i].locforce[j].mat[1][1]*T[i].locforce[j].mat[1][1];
			if (sqrt(t) / T[i].area[j] >= MAX_STRESS)
				fitness[i] += PENALTY;
		}
	}
	
	best_fitness = fitness[0]; 
	best_fitness_loc = 0;
	avg_fitness = fitness[0]; 
	for (i=1; i<POPULATION; i++)
	{
		if (fitness[i] < fitness[best_fitness_loc])
		{
			best_fitness_loc = i;
			best_fitness = fitness[i];
		}
		avg_fitness += fitness[i];
	}
	avg_fitness /= POPULATION;
}
	
	
// fitness of the last member which is actually vtrial
void DE::findFitnessLast()
{
	int i, j;
	i = POPULATION;
	// Solve and obtain the forces in the members
	T[i].findKLocal(); T[i].findKGlobal(); T[i].condense(); T[i].solve();
	
	fitness[i] = T[i].area[0];
	for (j=1; j<=T[i].members(); j++)
	{
		fitness[i] += T[i].area[j];
		double t = T[POPULATION].locforce[j].mat[0][0]*T[POPULATION].locforce[j].mat[0][0] + 
				T[POPULATION].locforce[j].mat[1][1]*T[POPULATION].locforce[j].mat[1][1];
		if (sqrt(t) / T[i].area[j] >= MAX_STRESS)
			fitness[i] += PENALTY;
	}
}



// Main algorithm for evolution
void DE::evolution()
{
	int i, j, k;
	
	Truss temp;
	temp.getData();
	int GENERATION = 0;
	avg_fitness = 1e50; best_fitness = 1e49;
	while (avg_fitness - best_fitness >= 0.0001 && GENERATION < 10000)
	{
		findFitness();
		GENERATION++;

		cout <<"Generation number: " <<GENERATION <<"\tBest fitness: " <<best_fitness
			<<"\tAvg. fitness: " <<avg_fitness <<endl;
		int r1 = randomBetween(0, POPULATION);
		int r2 = randomBetween(0, POPULATION);
		int r3 = randomBetween(0, POPULATION);
		int r4 = randomBetween(0, POPULATION);

		// Diffenential mutation		
		for (j=0; j<=T[i].members(); j++)
		{
			temp.area[j] = T[best_fitness_loc].area[j] + F * (T[r1].area[j] - T[r2].area[j]);
			if (temp.area[j] <= MIN_AREA) temp.area[j] = MIN_AREA;
			if (temp.area[j] >= MAX_AREA) temp.area[j] = MAX_AREA;
		}
		
		// Crossing over
		for (j=0; j<=T[0].members(); j++)
		{
			if ( randomBetween(0.0, 1.0) < CR )
				T[POPULATION].area[j] = temp.area[j];
			else
				T[POPULATION].area[j] = T[r3].area[j];
		}
		
		// Recombination
		findFitnessLast();
		if (fitness[POPULATION] < fitness[r4])
			for (j=0; j<T[0].members(); j++)
				T[r3].area[j] = T[POPULATION].area[j];
	}
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


