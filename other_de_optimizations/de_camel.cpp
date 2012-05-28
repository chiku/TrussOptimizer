// Optimization of six camel humpback using DE

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>

#ifndef __DE_H__
#define __DE_H__

using namespace std;

inline int randomBetween(int a, int b) // [a, b)
{
	return ( a + rand() % (b+1-a) );
}

inline double randomBetween(double a, double b) // [a, b)
{
	return ( rand() * (b-a) /(0.0+RAND_MAX) + a);
}



class DE
{
	private:
		double *xvector, *yvector;;
		double *fitness;
		double avg_fitness, best_fitness;
		int best_fitness_loc;
		
		double MIN_X, MAX_X, MIN_Y, MAX_Y; // constraints
		int POPULATION;  // DE parameter
		double CR_MIN, CR_MAX, F_MIN, F_MAX;	 // DE parameters
		double PENALTY;  // constant penalty

	protected:
		void getData();
		void findFitness();
		void findFitnessLast();		

	public:
		DE();
		void evolution();
		void printResult()
		{
			cout <<"\n\nRESULTS" <<endl <<"\n---------------------------"
				<<"\nAverage fitness: " <<avg_fitness <<"\nBest fitness: " <<best_fitness
				<<"\n\nSolution" <<"\n\tx = " <<xvector[best_fitness_loc]
				<<"\n\ty = " <<yvector[best_fitness_loc]
				<<"\n\tFunction = " <<fitness[best_fitness_loc] <<"\n";
		}
};

// Constructor
DE::DE()
{
	getData();

	xvector = new double[POPULATION+1];
	yvector = new double[POPULATION+1];
	fitness = new double[POPULATION+1];

	for (int i=0; i<POPULATION; i++)
	{
		xvector[i] = randomBetween(MIN_X, MAX_X);
		yvector[i] = randomBetween(MIN_Y, MAX_Y);
	}		
}


// Get the DE parameters
void DE::getData()
{
	ifstream File("de_camel.dat");
	File >>MIN_X >>MAX_X >>MIN_Y >>MAX_Y;
	File >>POPULATION >>F_MIN >>F_MAX >>CR_MIN >>CR_MAX;
	File.close();
}

// Find the fitnesses inclusive of penalty
// Also finds the best and average fitness
void DE::findFitness()
{
	int i;
	for (i=0; i<POPULATION; i++)
	{
		double x=xvector[i], y=yvector[i];
		fitness[i] = ( 4.0 - 2.1*x*x + x*x*x*x/3.0 )*x*x + x*y + ( -4.0 + 4*y*y)*y*y;
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
inline void DE::findFitnessLast()
{
	double x=xvector[POPULATION], y=yvector[POPULATION];
	fitness[POPULATION] = ( 4.0 - 2.1*x*x + x*x*x*x/3.0 )*x*x + x*y + ( -4.0 + 4*y*y )*y*y;
}


// Main algorithm for evolution
void DE::evolution()
{
	int GENERATION = 0;
	double temp_x, temp_y;
	double F, CR;

	avg_fitness = 1e50; best_fitness = 1e49;

	while (avg_fitness - best_fitness >= 0.000000001 && GENERATION < 10000)
	{
		findFitness();
		GENERATION++;

		cout <<"Gen.: " <<GENERATION <<"\tBest fit.: " <<best_fitness
			<<"\tAvg. fit.: " <<avg_fitness << "\tbest (x,y) = " <<xvector[best_fitness_loc] 
			<<"\t" << yvector[best_fitness_loc] <<endl;
		int r1 = randomBetween(0, POPULATION);
		int r2 = randomBetween(0, POPULATION);
		int r3 = randomBetween(0, POPULATION);
		int r4 = randomBetween(0, POPULATION);

		// Differential mutation
		F = randomBetween(F_MIN, F_MAX);
		temp_x = xvector[best_fitness_loc] + F * (xvector[r1] - xvector[r2]);
		temp_y = yvector[best_fitness_loc] + F * (yvector[r1] - yvector[r2]);
		if (temp_x <= MIN_X) temp_x = MIN_X;
		if (temp_y <= MIN_Y) temp_y = MIN_Y;

		// Crossing over
		CR = randomBetween(CR_MIN, CR_MAX);
		if ( randomBetween(0.0, 1.0) < CR )
		{
			xvector[POPULATION] = temp_x;
			yvector[POPULATION] = temp_y;
		}
		else
		{
			xvector[POPULATION] = xvector[r3];
			yvector[POPULATION] = yvector[r3];
		}

		// Recombination
		findFitnessLast();
		if (fitness[POPULATION] < fitness[r4])
		{
			xvector[r4] = xvector[POPULATION];
			yvector[r4] = yvector[POPULATION];
		}
	}
}

#endif

int main()
{
	srand(time(0));
	cout.precision(8);
	DE de;
	de.evolution();
	de.printResult();
}
