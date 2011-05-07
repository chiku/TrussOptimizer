// Optimization of weight of spring using DE

#include <iostream>
#include <fstream>
#include <conio.h>
#include <cmath>

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
		double *x1vector, *x2vector, *x3vector;
		double *fitness;
		double *g1, *g2, *g3, *g4;
		double avg_fitness, best_fitness;
		int best_fitness_loc;
		
		double MIN_X1, MAX_X1, MIN_X2, MAX_X2, MIN_X3, MAX_X3; // constraints
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
				<<"\n\nSolution" <<"\n\tx1 = " <<x1vector[best_fitness_loc]
				<<"\n\tx2 = " <<x2vector[best_fitness_loc]
				<<"\n\tx3 = " <<x3vector[best_fitness_loc]
				<<"\n\tFunction = " <<fitness[best_fitness_loc]
				<<"\nConstraints: \n\t" <<g1[best_fitness_loc] <<"\n\t" <<g2[best_fitness_loc]
				<<"\n\t" <<g3[best_fitness_loc] <<"\n\t" <<g4[best_fitness_loc] <<"\n";
		}
};

// Constructor
DE::DE()
{
	getData();
	
	x1vector = new double[POPULATION+1];
	x2vector = new double[POPULATION+1];
	x3vector = new double[POPULATION+1];
	fitness = new double[POPULATION+1];
	g1 = new double[POPULATION+1];
	g2 = new double[POPULATION+1];
	g3 = new double[POPULATION+1];
	g4 = new double[POPULATION+1];
	
	for (int i=0; i<POPULATION; i++)
	{
		x1vector[i] = randomBetween(MIN_X1, MAX_X1);
		x2vector[i] = randomBetween(MIN_X2, MAX_X2);
		x3vector[i] = randomBetween(MIN_X3, MAX_X3);
	}		
}


// Get the DE parameters
void DE::getData()
{
	ifstream File("de_springwt.dat");
	File >>MIN_X1 >>MIN_X2 >>MIN_X3;
	File >>MAX_X1 >>MAX_X2 >>MAX_X3;
	File >>POPULATION >>F_MIN >> F_MAX >>CR_MIN >>CR_MAX;
	File.close();
}

// Find the fitnesses inclusive of penalty
// Also finds the best and average fitness
void DE::findFitness()
{
	int i, j;
	for (i=0; i<POPULATION; i++)
	{
		double x1=x1vector[i], x2=x2vector[i], x3=x3vector[i];
		fitness[i] = ( x3 + 2.0 )*x2*x1*x1;
		
		// calculate constriants
		g1[i] = 1.0 - ( x2*x2*x2*x3 ) / ( 71785.0*x1*x1*x1*x1 );
		g2[i] = ( 4.0*x2*x2 - x1*x2 ) / ( 12566.0 * (x2*x1*x1*x1 - x1*x1*x1*x1 ) ) + 1.0 / ( 5108.0 * x1*x1 ) - 1.0;
		g3[i] = 1.0 - ( 140.45*x1 ) / ( x2*x2*x3);
		g4[i] = ( x2 + x1 ) / 1.5 - 1.0;

		// add penalties to the fitness function
		if (g1[i]>0) fitness[i] += 1000;
		if (g2[i]>0) fitness[i] += 1000;
		if (g3[i]>0) fitness[i] += 1000;
		if (g4[i]>0) fitness[i] += 1000;
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
	int i=POPULATION;
	double x1=x1vector[i], x2=x2vector[i], x3=x3vector[i];
	fitness[i] = ( x3 + 2.0 )*x2*x1*x1;
	g1[i] = 1.0 - ( x2*x2*x2*x3 ) / ( 71785.0*x1*x1*x1*x1 );
	g2[i] = ( 4.0*x2*x2 - x1*x2 ) / ( 12566.0 * (x2*x1*x1*x1 - x1*x1*x1*x1 ) ) + 1.0 / ( 5108.0 * x1*x1 ) - 1.0;
	g3[i] = 1.0 - ( 140.45*x1 ) / ( x2*x2*x3 );
	g4[i] = ( x2 + x1 ) / 1.5 - 1.0;
	if (g1[i]>0) fitness[i] += 1000;
	if (g2[i]>0) fitness[i] += 1000;
	if (g3[i]>0) fitness[i] += 1000;
	if (g4[i]>0) fitness[i] += 1000;
}


// Main algorithm for evolution
void DE::evolution()
{
	int i, j, k;
	int GENERATION = 0;
	double temp_x1, temp_x2, temp_x3;
	double F, CR;
	
	avg_fitness = 1e50; best_fitness = 1e49;
	while (avg_fitness - best_fitness >= 0.00000001 && GENERATION < 100000)
	{
		findFitness();
		GENERATION++;

		cout <<"Gen.: " <<GENERATION <<"\tBest fit.: " <<best_fitness
			<<"\tAvg. fit.: " <<avg_fitness << "\t(x1,x2,x3) = " <<x1vector[best_fitness_loc] 
			<<"\t" <<x2vector[best_fitness_loc] <<"\t" <<x3vector[best_fitness_loc] <<endl;
		int r1 = randomBetween(0, POPULATION);
		int r2 = randomBetween(0, POPULATION);
		int r3 = randomBetween(0, POPULATION);
		int r4 = randomBetween(0, POPULATION);

		// Differential mutation
		F = randomBetween(F_MIN, F_MAX);	
		temp_x1 = x1vector[best_fitness_loc] + F * (x1vector[r1] - x1vector[r2]);
		temp_x2 = x2vector[best_fitness_loc] + F * (x2vector[r1] - x2vector[r2]);
		temp_x3 = x3vector[best_fitness_loc] + F * (x3vector[r1] - x3vector[r2]);
		if (temp_x1 < MIN_X1) temp_x1 = MIN_X1;
		if (temp_x2 < MIN_X2) temp_x2 = MIN_X2;
		if (temp_x3 < MIN_X3) temp_x3 = MIN_X3;

		// Crossing over
		CR = randomBetween(CR_MIN, CR_MAX);
		if ( randomBetween(0.0, 1.0) < CR )
		{
			x1vector[POPULATION] = temp_x1;
			x2vector[POPULATION] = temp_x2;
			x3vector[POPULATION] = temp_x3;
		}
		else
		{
			x1vector[POPULATION] = x1vector[r3];
			x2vector[POPULATION] = x2vector[r3];
			x3vector[POPULATION] = x3vector[r3];
		}
		
		// Recombination
		findFitnessLast();
		if (fitness[POPULATION] < fitness[r4])
		{
			x1vector[r4] = x1vector[POPULATION];
			x2vector[r4] = x2vector[POPULATION];
			x3vector[r4] = x3vector[POPULATION];
		}			
	}
}

#endif

#include <ctime>

int main()
{
	srand(time(0));
	cout.precision(8);
	DE de;
	de.evolution();
	de.printResult();
	getch();
}

