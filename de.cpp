// Optimization of truss using DE

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "vendor/truss/include/truss.h"
#include "de.h"

inline int randomBetween(int a, int b) // [a, b)
{
	return ( a + rand() % (b+1-a) );
}

inline double randomBetween(double a, double b) // [a, b)
{
	return ( rand() * (b-a) /(0.0+RAND_MAX) + a);
}


// Constructor
DE::DE(const std::string de_file, const std::string truss_file)
{
	this->truss_file = truss_file;
	getData(de_file);

	T = new Truss[POPULATION+1];
	fitness = new double[POPULATION+1];
	change = new bool[POPULATION+1];
	for (int i=0; i<POPULATION+1; i++)
	{
		T[i].getData(truss_file.c_str());
		// Random initilization of area
		for (int j=0; j<T[i].members(); j++)
			T[i].area[j] = randomBetween(MIN_AREA, MAX_AREA);
		change[i] = true;
	}
	total_evals = 0;
}


DE::~DE()
{
	delete[] T;
	delete[] fitness;
	delete[] change;
}


// Get the DE parameters
void DE::getData(const std::string de_file)
{
	std::ifstream file(de_file);
	file >> MIN_AREA >> MAX_AREA;
	file >> MAX_STRESS >> MAX_DISP;
	file >> PENALTY;
	file >> POPULATION >> F_MIN >> F_MAX >> CR_MIN >> CR_MAX;
	file.close();
}

// Find the fitnesses inclusive of penalty
// Also finds the best and average fitness
void DE::findFitness()
{
	int i, j;
	for (i=0; i<POPULATION; i++)
	{
		// Solve and obtain the forces in the members
		if (change[i] == true)
		{
			T[i].findKLocal();
			T[i].findKGlobal();
			T[i].condense();
			T[i].solve();
			total_evals++;
			change[i] = false;

			fitness[i] = 0;
			for (j=0; j<T[0].members(); j++)
			{
				fitness[i] += T[i].area[j]*T[i].length[j];
				// stress penalty
				double t = sqrt( T[i].locforce[j](0, 0)*T[i].locforce[j](0, 0) +
					T[i].locforce[j](1, 0)*T[i].locforce[j](1, 0) );
				if ( fabs(t / T[i].area[j]) >= MAX_STRESS)
					fitness[i] += PENALTY;
			}

			// displacement penalty
			for (j=0; j<T[0].uglobal.rows(); j++)
				if ( fabs(T[i].uglobal(j, 0)) >= MAX_DISP )
					fitness[i] += PENALTY;
		}
	}

	best_fitness = fitness[0];
	best_fitness_loc = 0;
	avg_fitness = 0;
	for (i=0; i<POPULATION; i++)
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
	total_evals++;

	fitness[i] = 0;
	for (j=0; j<T[i].members(); j++)
	{
		fitness[i] += T[i].area[j]*T[i].length[j];
		// stress penalty
		double t = sqrt( T[i].locforce[j](0, 0)*T[i].locforce[j](0, 0) +
			T[i].locforce[j](1, 0)*T[i].locforce[j](1, 0) );

		if ( fabs(t / T[i].area[j]) >= MAX_STRESS)
			fitness[i] += PENALTY;
	}
	// displacement penalty
	for (j=0; j<T[i].uglobal.rows(); j++)
		if ( fabs(T[i].uglobal(j, 0)) >= MAX_DISP)
			fitness[i] += PENALTY;
}



// Main algorithm for evolution
void DE::evolution()
{
	int j;

	Truss temp;
	temp.getData(this->truss_file.c_str());
	long int GENERATION = 0;
	avg_fitness = 1e50; best_fitness = 1e49;
	while (avg_fitness - best_fitness >= 0.000001 && GENERATION < 50000)
	{
		findFitness();
		GENERATION++;

		std::cout << "Generation: " << GENERATION << "\tBest fit.: " << best_fitness
			<< "\tAvg. fit.: " << avg_fitness << "\n";
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

		// Crossing over
		for (j=0; j<=T[0].members(); j++)
		{
			double CR = randomBetween(CR_MIN, CR_MAX);
			if ( randomBetween(0.0, 1.0) < CR )
				T[POPULATION].area[j] = temp.area[j];
			else
				T[POPULATION].area[j] = T[r4].area[j];
		}


		// Recombination
		findFitnessLast();
		if (fitness[POPULATION] < fitness[r5])
		{
			for (j=0; j<T[0].members(); j++)
				T[r5].area[j] = T[POPULATION].area[j];
			change[r5] = true; // would be requiring a change
		}
	}
}


void DE::printResult()
{
	int i;
	std::cout << "\n\nThe stress in members are\n";
	for (i=0; i<T[0].members(); i++)
	{
		double t =  T[best_fitness_loc].locforce[i](0, 0) * T[best_fitness_loc].locforce[i](0, 0) +
					T[best_fitness_loc].locforce[i](1, 0) * T[best_fitness_loc].locforce[i](1, 0);
		std::cout << "\t" << i << ">  " << sqrt(t)/T[best_fitness_loc].area[i] << "\n";
	}
	std::cout <<"\n\nThe nodal displacements are\n";
	for (i=0; i<T[0].uglobal.rows(); i++) {
		std::cout << "\t" << i << "> " << T[best_fitness_loc].uglobal(i, 0) << "\n";
	}
	std::cout <<"\n\nThe member area are\n";
	for (i=0; i<T[0].members(); i++) {
		std::cout << "Member no " << i << " " << T[best_fitness_loc].area[i] << "\n";
	}

	std::cout << "\n\nTotal truss evaluations performed: " << total_evals << "\n";
}
