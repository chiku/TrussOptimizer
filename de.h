// Optimization of truss using DE

#include "vendor/truss/include/truss.h"

#ifndef __DE_H__
#define __DE_H__

int randomBetween(int a, int b); // [a, b)

double randomBetween(double a, double b); // [a, b)

class DE
{
    private:
        Truss *T;
        double *fitness;
        bool *change;    // update the fitness of only if change equals true
        double avg_fitness, best_fitness;
        int best_fitness_loc;
        int total_evals; // total truss evaluations required in the whole process

        double MIN_AREA, MAX_AREA;             // geometric constraints
        double MAX_STRESS, MAX_DISP;           // behaviour constraints
        int POPULATION;                        // DE parameter
        double F_MIN, F_MAX, CR_MIN, CR_MAX;   // DE parameters
        double PENALTY;                        // constant penalty

    protected:
        void findFitness();
        void findFitnessLast();

    public:
        DE();
        ~DE();
        void getData();
        void evolution();
        void printResult();
};


#endif
