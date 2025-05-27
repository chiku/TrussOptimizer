// Optimization of truss using DE

#include <string>

#include "truss.h"

#ifndef __DE_H__
#define __DE_H__

int randomBetween(int a, int b); // [a, b)

double randomBetween(double a, double b); // [a, b)

class DE
{
    private:
        std::string truss_file;
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
        void getData(const std::string de_file);

    public:
        DE(const std::string de_file, const std::string truss_file);
        ~DE();
        void evolution();
        void printResult();
};


#endif
