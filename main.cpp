#include <string>
#include <ctime>

#include "de.h"

int main()
{
	time_t t; time(&t); srand(t); // random seed
	cout.precision(8);
	DE opti("data/de.dat", "data/truss1.dat");
	opti.evolution();
	opti.printResult();
}
