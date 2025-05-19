#include <ctime>

#include "de.h"

int main()
{
	time_t t; time(&t); srand(t); // random seed
	cout.precision(8);
	DE opti;
	opti.getData();
	opti.evolution();
	opti.printResult();
}
