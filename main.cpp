#include "de.cpp"

int main()
{
	cout.precision(7);
	DE opti;
	opti.getData();
	opti.evolution();
	opti.printResult();
	getch();
}
