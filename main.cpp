#include "de.cpp"

int main()
{
	cout.precision(8);
	DE opti;
	opti.getData();
	opti.evolution();
	getch();
}
