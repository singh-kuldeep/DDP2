#include <iostream>
#include <math.h>
using namespace std;


int main()
{
	// arae ratio = 32
	// Mexit = 5.3 

	double gamma = 1.4;
	double totalpressure = 79276600.8;
	double inletpressure = 271767.0*(gamma-1);
	double Mach = sqrt( (2/(gamma-1)) * ( pow( (totalpressure/inletpressure),((gamma-1)/gamma) ) -1 ) ) ;
	cout << "Mach" << Mach << endl;
	return 0;
}