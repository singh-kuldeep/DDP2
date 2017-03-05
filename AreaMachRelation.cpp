#include "iostream"
#include "math.h"
using namespace std;

double getMach(double areaRatio)
{
	double Mach = 100;
	return Mach;
}

double getAreaRatio(double Mach)
{
	double areaRatio = 100;
	double gamma = 1.4;
	 areaRatio = {[(gam+1)/2]^-[(gam+1)/(gam-1)/2]} / M * [1 + M^2 * (gam-1)/2]^[(gam+1)/(gam-1)/2]
	return areaRatio;
}

int main()
{
	double Mach = 5;
	cout <<"Area Ratio for Mach " << Mach << "is :" <<getAreaRatio(Mach)<<endl;
	return 0;
}