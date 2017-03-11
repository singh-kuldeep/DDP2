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
	double areaRatio ;
	double gamma = 1.4;
	areaRatio = pow((gamma+1)/2,-((gamma+1)/(2*(gamma-1))))*pow((1+0.5*(gamma-1)*Mach*Mach),((gamma+1)/(2*(gamma-1))))/Mach;
	return areaRatio;
}

int main()
{
	double Mach = 5;
	cout <<"Area Ratio for Mach " << Mach << " is : " <<getAreaRatio(Mach)<<endl;
	return 0;
}