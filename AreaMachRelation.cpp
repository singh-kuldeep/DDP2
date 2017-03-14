#include "iostream"
#include "math.h"
using namespace std;

double getMach(double areaRatio, int &a)
{
	double mach_tolaerence = 1 ; // 0.001;
	double Mach_down = 1.0001;
	double Mach_up = 100;
	double Mach = 0.5*(Mach_up + Mach_down);
	double gamma = 1.4;
	while(mach_tolaerence > 0.001)
	{
		Mach = 0.5*(Mach_up + Mach_down);
		double area_difference = areaRatio - pow((gamma+1)/2,-((gamma+1)/(2*(gamma-1))))*pow((1+0.5*(gamma-1)*Mach*Mach),((gamma+1)/(2*(gamma-1))))/Mach;

		if(area_difference<0)
		{
			Mach_up = Mach;
		}
		else if(area_difference>=0)
		{
			Mach_down = Mach;
		}
		mach_tolaerence = Mach_up - Mach_down;		
	}
	a  = 5;
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
	int a =2 ;
	double Mach = 5;
	double areaRatio = 25;
	Mach = getMach(areaRatio,a) ;
	// cout <<"Area Ratio for Mach " << Mach << " is : " <<getAreaRatio(Mach)<<endl;
	cout <<"Mach for area ratio " << areaRatio << " is : " << Mach << "  a  "<< a << endl;
	// cout << "atan(1)  " << atan((1)/(1))*180/acos(-1) << endl;
	return 0;
}