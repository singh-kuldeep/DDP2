#include <iostream>
#include <math.h>
using namespace std;

// void mirror(double &x,double&y, double a,double b,double c,double d,double l,double m)
// {
// 	double slop = (d-b)/(c-a);
// 	y = (slop*(2*l-2*a)+m*slop*slop+2*b-m)/(1+slop*slop);
// 	x =  (m-y)*slop+l;
// }

int main()
{

	// double gamma = 1.4;
	// double totalpressure = 79276600.8;
	// double inletpressure = 271767.0*(gamma-1);
	// double Mach = sqrt( (2/(gamma-1)) * ( pow( (totalpressure/inletpressure),((gamma-1)/gamma) ) -1 ) ) ;
	// cout << "Mach" << Mach << endl;
	
	double x=1.0,y=1.121512;
	// double a = 0.0, b = 0.0, c=2.0, d=0.0, l = 1.0, m=1.0;
	// double a = 0.0, b = 0.0, c=1.0, d=0.0, l = 1.0, m=0.0;
	// double a = 0.0, b = 0.0, c=1.0, d=0.0, l = 1.0, m=0.0;
	// mirror(x,y,a,b,c,d,l,m);
	cout << "mirrored points are:  " << x << ","  << y << endl;
	return 0;
}