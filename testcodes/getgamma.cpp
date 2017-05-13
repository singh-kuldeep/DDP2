#include "iostream"
#include "vector"
#include "math.h"
using namespace std;
double getgamma(std::vector<double> U)
{
	double gamma_p = 1.4;
	double T_theta = 3055;
	double gamma_old = 1.0; 
	double gamma_new = 1.4 ; 
	double gammaTolarance = 1e-5;

	double Temperature, TemperatureRatio;

	while(fabs(gamma_old - gamma_new)>gammaTolarance)
	{	
		gamma_old = gamma_new;

		double R = 287.1;
		double Density = U[0];
		double pressure = (gamma_old-1)*(U[4]-(0.5*(pow(U[1],2)+pow(U[2],2)+pow(U[3],2))/Density)); 

		Temperature = pressure/(Density*R);
		cout << "Temperature " <<  Temperature << endl ;
		TemperatureRatio = T_theta/Temperature;
		gamma_new = 1 + (gamma_p - 1)/(1 + (gamma_p - 1)*(pow(TemperatureRatio,2)*
			exp(TemperatureRatio)/pow((exp(TemperatureRatio)-1),2)));
		cout << gamma_new << endl;
	}
	return gamma_new;
}
int main()
{
	double gamma;
	std::vector<double> U(5);
	double Density = 1.23;
	double XVelocity = 1000;
	double YVelocity = 00;
	double ZVelocity = 00;
	double pressure = 100e5;
	U[0] = Density;
	U[1] = Density*XVelocity;
	U[2] = Density*YVelocity;
	U[3] = Density*ZVelocity;
	U[4] =0.4*(pressure-0.5*Density*(pow(XVelocity,2)+pow(YVelocity,2)+pow(ZVelocity,2)));
	cout << getgamma(U) << endl;
	return 0;
}