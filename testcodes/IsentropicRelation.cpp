#include <iostream>
#include "math.h"
using namespace std;

int main()
{
	double SpecificHeatRatio = 1.4;
	double ExitMach = 5.0;
	double ExitPressure = 1e5; // N/m^2
	double ExitTemperature = 300; // K
	double TotalTemperature;
	double TotalPressure;

	TotalPressure = ExitPressure*
	pow((1+((SpecificHeatRatio-1)*pow(ExitMach,2)/2)),(SpecificHeatRatio/(SpecificHeatRatio-1)));

	TotalTemperature = ExitTemperature*(1+((SpecificHeatRatio-1)*pow(ExitMach,2)/2)); 

	double InletMach = 0.3174;
	double InletPressure, InletTemperature, InletVelocity, InletDensity;

	double ExitRP = 1 + ((SpecificHeatRatio-1)*pow(ExitMach,2)/2); 
	double InletRP = 1 + ((SpecificHeatRatio-1)*pow(InletMach,2)/2); 
	
	InletPressure =  ExitPressure*pow((ExitRP/InletRP),(SpecificHeatRatio)/(SpecificHeatRatio-1));
	InletTemperature = ExitTemperature*ExitRP/InletRP;
	InletDensity = InletPressure/(InletTemperature*287.14);
	InletVelocity = InletMach*sqrt(SpecificHeatRatio*InletPressure/InletDensity);

	// cout <<"InletDensity " << InletDensity << endl;
	// cout <<"InletVelocity " << InletVelocity <<endl;
	// cout << "InletPressure " << InletPressure << endl;

	// cout << "TotalPressure " << TotalPressure << endl;
	// cout << "TotalTemperature " << TotalTemperature << endl;

	double Mach = sqrt((2/(SpecificHeatRatio-1))*
	(pow((ExitPressure/TotalPressure),(-(SpecificHeatRatio-1)/SpecificHeatRatio))-1));
	cout << Mach << endl;
}