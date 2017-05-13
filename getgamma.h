#ifndef GETGAMMA_H
#define GETGAMMA_H

// this function to calculate the gamma  
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
		// cout << "Temperature " <<  Temperature << endl ;
		TemperatureRatio = T_theta/Temperature;
		gamma_new = 1 + (gamma_p - 1)/(1 + (gamma_p - 1)*(pow(TemperatureRatio,2)*
			exp(TemperatureRatio)/pow((exp(TemperatureRatio)-1),2)));
		// cout << gamma_new << endl;
	}
	return gamma_new;
}
#endif // getgamma.h ends here 