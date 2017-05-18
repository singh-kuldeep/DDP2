/*! \file GetGamma.h 
\brief Contains the getgamma() function which calculates the specific heats 
ratio for given conserved quantities
*/ 
#ifndef GETGAMMA_H
#define GETGAMMA_H

/*! \fn double getgamma(std::vector<double> U)
\brief Function getgamma() calculates the gamma(T) for given conserved 
quantities in the cell  
\param [in] U Conserved quantities in the cell
*/
/*!
@image html gamma.png "Flow chart of varying gamma implementation" width=5cm
*/   
double getgamma(std::vector<double> U)
{
	double gamma_p = 1.4;
	/**\param gamma_p Reference specific heat ratio*/
	double T_theta = 3055;
	/**\param T_theta Reference temperature for gamma calculation*/
	double gamma_old = 1.0; 
	/**\param gamma_old Lower limit of specific heat ratio in the iteration*/
	double gamma_new = 1.4 ; 
	/**\param gamma_new Upper limit of specific heat ratio in the iteration*/
	double gammaTolerance = 1e-5;
	/**\param gammaTolerance Tolerance value for specific heat ratio in the 
	iteration*/
	double Temperature, TemperatureRatio;
	/**\param Temperature Temperature where gamma needs to be calculated*/

	// iterations for gamma calculation
	while(fabs(gamma_old - gamma_new)>gammaTolerance)
	{	
		gamma_old = gamma_new;

		double R = 287.1;
		double Density = U[0];
		double pressure = (gamma_old-1)*(U[4]-(0.5*(pow(U[1],2)+pow(U[2],2)+
			pow(U[3],2))/Density)); 

		Temperature = pressure/(Density*R);
		// cout << "Temperature " <<  Temperature << endl ;
		TemperatureRatio = T_theta/Temperature;
		gamma_new = 1 + (gamma_p - 1)/(1 + (gamma_p - 1)*
			(pow(TemperatureRatio,2)*exp(TemperatureRatio)/
				pow((exp(TemperatureRatio)-1),2)));
		// cout << gamma_new << endl;
	}
	return gamma_new;
}
#endif // getgamma.h ends here 