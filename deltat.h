/*! \file deltat.h
    \brief Contains the functions getLocalDeltaT() and getGlobalDeltaT(), which
    calculates the Local and global time step values respectively. 
    
    \date 18-May-2017 
*/
#ifndef DELTA_H
#define DELTA_H

#include <fstream>
#include "getgamma.h"

/*!\fn double getLocalDeltaT(vector<double> ConservedVariables, 
double MinimumDistance,double CFL, string gamma, double SpecificHeatRatio)
\brief Calculates the local time step value for a given cell
\param [IN] ConservedVariables Conserved variables at cell center of a 
particular cell
\param [IN] MinimumDistance Minimum distance of a given cell
\param [IN] CFL 
\param [IN] gamma String which tells that specific heat ratio is constant or 
function of temperature
\param SpecificHeatRatio Specific heat ratio 
*/
double getLocalDeltaT(vector<double> ConservedVariables, double MinimumDistance,
		double CFL, string gamma, double SpecificHeatRatio)
{
	double deltat = 1000;
	/**\param deltat time step*/ 

	double Density ;
	/**\param Density Density in the cell*/ 
	double Pressure ;
	/**\param Pressure Pressure in the cell*/ 
	double Velocity ;
	/**\param Velocity Velocity in the cell*/ 
	double VelocitySound;
	/**\param VelocitySound Sound velocity in the cell*/ 

	Density = ConservedVariables[0];

	// gamma is the function of temperature
	if(gamma == "Gamma(T)")
	{
		Pressure = (getgamma(ConservedVariables) -1)*
		(ConservedVariables[4] - 0.5*
		(pow(ConservedVariables[1],2)+
		pow(ConservedVariables[2],2)+
		pow(ConservedVariables[3],2))/Density); 
	}	
	else //constant 
	{
		Pressure = (SpecificHeatRatio -1)*
		(ConservedVariables[4] - 0.5*
		(pow(ConservedVariables[1],2)+
		pow(ConservedVariables[2],2)+
		pow(ConservedVariables[3],2))/Density);  
	}

	Velocity = sqrt(pow(ConservedVariables[1],2)+
	pow(ConservedVariables[2],2)+
	pow(ConservedVariables[3],2))/Density;
	
	if (gamma == "Gamma(T)")
	{
		// std::cout << gamma << endl;
		VelocitySound = sqrt(getgamma(ConservedVariables)*Pressure/Density);
	}
	else
	{
		// std::cout << gamma << endl;
		VelocitySound = sqrt(SpecificHeatRatio*Pressure/Density);
	}

	if(deltat > (CFL*MinimumDistance)/
		(Velocity+VelocitySound))
	{
		deltat = (CFL*MinimumDistance)/
		(Velocity+VelocitySound);
	}
	return deltat;
}	

/*!\fn double getGlobalDeltaT(
vector<vector<vector<vector<double> > > > ConservedVariables,
vector<vector<vector<double> > > MinimumDistance, double CFL, int Ni, int Nj, 
int Nk, string gamma, double SpecificHeatRatio)
\brief Calculates the global time step value at every time iteration
\param [IN] ConservedVariables Conserved variables at cell center of a 
particular cell
\param [IN] MinimumDistance Minimum distance of a given cell
\param [IN] CFL 
\param [IN] Ni Number of cells in in "i" direction.  
\param [IN] Nj Number of cells in in "j" direction.  
\param [IN] Nk Number of cells in in "k" direction.  
\param [IN] gamma String which tells that specific heat ratio is constant or 
function of temperature
\param SpecificHeatRatio Specific heat ratio 
*/
double getGlobalDeltaT(
vector<vector<vector<vector<double> > > > ConservedVariables,
vector<vector<vector<double> > > MinimumDistance, double CFL, int Ni, int Nj, 
int Nk, string gamma, double SpecificHeatRatio)
{
	double deltat = 1000.0;
	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{
				double LocalDeltaT = getLocalDeltaT(ConservedVariables[i][j][k],
				MinimumDistance[i][j][k],CFL, gamma, SpecificHeatRatio);
				if(deltat > LocalDeltaT)
				{
					deltat = LocalDeltaT;
				}
			}
		}
	}
}	
#endif // deltat.h ends here