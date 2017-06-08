#ifndef CONVECTIVEFLUXAUSM_H
#define CONVECTIVEFLUXAUSM_H
/*! \file ConvectiveFluxAUSM.h 
\brief Calculates the convective flux of the AUSM scheme
\date 18-May-2017
*/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "math.h"
#include "iostream"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#include "GetGamma.h"

using namespace std ;

/*! 
	\class eulerfluxAUSM
	\brief This class calculates the Euler flux vectors(Ee,Fe,Ge) at the 
	interface for the AUSM scheme.
	\date 18-May-2017
*/

class eulerfluxAUSM
{
	public:
	/**\param Flux Convective flux vector in the cell*/	
	double Flux[5] ;
	/**\param MachPlus MachPlus needed in AUSM scheme, at a the cell*/	
	double MachPlus;
	/**\param MachMinus MachMinus needed in AUSM scheme, at a the cell*/	
	double MachMinus;
	/**\param PressurePlus PressurePlus needed in AUSM scheme, at a the cell*/		
	double PressurePlus;	
	/**\param PressureMinus PressureMinus needed in AUSM scheme, at a the cell*/	
	double PressureMinus;
	/**\param Mach Mach number at the cell*/	
	double Mach;

	/*! A constructor to calculate the convective flux and initiates the 
	parameter which are required in the flux calculation.
	\param ConservedVariable All the conserved variable in the cell
	\param AreaVector Face area vector
	\param gamma String tells whether to consider specific heat ratio is 
	constant or varying	with temperature 
	\param SpecificHeatRatio Specific heat ratio in case of constant gamma*/
	eulerfluxAUSM(vector<double> ConservedVariable, vector<double> AreaVector, 
		string gamma, double SpecificHeatRatio)
	{
		/**\param AreaVectorNormal Interface unit area vector*/	
		double AreaVectorNormal[3];
		
		/**\param AreaVectorMagnitude Magnitude of area vector*/	
		double AreaVectorMagnitude = sqrt(pow(AreaVector[0],2) + 
			pow(AreaVector[1],2) + pow(AreaVector[2],2));
		
		AreaVectorNormal[0] = AreaVector[0]/AreaVectorMagnitude;
		AreaVectorNormal[1] = AreaVector[1]/AreaVectorMagnitude;
		AreaVectorNormal[2] = AreaVector[2]/AreaVectorMagnitude;

		double Density = ConservedVariable[0];
		
		/**\param VelocityNormal Magnitude of velocity vector normal to the 
		interface, or <B>contravarient velocity*/	
		double VelocityNormal = (ConservedVariable[1]*AreaVectorNormal[0]+
			ConservedVariable[2]*AreaVectorNormal[1]+
			ConservedVariable[3]*AreaVectorNormal[2])/Density;
		
		double Pressure, VelocitySound;

		if (gamma == "Gamma(T)")
		{
			// gamma change
			Pressure = (getgamma(ConservedVariable) -1)*( ConservedVariable[4]- 
			0.5*(pow(ConservedVariable[1],2)+pow(ConservedVariable[2],2)+
			pow(ConservedVariable[3],2))/Density ) ;  
		}
		else // constant gamma
		{
			Pressure = (SpecificHeatRatio -1)*( ConservedVariable[4] - 0.5*(
			pow(ConservedVariable[1],2)+pow(ConservedVariable[2],2)+
			pow(ConservedVariable[3],2))/Density ) ;  
		}

		if (gamma == "Gamma(T)")
		{
			// gamma change
			VelocitySound = sqrt(getgamma(ConservedVariable)*Pressure/Density);
		}
		else
		{
			VelocitySound = sqrt(SpecificHeatRatio*Pressure/Density);
		}
		
		
		Mach = VelocityNormal/VelocitySound; 

		if(fabs(Mach)<=1)
		{
			MachPlus = 0.25*(Mach+1)*(Mach+1);
			MachMinus = -0.25*(Mach-1)*(Mach-1);

			PressurePlus = 0.25*Pressure*(Mach+1)*(Mach+1)*(2-Mach);
			PressureMinus = 0.25*Pressure*(Mach-1)*(Mach-1)*(2+Mach);
		}
		else
		{
			MachPlus = 0.5*(Mach+fabs(Mach));
			MachMinus = 0.5*(Mach-fabs(Mach));
			
			PressurePlus = 0.5*Pressure*(Mach+fabs(Mach))/Mach; 
			PressureMinus = 0.5*Pressure*(Mach-fabs(Mach))/Mach; 
		}


		//  Euler flux
		Flux[0] = Density*VelocitySound ;
		
		Flux[1] = VelocitySound*ConservedVariable[1];     
		
		Flux[2] = VelocitySound*ConservedVariable[2] ; 
		
		Flux[3] = VelocitySound*ConservedVariable[3] ;

		// Flux[3] = 0; // Deliberately making z component zero

		Flux[4] = VelocitySound*(ConservedVariable[4] + Pressure) ;
		// rho*V*H = rho*V*(e+p/rho) = V*(rho*e+p) = V*(E+p)
	 };

};
#endif // CONVECTIVEFLUXAUSM_H