#include "math.h"
#include "iostream"
#include "getgamma.h"
#define SpecificHeatRatio 1.4 /*!< This is gas constant (Gamma). For air at room
 temperature it is almost equal to 1.4. If you are using some other 
 gas at some other temperature then change it*/

using namespace std ;
/*! 
 *	\file eulerflux.h
 *  \brief     This class calculates the euler flux vectors(Ee,Fe,Ge) at the 
 *interface.
 *  \author    Kuldeep Singh
 *  \date      2017
 *  \bug       Not all memory is freed when deleting an object of this class.
 *  \copyright GNU Public License.
*	\param Flux Euler flux vector at interface	
*	\param [in] AreaVectorNormal Interface area vector	
*	\param [in] ConservedVariable Conserved variable vector ([Density , 
* x-momentum, y-momentum, z-momentum, Energy])
*\param Pressure Satic pressure (p)
*/
class eulerfluxAUSM
{
	public:
	double Flux[5] ;
	double MachPlus;
	double MachMinus;
	double PressurePlus;
	double PressureMinus;
	double Mach;
	eulerfluxAUSM(vector<double> ConservedVariable, vector<double> AreaVector, string gamma)
	{
		double AreaVectorNormal[3];
		
		double AreaVectorMagnitude = sqrt(pow(AreaVector[0],2) + 
			pow(AreaVector[1],2) + pow(AreaVector[2],2));
		
		AreaVectorNormal[0] = AreaVector[0]/AreaVectorMagnitude;
		AreaVectorNormal[1] = AreaVector[1]/AreaVectorMagnitude;
		AreaVectorNormal[2] = AreaVector[2]/AreaVectorMagnitude;

		double Density = ConservedVariable[0];
		
		double VelocityNormal = (ConservedVariable[1]*AreaVectorNormal[0]+
			ConservedVariable[2]*AreaVectorNormal[1]+
			ConservedVariable[3]*AreaVectorNormal[2])/Density;
		
		double Pressure, VelocitySound;

		if (gamma == "Gamma(T)")
		{
			// gamma change
			Pressure = (getgamma(ConservedVariable) -1)*( ConservedVariable[4] - 0.5*(
			pow(ConservedVariable[1],2)+pow(ConservedVariable[2],2)+
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

		// Flux[3] = 0; // Delibreatly making z component zero

		Flux[4] = VelocitySound*(ConservedVariable[4] + Pressure) ;

		// Flux[4] = Density*VelocitySound*(ConservedVariable[4] + Pressure)/Density ; // H = (e+p)/rho;

	 };

};