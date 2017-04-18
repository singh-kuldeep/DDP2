#include "math.h"
#include "iostream"
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
*	\param EulerFlux Euler flux vector at interface	
*	\param [in] AreaVector Interface area vector	
*	\param [in] ConservedVariable Conserved variable vector ([Density , 
* x-momentum, y-momentum, z-momentum, Energy])
*\param Pressure Satic pressure (p)
*/
class diffusionfluxAUSM
{
	public:
	double Flux[5] ;
		
	diffusionfluxAUSM(vector<double>& ConservedVariable, vector<double>& AreaVector)
	{
		double Pressure = (SpecificHeatRatio -1)*( ConservedVariable[4] - 0.5*(
		pow(ConservedVariable[1],2)+pow(ConservedVariable[2],2)+
		pow(ConservedVariable[3],2))/ConservedVariable[0] ) ;  

		double VelocityNormal = fabs((ConservedVariable[1]*AreaVector[0]+
					ConservedVariable[2]*AreaVector[1]+
					ConservedVariable[3]*AreaVector[2])/ConservedVariable[0]);
		double MachNormal = VelocityNormal/
		sqrt(SpecificHeatRatio*Pressure/ConservedVariable[0]);
		
		//  Euler flux
		Flux[0] = VelocityNormal*ConservedVariable[0] + 0 ;
		
		Flux[1] = VelocityNormal*ConservedVariable[1] + 
		Pressure*MachNormal*AreaVector[0];     
		
		Flux[2] = VelocityNormal*ConservedVariable[2] + 
		Pressure*MachNormal*AreaVector[1]; 
		
		Flux[3] = VelocityNormal*ConservedVariable[3] +
		Pressure*MachNormal*AreaVector[2];

		Flux[4] = VelocityNormal*ConservedVariable[4] + 
		Pressure*sqrt(SpecificHeatRatio*Pressure/ConservedVariable[0]) ;

	 };
	 // ~eulerflux();	
};