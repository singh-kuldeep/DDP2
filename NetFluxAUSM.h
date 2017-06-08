/*! \file NetFluxAUSM.h
\brief Calculates the AUSM net flux(convective and pressure) at a cell 
interface
\date 18-May-2017
 */

#ifndef NETFLUXAUSM
#define NETFLUXAUSM

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "math.h"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#include "ConvectiveFluxAUSM.h"

using namespace std ;
/*! \class netfluxAUSM
  \brief Calculates the net AUSM flux at the interface
 */
class netfluxAUSM
{
	public:
	/**\param NeTFlux[5] Net AUSM flux vector at the cell interface*/	
	double NetFlux[5] ;	

	/*! A constructor to calculate the net AUSM flux and to initiates the 
	parameter which are required in the flux calculation.
	\param LeftConservedVariable Conserved variable of the left cell
	\param RightConservedVariable Conserved variable of the right cell
	\param AreaVector Face area vector
	\param gamma String tells whether to consider specific heat ratio is 
	constant or varying	with temperature 
	\param SpecificHeatRatio Specific heat ratio in case of constant gamma
	*/	
	netfluxAUSM(
	 	vector<double> LeftConservedVariable,
	 	vector<double> RightConservedVariable,
		vector<double> AreaVector, string gamma, double SpecificHeatRatio)
	{
		double MachHalf;
		double PressureHalf;

		double AreaVectorNormal[5];
		
		double AreaVectorMagnitude = sqrt(pow(AreaVector[0],2) + 
			pow(AreaVector[1],2) + pow(AreaVector[2],2));
		
		AreaVectorNormal[0] = 0;
		AreaVectorNormal[1] = AreaVector[0]/AreaVectorMagnitude;
		AreaVectorNormal[2] = AreaVector[1]/AreaVectorMagnitude;
		AreaVectorNormal[3] = AreaVector[2]/AreaVectorMagnitude;
		AreaVectorNormal[4] = 0;

		eulerfluxAUSM left(LeftConservedVariable, AreaVector, gamma, 
			SpecificHeatRatio);
		eulerfluxAUSM right(RightConservedVariable, AreaVector, gamma, 
			SpecificHeatRatio);
		
		MachHalf = left.MachPlus + right.MachMinus;
		PressureHalf = left.PressurePlus + right.PressureMinus;
		

		for (int i = 0; i < 5; ++i)
		{
			NetFlux[i] = AreaVectorMagnitude*( 0.5*MachHalf*
				(left.Flux[i]+right.Flux[i]) -
				0.5*fabs(MachHalf)*(right.Flux[i] -
				left.Flux[i]) + PressureHalf*AreaVectorNormal[i]); 
		}
		/* cout << NetFlux[0] << " " <<  NetFlux[1] << " " <<
		 NetFlux[2] << " , " << NetFlux[3] << " " << NetFlux[4]<< endl;*/
	};
	// ~netfluxinterface();
	
};
#endif /* !NETFLUXAUSM */
