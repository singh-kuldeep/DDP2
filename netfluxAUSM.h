/* File foo.  */
#ifndef NETFLUXAUSM
#define NETFLUXAUSM

#include "math.h"
#include "eulerfluxAUSM.h"
// #include "netfluxBase.h"
// #include "diffusionfluxAUSM.h"

using namespace std ;
/*! \file 	   netfluxinterface.h
 *  \brief	   Calculates the net flux vector(numerical diffusion and euler 
 * flux) at the interface.
 *  \details   This class uses the two other class. One Euler for euler fulx
 * calculation and second for numerical diffusion flux calculation.
 *  \author    Kuldeep Singh
 *  \date      2017
 *  \copyright GNU Public License(GPL).
*  \param DiffusionFluxVector Numerical diffusion flux vector at the 
	interface
*  \param [in] ConservedVariable Conserved variable vector ([Density , 
	x-momentum, y-momentum, z-momentum, Energy])
*  \param [in] CellVulume Pointer to the cell volume vector
*  \param [in] LeftMinus Cell just previous to the left	
*  \param [in] RightPlus Cell just Next to the right
*  \param [in] DeltaT Time step		
 */
class netfluxAUSM
// : public netfluxBase
{
	public:
	double NetFlux[5] ;	
		
	netfluxAUSM(
	 	vector<double> LeftConservedVariable,
	 	vector<double> RightConservedVariable,
		vector<double> AreaVector, string gamma, double SpecificHeatRatio)
	{
		// double NetFlux[5] ; 
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

		eulerfluxAUSM left(LeftConservedVariable, AreaVector, gamma, SpecificHeatRatio);
		eulerfluxAUSM right(RightConservedVariable, AreaVector, gamma, SpecificHeatRatio);
		
		MachHalf = left.MachPlus + right.MachMinus;
		PressureHalf = left.PressurePlus + right.PressureMinus;
		

		for (int i = 0; i < 5; ++i)
		{
			NetFlux[i] = AreaVectorMagnitude*( 0.5*MachHalf*(left.Flux[i]+right.Flux[i]) - 
						0.5*fabs(MachHalf)*(right.Flux[i] - left.Flux[i]) + PressureHalf*AreaVectorNormal[i]); 
		}
		// cout << "new" << endl;
		// cout << NetFlux[0] << " " <<  NetFlux[1] << " " << NetFlux[2] << " , " << NetFlux[3] << " " << NetFlux[4]<< endl;
	};
	// ~netfluxinterface();
	
};
#endif /* !NETFLUXAUSM */
