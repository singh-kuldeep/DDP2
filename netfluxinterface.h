#include "math.h"
#include "eulerflux.h"
#include "diffusionfluxinterface.h"

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
class netfluxinterface
{
	public:
	double NetFlux[5] ;
	
	netfluxinterface(
		vector<double>& ConservedVariableLeftMinus,
	 	vector<double>& ConservedVariableLeft,
	 	vector<double>& ConservedVariableRight,
		vector<double>& ConservedVariableRightPlus,
		vector<double>& FaceAreaLeft,
		vector<double>& FaceAreaVectorRight,
		vector<double>& FaceAreaVectorRightplus,
		double CellVolumeLeftMins,
		double CellVolumeLeft,
		double CellVolumeRight,
		double CellVolumeRightPlus,
		double DeltaT)
	{
		eulerflux left(ConservedVariableLeft);  
		/*This object is euler flux calculated using the left cell conserved variables*/
		eulerflux right(ConservedVariableRight); /* This object is
		euler flux calculated using the right cell conserved variables*/

		double CellVolumeInterface = (CellVolumeLeft + CellVolumeRight)/2 ; 
		/**\param CellVolumeInterface Average of left and right cell volume*/

		//Look at the documentation if diffusionfluxinterface class 
		//and eulerflux class 
		/*! 
			\sa diffusionfluxinterface()
			\sa eulerflux()
		*/
		diffusionfluxinterface diffusion(
			ConservedVariableLeftMinus,
			ConservedVariableLeft,
			ConservedVariableRight,
			ConservedVariableRightPlus,
			FaceAreaLeft,FaceAreaVectorRight,
			FaceAreaVectorRightplus,
			CellVolumeLeftMins,
			CellVolumeLeft,
			CellVolumeRight,
			CellVolumeRightPlus,
			DeltaT); 

		// Averaged interface Euler flux 
		double hx = FaceAreaVectorRight[0] / CellVolumeInterface ;
		double hy = FaceAreaVectorRight[1] / CellVolumeInterface ;	
		double hz = FaceAreaVectorRight[2] / CellVolumeInterface ;

		double NetFluxX[5] ; 
		double NetFluxY[5] ; 
		double ZNetFluxZ[5] ; 
		for (int i = 0; i < 5; ++i)
		{
			NetFluxX[i] = 0.5*(left.EulerFluxX[i]+right.EulerFluxX[i]) ; 
			NetFluxY[i] = 0.5*(left.EulerFluxY[i]+right.EulerFluxY[i]) ; 
			ZNetFluxZ[i] = 0.5*(left.EulerFluxZ[i]+right.EulerFluxZ[i]) ;

			NetFlux[i] = (NetFluxX[i]*hx + NetFluxY[i]*hy + ZNetFluxZ[i]*hz)*
			CellVolumeInterface + 0.5*diffusion.DiffusionFluxVector[i];
		}
	};
	// ~netfluxinterface();
	
};