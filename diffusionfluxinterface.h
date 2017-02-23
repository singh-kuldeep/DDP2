// checked
#include "math.h"
#include "interface.h"

using namespace std ;
/*! 
 *  \file      diffusionfluxinterface.h 
 *  \brief     This class calculates the numerical diffusion flux.
 *  \author    Kuldeep Singh
 *  \date      2017
 *  \copyright GNU Public License.
 */
class diffusionfluxinterface
{		
	public:
	double DiffusionFluxVector[5] ; 
	/**@param DiffusionFluxVector Numerical diffusion flux vector 
	at the interface*/
	/**\param [in] ConservedVariable Conserved variable vector ([Density , 
	x-momentum, y-momentum, z-momentum, Energy])*/
	/**@param [in] CellVulume Pointer to the cell volume vector*/
	/**@param [in] LeftMinus Cell just previous to the left*/	
	/**@param [in] RightPlus Cell just Next to the right*/
	/**\param [in] DeltaT Time step*/	
	diffusionfluxinterface(
		vector<double>& ConservedVariableLeftMinus, 
		vector<double>& ConservedVariableLeft,
		vector<double>& ConservedVariableRight,
		vector<double>& ConservedVariableRightPlus,
		vector<double>& FaceAreaVectorLeft,
		vector<double>& FaceAreaVectorRight,
		vector<double>& FaceAreaVectorRightPlus,
		double CellVolumeLeftMins,
		double CellVolumeLeft,
		double CellVolumeRight,
		double CellVolumeRightPlus,
		double DeltaT)
	{
		double ThetaI[5]; 
		double ThetaIPlus[5];
		double GVector[5] ;
		double GVectorPlus[5] ;
		double Signale ;
		double BetaInterface[5] ;
		double PhiInterface[5] ;
		double ZVectorInterface[5];
		double PshiInterface[5] ;
	
		interface left(
			ConservedVariableLeftMinus,
			ConservedVariableLeft,
			FaceAreaVectorLeft,
			CellVolumeLeftMins,
			CellVolumeLeft,
			DeltaT);

		interface right(
			ConservedVariableLeft,
			ConservedVariableRight,
			FaceAreaVectorRight,
			CellVolumeLeft,
			CellVolumeRight,
			DeltaT);
		interface rightplus(
			ConservedVariableRight,
			ConservedVariableRightPlus,
			FaceAreaVectorRightPlus,
			CellVolumeRight,
			CellVolumeRightPlus,
			DeltaT);

		// Defining the gi and gi+1
		for (int i = 0; i < 5; ++i)
		{
			if (right.GVectorInterface[i] >= 0)
			{
				Signale = 1.0;
			}
			else
			{
				Signale = -1.0 ;
			}
			/**@bug Here syntax needs to be changed for gvactor[i] calculation*/
			double tempB = (left.GVectorInterface[i]*Signale)  ;
			double tempA = fabs(right.GVectorInterface[i]) ;
			tempA = min(tempA,tempB) ;
			double temp_zero=0.0;
			tempB = max(temp_zero, tempA) ;
			GVector[i] = Signale *  tempB ;

			//GVector[i]=Signale*max(0.0, min(fabs(right.GVectorInterface[i]),
			//left.GVectorInterface[i]*Signale ) ) ;
		}

		 for (int i = 0; i < 5; ++i)
		{
			if (rightplus.GVectorInterface[i] >= 0)
			{
				Signale = 1.0;
			}
			else
			{
				Signale = -1.0 ;
			}

			double tempA = fabs(rightplus.GVectorInterface[i]) ;
			double tempB = (right.GVectorInterface[i]*Signale)  ;
			tempA = min(tempA,tempB) ;
			double temp_zero=0.0;
			tempB = max(temp_zero, tempA) ;
			GVectorPlus[i] = Signale * tempB ;
			//GVector[i]=Signale*max(0.0, min( fabs(right.GVectorInterface[i]),
			//left.GVectorInterface[i]*Signale ) ) ;
			
		 }
		//Defying theta i term 
		for (int i = 0; i < 5; ++i)
			{
				/**@bug Here I have doubt about "not equal to sign" because it 
				can't be exactly equal to 0.00000 so most of the 
				time we end up choosing theta i = 0.0*/
				if((fabs(right.AlphaVectorInterface[i]) + 
					fabs(left.AlphaVectorInterface[i]))!=0.0)
				{
					ThetaI[i] = fabs(right.AlphaVectorInterface[i] - 
					left.AlphaVectorInterface[i])/
					(fabs(right.AlphaVectorInterface[i]) + 
					fabs(left.AlphaVectorInterface[i])) ;
				}
				else
				{
					ThetaI[i] = 0.0;
				}
			}

		// theta i plus term
		for (int i = 0; i < 5; ++i)
		{
			// here I have doubt about not equal to soymbol because it cant be 
			//exect 0.00000 so most of the time we end up choosing theta i = 0.0

			if((fabs(rightplus.AlphaVectorInterface[i])+
				fabs(right.AlphaVectorInterface[i]))!=0.0)
			{
				ThetaIPlus[i]=fabs(rightplus.AlphaVectorInterface[i]-
				right.AlphaVectorInterface[i])/
				(fabs(rightplus.AlphaVectorInterface[i]) + 
				fabs(right.AlphaVectorInterface[i])) ;
			}
			else
			{
				ThetaIPlus[i] = 0.0 ;
			}
		}

		// Defying beta 
		// before that we need to define the omega values which are constant 
		double omega[5] = {0.25, 1.0, 1.0, 1.0, 0.25} ;
		for (int i = 0; i < 5; ++i)
		{
			BetaInterface[i] = 1.0 + omega[i]*(max(ThetaI[i],ThetaIPlus[i]));
		}
		/////////////////////////////////
		// Defying phi 
		for (int i = 0; i < 5; ++i)
		{
			if (right.AlphaVectorInterface[i] != 0)
			{
				PhiInterface[i] = (GVectorPlus[i] - GVector[i]) /
				right.AlphaVectorInterface[i] ;
			}
			else
			{
				PhiInterface[i] = 0.0 ;
			}
		}

		// now re-defying the Zinterface
		for (int i = 0; i < 5; ++i)
		{
			ZVectorInterface[i] = right.ZVectorInterface[i] + 
			BetaInterface[i]*PhiInterface[i] ;
		}

		// PshiInterface re-defying (for inviscid wall flow)
		for (int i = 0; i < 5; ++i)
		{
			PshiInterface[i] = pow(ZVectorInterface[i],2) + 0.25 ;
		}

		// re-defying the PshiInterface
		// // again we need to specify the deltaf = 0.2
		// 	double deltaf = 0.2 ;
		// 	for (int i = 0; i < 5; ++i)
		// 	{
		// 		if (ZVectorInterface[i] >= deltaf)
		// 		{
		// 			PshiInterface[i] = fabs(ZVectorInterface[i]) ;
		// 		}
		// 		else
		// 		{
		// 			PshiInterface[i] = 0.5*(pow(ZVectorInterface[i],2) + 
		//          pow(deltaf,2))/ deltaf ;
		// 		}
		// 	}

		// finally calculating the interface diffusion term =DiffusionFluxVector
		// before that we can calculate tempvector as
		double tempvector[5];
		for (int i = 0; i < 5; ++i)
		{
			tempvector[i] = ( BetaInterface[i]*(GVector[i]+GVectorPlus[i]) -
			 PshiInterface[i]*right.AlphaVectorInterface[i] )/DeltaT ;
		}
		// now the DiffusionFluxVector
		for (int i = 0; i < 5; ++i)
		{
			DiffusionFluxVector[i] = 0.0 ;
			for (int l = 0; l < 5; ++l)
			{
				DiffusionFluxVector[i] += right.EigenVectorMatrix[i][l] * 
				tempvector[l] ;
			}
		}
	};
	// ~diffusionfluxinterface();
};