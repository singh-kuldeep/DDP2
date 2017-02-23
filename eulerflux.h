#include "math.h"
#include "iostream"
#define SpecificHeatRatio 1.4 /*!< This is gas constant (Gamma). For air at room
 temperature it is almost equal to 1.4. If you are using some other 
 gas at some other temperature then change it*/

using namespace std ;
/*! 
 *	\file eulerflux.h
 *  \brief     This class calculates the euler flux vectors(Ee,Fe,Ge) at the interface.
 *  \author    Kuldeep Singh
 *  \date      2017
 *  \bug       Not all memory is freed when deleting an object of this class.
 *  \copyright GNU Public License.
 */
class eulerflux
{
	public:
	double EulerFluxX[5] ;
	double EulerFluxY[5] ;
	double EulerFluxZ[5] ;
	/**\param EulerFluxX x direction euler flux vector (Ee) at interface*/	
	/**\param EulerFluxY y direction euler flux vector (Fe) at interface*/
	/**\param EulerFluxZ z direction euler flux vector (Ge) at interface*/
	/**\param [in] ConservedVariable Conserved variable vector ([Density , x-momentum, y-momentum, z-momentum, Energy])*/
	/**\param Pressure Satic pressure (p)*/
	eulerflux(vector<double>& ConservedVariable)
	{
		double Pressure = (SpecificHeatRatio -1)*( ConservedVariable[4] - 0.5*(pow(ConservedVariable[1],2)+pow(ConservedVariable[2],2)+
		pow(ConservedVariable[3],2))/ConservedVariable[0] ) ;  

		//  Euler flux
		EulerFluxX[0] = ConservedVariable[1] ;
		EulerFluxX[1] = (pow(ConservedVariable[1],2)/ConservedVariable[0]) + Pressure ;     
		EulerFluxX[2] = ConservedVariable[1]*ConservedVariable[2]/ConservedVariable[0] ; 
		EulerFluxX[3] = ConservedVariable[1]*ConservedVariable[3] /ConservedVariable[0] ;
		EulerFluxX[4] = ((ConservedVariable[4] + Pressure)*ConservedVariable[1]) / ConservedVariable[0] ;

		EulerFluxY[0] = ConservedVariable[2] ;
		EulerFluxY[1] = ConservedVariable[1]*ConservedVariable[2]/ConservedVariable[0] ; 
		EulerFluxY[2] = (pow(ConservedVariable[2],2)/ConservedVariable[0]) + Pressure ;     
		EulerFluxY[3] = ConservedVariable[2]*ConservedVariable[3] /ConservedVariable[0] ;
		EulerFluxY[4] = ((ConservedVariable[4] + Pressure)*ConservedVariable[2]) / ConservedVariable[0] ;

		EulerFluxZ[0] = ConservedVariable[3] ;
		EulerFluxZ[1] = ConservedVariable[1]*ConservedVariable[3]/ConservedVariable[0] ; 
		EulerFluxZ[2] = ConservedVariable[2]*ConservedVariable[3] /ConservedVariable[0] ;
		EulerFluxZ[3] = (pow(ConservedVariable[3],2)/ConservedVariable[0]) + Pressure ;     
		EulerFluxZ[4] = ((ConservedVariable[4] + Pressure)*ConservedVariable[3]) / ConservedVariable[0] ;

	 };
	 // ~eulerflux();	
};