// checked
#include "math.h"
#include "iostream"

#define SpecificHeatRatio 1.4 /*!< This is gas constant (Gamma). For air at room
 temperature it is almost equal to 1.4. If you are using some other 
 gas at some other temperature then change it*/

using namespace std ;

/*! 
 *  \file      interface.h 
 *  \brief     This class calculates the interface parameters using 
 Reo scheme flux.
 *  \author    Kuldeep Singh
 *  \date      2017
 *  \copyright GNU Public License.
 *	\param DensityInterface Roe density at interface
 *	\param VelocityXInterface x velocity at interface
 *  \param VelocityYInterface y velocity at interface
 *  \param VelocityZInterface z velocity at interface
 *  \param EnthalpyInterface Enthalpy at interface
 *  \param EnthalpyInterface Enthalpy at interface
 *  \param VectorJumpInterface Change in the
	conserved parameters at the interface
 *  \param EigenValue Eigenvalue of the Jacobian matrix
 *  \param EigenVectorMatrix Eigenvector of the Jacobian matrix 
 *  \param EigenVectorMatrixInverse Inverse of the Jacobian matrix
 *  \param AlphaVectorInterface[5] 
	EigenVectorMatrixInverse[5][5]*VectorJumpInterface
 *  \param MuVectorInterface = delta t * EigenValue
 *  \param ZVectorInterface This is same as MuVectorInterface
 *  \param [in] ConservedVariables This is the pointer to the 4D vector 
	where all the conserved variables of previous time step are stored.
 *  \param [in] FaceAreaVectorInterface This is the pointer to the area 
	vector the cell interface
 *  \param CellVolume 3D vector which has the cell volume of all cells inside
	the domain	 	
 *  	 	
 *  	 	
 *  	 	
 */
class interface
{
	public:
	double DensityInterface ;
	double VelocityXInterface ;
	double VelocityYInterface ;
	double VelocityZInterface ; 
	double EnthalpyInterface ;
	

	double VectorJumpInterface[5] ; 
	double EigenValue[5] ;
	

	double EigenVectorMatrix[5][5] ; 
	
	double EigenVectorMatrixInverse[5][5] ; 
	

	double AlphaVectorInterface[5] ; 

	double MuVectorInterface[5] ; 
	
	double ZVectorInterface[5] ; 
	

	double PshiVectorInterface[5] ;
	double GVectorInterface[5];
	
	interface(
		vector<double>& ConservedVariableLeft,
		vector<double>& ConservedVariableRight,
		vector<double>& FaceAreaVectorInterface, 
		double CellVolumeLeft,	
		double CellVolumeRight, 
		double DeltaT) 
	{
		double CellVolumeInterface = 0.5*(CellVolumeLeft + CellVolumeRight) ;

		double PressureLeft = (SpecificHeatRatio-1)*(ConservedVariableLeft[4] - 
		0.5*(pow(ConservedVariableLeft[1],2)+pow(ConservedVariableLeft[2],2)+
		pow(ConservedVariableLeft[3],2))/ConservedVariableLeft[0] ) ;
		
		double PressureRight = (SpecificHeatRatio-1)*(ConservedVariableRight[4]-
		0.5*(pow(ConservedVariableRight[1],2)+pow(ConservedVariableRight[2],2)+
		pow(ConservedVariableRight[3],2))/ConservedVariableRight[0] ) ;
		
		double EnthalpyLeft = (ConservedVariableLeft[4] + PressureLeft)/
		ConservedVariableLeft[0];

		double EnthalpyRight = (ConservedVariableRight[4] + PressureRight)/
		ConservedVariableRight[0];

		DensityInterface = sqrt(ConservedVariableLeft[0]*
		ConservedVariableRight[0]) ; 
		
		VelocityXInterface = (ConservedVariableLeft[1] + 
		ConservedVariableRight[1]*sqrt(ConservedVariableLeft[0]/
		ConservedVariableRight[0]))/(ConservedVariableLeft[0] + sqrt(
		ConservedVariableRight[0]*ConservedVariableLeft[0])) ;
		
		VelocityYInterface = (ConservedVariableLeft[2] + 
		ConservedVariableRight[2]*sqrt(ConservedVariableLeft[0]/
		ConservedVariableRight[0]))/(ConservedVariableLeft[0]+sqrt(
		ConservedVariableRight[0]*ConservedVariableLeft[0])) ;
		
		VelocityZInterface = (ConservedVariableLeft[3] + 
		ConservedVariableRight[3]*sqrt(ConservedVariableLeft[0]/
		ConservedVariableRight[0]))/(ConservedVariableLeft[0] + sqrt(
		ConservedVariableRight[0]*ConservedVariableLeft[0])) ;
		
		EnthalpyInterface = (EnthalpyLeft + EnthalpyRight*sqrt(
		ConservedVariableRight[0]/ConservedVariableLeft[0])) /
		(1+sqrt(ConservedVariableRight[0]/ConservedVariableLeft[0])) ; 


		// Eigenvalue 
		double hx = FaceAreaVectorInterface[0] /CellVolumeInterface ;
		double hy = FaceAreaVectorInterface[1] /CellVolumeInterface ;
		
		double hz = FaceAreaVectorInterface[2] /CellVolumeInterface ;
		double hn = sqrt(pow(FaceAreaVectorInterface[0],2) + 
		pow(FaceAreaVectorInterface[1],2) + 
		pow(FaceAreaVectorInterface[2],2)) / CellVolumeInterface ;
		
		double Ucont = VelocityXInterface*hx + VelocityYInterface*hy +
		VelocityZInterface*hz ;
		
		double VelocitySoundInterface = sqrt((SpecificHeatRatio -1)*
		(EnthalpyInterface - 0.5*(pow(VelocityXInterface,2)+
		pow(VelocityYInterface,2) + pow(VelocityZInterface,2)))) ;

		EigenValue[0] = Ucont - VelocitySoundInterface*hn ;
		EigenValue[1] = Ucont ;
		EigenValue[2] = Ucont ;
		EigenValue[3] = Ucont ;
		EigenValue[4] = Ucont + VelocitySoundInterface*hn ;

		// Jump vector at interface 
		VectorJumpInterface[0] = CellVolumeInterface*(ConservedVariableRight[0]-
		ConservedVariableLeft[0]) ;
		
		VectorJumpInterface[1] = CellVolumeInterface*(ConservedVariableRight[1]-
		ConservedVariableLeft[1]) ;
		
		VectorJumpInterface[2] = CellVolumeInterface*(ConservedVariableRight[2]-
		ConservedVariableLeft[2]) ;
		
		VectorJumpInterface[3] = CellVolumeInterface*(ConservedVariableRight[3]-
		ConservedVariableLeft[3]) ;
		
		VectorJumpInterface[4] = CellVolumeInterface*(ConservedVariableRight[4]-
		ConservedVariableLeft[4]) ;


		// Defying the eigenvector matrix "R"
		double hdesx = hx/hn ; 
		double hdesy = hy/hn ; 
		double hdesz = hz/hn ; 
		double Phi = VelocityXInterface*hdesx + VelocityYInterface*hdesy + 
		VelocityZInterface*hdesz ;
		
		double q = sqrt(pow(VelocityXInterface,2)+pow(VelocityYInterface,2)+
		pow(VelocityZInterface,2)) ;

		EigenVectorMatrix[0][0] = 1 ; 
		EigenVectorMatrix[0][1] = 1 ;
		EigenVectorMatrix[0][2] = 0 ;
		EigenVectorMatrix[0][3] = 0 ; 
		EigenVectorMatrix[0][4] = 1 ;

		EigenVectorMatrix[1][0] = VelocityXInterface - hdesx * 
		VelocitySoundInterface ; 
		EigenVectorMatrix[1][1] = VelocityXInterface ;
		EigenVectorMatrix[1][2] = hdesy ;
		EigenVectorMatrix[1][3] = hdesz ; 
		EigenVectorMatrix[1][4] = VelocityXInterface + hdesx * 
		VelocitySoundInterface ;

		EigenVectorMatrix[2][0] = VelocityYInterface - hdesy * 
		VelocitySoundInterface ; 
		EigenVectorMatrix[2][1] = VelocityYInterface ;
		EigenVectorMatrix[2][2] = hdesz ;
		EigenVectorMatrix[2][3] = hdesx ; 
		EigenVectorMatrix[2][4] = VelocityYInterface + hdesy * 
		VelocitySoundInterface ;

		EigenVectorMatrix[3][0] = VelocityZInterface - hdesz * 
		VelocitySoundInterface ; 
		EigenVectorMatrix[3][1] = VelocityZInterface ;
		EigenVectorMatrix[3][2] = hdesx ;
		EigenVectorMatrix[3][3] = hdesy ; 
		EigenVectorMatrix[3][4] = VelocityZInterface + hdesz * 
		VelocitySoundInterface ;

		EigenVectorMatrix[4][0] = EnthalpyInterface - (hdesx*VelocityXInterface+
		hdesy*VelocityYInterface + hdesz*VelocityZInterface) * 
		VelocitySoundInterface ; 
		EigenVectorMatrix[4][1] = 0.5 * pow(q,2);
		EigenVectorMatrix[4][2] = hdesx*VelocityZInterface + 
		hdesz*VelocityYInterface + hdesy*VelocityXInterface;
		EigenVectorMatrix[4][3] = hdesy*VelocityZInterface + 
		hdesx*VelocityYInterface + hdesz*VelocityXInterface;
		EigenVectorMatrix[4][4] = EnthalpyInterface + 
		(hdesx*VelocityXInterface + hdesy*VelocityYInterface + 
		hdesz*VelocityZInterface) * VelocitySoundInterface ;

		// EigenVectorMatrixInverse defying
		EigenVectorMatrixInverse[0][0] = 0.5*(0.5*(pow(q,2))*
		((SpecificHeatRatio-1)/pow(VelocitySoundInterface,2))+ 
		(Phi/VelocitySoundInterface));
		EigenVectorMatrixInverse[0][1] = -0.5*(VelocityXInterface*
		((SpecificHeatRatio-1)/pow(VelocitySoundInterface,2))+ 
		(hdesx/VelocitySoundInterface));
		EigenVectorMatrixInverse[0][2] = -0.5*(VelocityYInterface*
		((SpecificHeatRatio-1)/pow(VelocitySoundInterface,2))+ 
		(hdesy/VelocitySoundInterface));
		EigenVectorMatrixInverse[0][3] = -0.5*(VelocityZInterface*
		((SpecificHeatRatio-1)/pow(VelocitySoundInterface,2))+ 
		(hdesz/VelocitySoundInterface));
		EigenVectorMatrixInverse[0][4] = 0.5*(SpecificHeatRatio -1)/ 
		pow(VelocitySoundInterface,2) ;

		EigenVectorMatrixInverse[1][0] = 1 - 0.5*pow(q,2)*((SpecificHeatRatio-
		1)/pow(VelocitySoundInterface,2)) ; 
		EigenVectorMatrixInverse[1][1] = VelocityXInterface*(SpecificHeatRatio-
		1)/pow(VelocitySoundInterface,2) ; 
		EigenVectorMatrixInverse[1][2] = VelocityYInterface*(SpecificHeatRatio-
		1)/pow(VelocitySoundInterface,2) ;   
		EigenVectorMatrixInverse[1][3] = VelocityZInterface*(SpecificHeatRatio-
		1)/pow(VelocitySoundInterface,2) ; 
		EigenVectorMatrixInverse[1][4] = -(SpecificHeatRatio-1)/
		pow(VelocitySoundInterface,2) ;

		EigenVectorMatrixInverse[2][0] = -(hdesy*VelocityXInterface + hdesz *
		VelocityYInterface + hdesx*VelocityZInterface);
		EigenVectorMatrixInverse[2][1] = hdesy ; 
		EigenVectorMatrixInverse[2][2] = hdesz ;
		EigenVectorMatrixInverse[2][3] = hdesx ; 
		EigenVectorMatrixInverse[2][4] = 0 ;

		EigenVectorMatrixInverse[3][0] = -(hdesz*VelocityXInterface + hdesx * 
		VelocityYInterface + hdesy*VelocityZInterface);
		EigenVectorMatrixInverse[3][1] = hdesz ; 
		EigenVectorMatrixInverse[3][2] = hdesx ;
		EigenVectorMatrixInverse[3][3] = hdesy ; 
		EigenVectorMatrixInverse[3][4] = 0 ;

		EigenVectorMatrixInverse[4][0] = 0.5*(0.5*(pow(q,2))*(
		(SpecificHeatRatio-1)/pow(VelocitySoundInterface,2)) - 
		(Phi/VelocitySoundInterface));
		EigenVectorMatrixInverse[4][1] = 0.5*( -VelocityXInterface*(
		(SpecificHeatRatio-1)/pow(VelocitySoundInterface,2))+ 
		(hdesx/VelocitySoundInterface));
		EigenVectorMatrixInverse[4][2] = 0.5*( -VelocityYInterface*(
		(SpecificHeatRatio-1)/pow(VelocitySoundInterface,2))+ 
		(hdesy/VelocitySoundInterface));
		EigenVectorMatrixInverse[4][3] = 0.5*( -VelocityZInterface*(
		(SpecificHeatRatio-1)/pow(VelocitySoundInterface,2))+ 
		(hdesz/VelocitySoundInterface));
		EigenVectorMatrixInverse[4][4] = 0.5*(SpecificHeatRatio -1)/
		pow(VelocitySoundInterface,2) ;

//AlphaVectorInterface defying 
 		for (int i = 0; i < 5; ++i)
		{
			AlphaVectorInterface[i] = 0.0 ;
			for (int l = 0; l < 5; ++l)
			{
				AlphaVectorInterface[i] += EigenVectorMatrixInverse[i][l] * 
				VectorJumpInterface[l] ;
			}
				// cout <<   AlphaVectorInterface[i] << endl ;
		}

// Defying the MuVectorInterface or ZVectorInterface
	for (int i = 0; i < 5; ++i)
		{
			MuVectorInterface[i] = DeltaT * EigenValue[i] ;
			ZVectorInterface[i] = DeltaT * EigenValue[i] ;
		}
			
// // PshiVectorInterface defying (for invisid flow)
		for (int i = 0; i < 5; ++i)
		{
			PshiVectorInterface[i] = pow(ZVectorInterface[i],2) + 0.25 ;
		}
// // PshiVectorInterface defying (for viscus flow)
// 	double deltaf = 0.2 ; 
// this is given constant value (it is between 0.1 to 0.5)
// 		for (int i = 0; i < 5; ++i)
// 		{
// 			if(fabs(ZVectorInterface[i]) >= deltaf ){
// 				PshiVectorInterface[i] = fabs(ZVectorInterface[i]);
// 			}
// 			else
// 			{
// 				PshiVectorInterface[i] = 0.5*(pow(ZVectorInterface[i],2)+
//				pow(deltaf,2)) / deltaf ;
// 			}
// 			cout<<"PshiVectorInterface ["<<i+1<<"]"<<PshiVectorInterface[i];
// 		}

// GVectorInterface defying
		for (int i = 0; i < 5; ++i)
		{
			GVectorInterface[i] = 0.5*(PshiVectorInterface[i] - 
				pow(ZVectorInterface[i],2))*AlphaVectorInterface[i];
		}
	 };
	// ~interface();
};