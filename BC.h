/*! \file BC.h
    \brief This header file implements all three boundary conditions. 
    - Inlet
    - Exit and 
    - Wall boundary

    \author Kuldeep Singh
    \date 2015 
*/

// #include "iostream"
#include "math.h"
#include <vector>

#define SpecificHeatRatio 1.4 /*!< This is gas constant (Gamma). For air at 
room temperature it is almost equal to 1.4. If you are using some other gas at
some other temperature then change it*/
#define IdealGasConstant 287.14 
/*!< This is ideal gas constant \f$ R(J Kg^{-1}K^{-1}) = (c_p - c_v)\f$ */
#define SpecificHeatVolume 717.5 
/*!< This is specific heat at constant volume  of air (\f$ C_v \f$)*/

/** \brief Changes the input vector into the unit normal vector.
*\param areaVector A 3D vector.
*\param vectorMagnitude Magnitude of the 3D vector.
*\return void
*/
void getNormal(vector<double> & areaVector)
{
	double vectorMagnitude = sqrt(areaVector[0]*areaVector[0] + areaVector[1]*
		areaVector[1] + areaVector[2]*areaVector[2]);
	areaVector[0] = areaVector[0]/vectorMagnitude;
	areaVector[1] = areaVector[1]/vectorMagnitude;
	areaVector[2] = areaVector[2]/vectorMagnitude;
}

/** \brief Function BC() implements the boundary condition. 
* Here two ghost cell are used to implement the boundary condition. In simple 
words this function calculates the conserved variables for all ghost cells.
*For inlet it uses the stagnation parameters, for exit it simply uses the live 
cell parameters and copies them into the ghost cells, and for wall boundary
* it uses the fact that flow should be parallel to the wall. 
*\param [in] ConservedVariables This is the pointer to the 4D vector where all 
the conserved variables of previous time step are stored.
*\param [in] &iFaceAreaVector This is a pointer to the 4D vector which has the 
area vector of all faces which are in "i" direction.  
*\param [in] &jFaceAreaVector This is a pointer to the 4D vector which has the 
area vector of all faces which are in "j" direction.  
*\param [in] &kFaceAreaVector This is a pointer to the 4D vector which has the 
area vector of all faces which are in "k" direction.  
*\param [in] Ni Number of cells in in "i" direction.  
*\param [in] Nj Number of cells in in "j" direction.  
*\param [in] Nk Number of cells in in "k" direction.  
*\return void
*/
void BC(
	vector<vector<vector<vector<double> > > > & ConservedVariables,
	vector<vector<vector<vector<double> > > > & jFaceAreaVector,
	vector<vector<vector<vector<double> > > > & kFaceAreaVector,
	int Ni, int Nj, int Nk)
{
	/**
     * Inlet conditions are user given data.
     * one has to mention the stagnation parameters at inlet (ex. stagnation 
     pressure (\f$ P_0 \f$), temperature(\f$ T_0 \f$))
     */
	double TemperatureStagnation = 518.76 ; 
	/**\param TemperatureStagnation Stagnation temperature at inlet */  
	double PressureStagnation = 792766.8; 
	/**\param PressureStagnation Stagnation pressure at inlet */
	double DensityStagnation = PressureStagnation /
		(IdealGasConstant*TemperatureStagnation) ; 
	/**\param DensityStagnation Stagnation density at inlet */
	double theta = 3.14159265 * 0 / 180 ;
	/**\param Geometry rotation angle */
	 	

	/* Inlet ghost cells are being updated using the stagnation quantities
	(\f$ P_0, T_0 \f$) and flow direction */
	for (int j =2; j < Nj-2; ++j)
	{
		for (int k =2; k < Nk-2; ++k)
		{
			/**\param InletPressure Static pressure at inlet */
			double InletPressure = 
			(SpecificHeatRatio-1)*(ConservedVariables[2][j][k][4] - 0.5*(
				pow(ConservedVariables[2][j][k][1],2) +
				pow(ConservedVariables[2][j][k][2],2) + 
				pow(ConservedVariables[2][j][k][3],2)) /
				ConservedVariables[2][j][k][0]) ;

			/**\param Mach Mach number at inlet */
			double Mach=sqrt((2/(SpecificHeatRatio-1))*(pow((PressureStagnation/
			InletPressure),((SpecificHeatRatio-1)/SpecificHeatRatio) ) -1 ) ) ;
			/**\param InletTemperature Static temperature at inlet */
			double InletTemperature = TemperatureStagnation/(1+
			((SpecificHeatRatio-1)*Mach*Mach)/2);
			
			/**\param InletVelocity Flow velocity at inlet */
			double InletVelocity = Mach * sqrt(SpecificHeatRatio*
			IdealGasConstant*InletTemperature) ;
			
			/**\param InletDensity Flow density at inlet */
			double InletDensity = DensityStagnation / 
			pow((PressureStagnation/InletPressure),(1/SpecificHeatRatio)) ;

			ConservedVariables[0][j][k][0] = InletDensity ;
			ConservedVariables[0][j][k][1] = cos(theta)*
			InletDensity*InletVelocity ;
			ConservedVariables[0][j][k][2] = 0 ;
			ConservedVariables[0][j][k][3] = -sin(theta)*InletDensity*
			InletVelocity ;
			ConservedVariables[0][j][k][4] = InletPressure/(SpecificHeatRatio-1)
			+0.5*InletDensity*InletVelocity*InletVelocity ;

			ConservedVariables[1][j][k][0] = InletDensity ;
			ConservedVariables[1][j][k][1] = cos(theta)*InletDensity*
			InletVelocity;
			ConservedVariables[1][j][k][2] = 0 ;
			ConservedVariables[1][j][k][3] = -sin(theta)*InletDensity*
			InletVelocity;
			ConservedVariables[1][j][k][4] = InletPressure/(SpecificHeatRatio-1)
			+0.5*InletDensity*InletVelocity*InletVelocity ;
		}
	}

	/* At exit updating the i ghost cells (this is true where flow is 
	supersonic)*/
	for (int j =2; j < Nj-2; ++j)
	{
		for (int k =2; k < Nk-2; ++k)
		{
			for (int l = 0; l < 5; ++l)
			{
				ConservedVariables[Ni-2][j][k][l] = 
				ConservedVariables[Ni-3][j][k][l] ;
				ConservedVariables[Ni-1][j][k][l] = 
				ConservedVariables[Ni-4][j][k][l] ;
			}
		}
	}

	
	/* Updating the ghost cell conserved parameters value at j - wall */
	for (int i = 2; i < Ni-2; ++i)
	{
		for (int k = 2; k < Nk-2; ++k)
		{
			double jFaceAreaVector_unit_vector[3] = {0.0,2.0,0.0} ;

			//using leftcell (j=2) filling j=1
			jFaceAreaVector_unit_vector[0] = jFaceAreaVector[i][2][k][0] / sqrt(
			pow(jFaceAreaVector[i][2][k][0] ,2) + 
			pow(jFaceAreaVector[i][2][k][1] ,2) + 
			pow(jFaceAreaVector[i][2][k][2] ,2) ) ;

			jFaceAreaVector_unit_vector[1] = jFaceAreaVector[i][2][k][1] / sqrt(
			pow(jFaceAreaVector[i][2][k][0] ,2) + 
			pow(jFaceAreaVector[i][2][k][1] ,2) + 
			pow(jFaceAreaVector[i][2][k][2] ,2) ) ;

			jFaceAreaVector_unit_vector[2] = jFaceAreaVector[i][2][k][2] / sqrt(
			pow(jFaceAreaVector[i][2][k][0] ,2) +
			pow(jFaceAreaVector[i][2][k][1] ,2) +
			pow(jFaceAreaVector[i][2][k][2] ,2) ) ;

			ConservedVariables[i][1][k][1] = (1-2*
				pow(jFaceAreaVector_unit_vector[0],2))*
				ConservedVariables[i][2][k][1]-
				(2*jFaceAreaVector_unit_vector[0]*
				jFaceAreaVector_unit_vector[1])*
				ConservedVariables[i][2][k][2]-
				(2*jFaceAreaVector_unit_vector[0]*
				jFaceAreaVector_unit_vector[2])*
				ConservedVariables[i][2][k][3];

			ConservedVariables[i][1][k][2] = -(2*jFaceAreaVector_unit_vector[1]*
			jFaceAreaVector_unit_vector[0])*ConservedVariables[i][2][k][1]+
			(1-2*pow(jFaceAreaVector_unit_vector[1],2))*
			ConservedVariables[i][2][k][2] - (2*jFaceAreaVector_unit_vector[1]*
			jFaceAreaVector_unit_vector[2])*ConservedVariables[i][2][k][3] ;

			ConservedVariables[i][1][k][3] = -(2*jFaceAreaVector_unit_vector[2]*
			jFaceAreaVector_unit_vector[0])*ConservedVariables[i][2][k][1]-
			(2*jFaceAreaVector_unit_vector[2]*jFaceAreaVector_unit_vector[1])*
			ConservedVariables[i][2][k][2] + (1-2*pow(
			jFaceAreaVector_unit_vector[2],2))*ConservedVariables[i][2][k][3] ;

			// using zero order extrapolation
			ConservedVariables[i][1][k][0] = ConservedVariables[i][1][k][0] ;
			ConservedVariables[i][1][k][4] = ConservedVariables[i][2][k][4] ;


			//using rightcell(j=3) filling j=0
			jFaceAreaVector_unit_vector[0] = jFaceAreaVector[i][3][k][0] / sqrt(
			pow(jFaceAreaVector[i][3][k][0] ,2) +
			pow(jFaceAreaVector[i][3][k][1] ,2) +
			pow(jFaceAreaVector[i][3][k][2] ,2) ) ;  
			
			jFaceAreaVector_unit_vector[1] = jFaceAreaVector[i][3][k][1] / sqrt(
			pow(jFaceAreaVector[i][3][k][0] ,2) +
			pow(jFaceAreaVector[i][3][k][1] ,2) +
			pow(jFaceAreaVector[i][3][k][2] ,2) ) ;
			
			jFaceAreaVector_unit_vector[2] = jFaceAreaVector[i][3][k][2] / sqrt(
			pow(jFaceAreaVector[i][3][k][0] ,2) +
			pow(jFaceAreaVector[i][3][k][1] ,2) +
			pow(jFaceAreaVector[i][3][k][2] ,2) ) ;

			ConservedVariables[i][0][k][1] = (1-2*pow(
			jFaceAreaVector_unit_vector[0],2)) * ConservedVariables[i][3][k][1]-
			(2*jFaceAreaVector_unit_vector[0]*jFaceAreaVector_unit_vector[1])*
			ConservedVariables[i][3][k][2] - (2*jFaceAreaVector_unit_vector[0]*
			jFaceAreaVector_unit_vector[2])*ConservedVariables[i][3][k][3] ;

			ConservedVariables[i][0][k][2] = -(2*jFaceAreaVector_unit_vector[1]*
			jFaceAreaVector_unit_vector[0])*ConservedVariables[i][3][k][1]+(1-2*
			pow(jFaceAreaVector_unit_vector[1],2))*ConservedVariables[i][3][k][2]
			-(2*jFaceAreaVector_unit_vector[1]*jFaceAreaVector_unit_vector[2])*
			ConservedVariables[i][3][k][3] ;

			ConservedVariables[i][0][k][3] = -(2*jFaceAreaVector_unit_vector[2]*
			jFaceAreaVector_unit_vector[0])*ConservedVariables[i][3][k][1] -
			(2*jFaceAreaVector_unit_vector[2]*jFaceAreaVector_unit_vector[1])*
			ConservedVariables[i][3][k][2] + ( 1- 2* pow(
			jFaceAreaVector_unit_vector[2],2))*ConservedVariables[i][3][k][3] ;

			// using zero order extrapolation
			ConservedVariables[i][0][k][0] = ConservedVariables[i][3][k][0] ;
			ConservedVariables[i][0][k][4] = ConservedVariables[i][3][k][4] ;


			//using rightcell(j=Nj-3) filling j=Nj-2
			//here every term has been multiplied by -1 because 
			//cell area vector is -A
			jFaceAreaVector_unit_vector[0] = - jFaceAreaVector[i][Nj-3][k][0] / 
			sqrt( pow(jFaceAreaVector[i][Nj-3][k][0] ,2) + 
			pow(jFaceAreaVector[i][Nj-3][k][1] ,2) + 
			pow(jFaceAreaVector[i][Nj-3][k][2] ,2) ) ;  
			jFaceAreaVector_unit_vector[1] = - jFaceAreaVector[i][Nj-3][k][1] / 
			sqrt( pow(jFaceAreaVector[i][Nj-3][k][0] ,2) + 
			pow(jFaceAreaVector[i][Nj-3][k][1] ,2) + 
			pow(jFaceAreaVector[i][Nj-3][k][2] ,2) ) ;
			jFaceAreaVector_unit_vector[2] = - jFaceAreaVector[i][Nj-3][k][2] / 
			sqrt( pow(jFaceAreaVector[i][Nj-3][k][0] ,2) + 
			pow(jFaceAreaVector[i][Nj-3][k][1] ,2) + 
			pow(jFaceAreaVector[i][Nj-3][k][2] ,2) ) ;
			
			ConservedVariables[i][Nj-2][k][1] = (1-2 * 
			pow(jFaceAreaVector_unit_vector[0],2)) * 
			ConservedVariables[i][Nj-3][k][1] -
			(2*jFaceAreaVector_unit_vector[0] *  
			jFaceAreaVector_unit_vector[1])*
			ConservedVariables[i][Nj-3][k][2] -
			(2*jFaceAreaVector_unit_vector[0]*
			jFaceAreaVector_unit_vector[2])*ConservedVariables[i][Nj-3][k][3] ;

			ConservedVariables[i][Nj-2][k][2] = -(2*
			jFaceAreaVector_unit_vector[1] * 
			jFaceAreaVector_unit_vector[0]) * 
			ConservedVariables[i][Nj-3][k][1] +
			(1-2*pow(jFaceAreaVector_unit_vector[1],2)) * 
			ConservedVariables[i][Nj-3][k][2] -
			(2*jFaceAreaVector_unit_vector[1]*
			jFaceAreaVector_unit_vector[2]) * 
			ConservedVariables[i][Nj-3][k][3] ;

			ConservedVariables[i][Nj-2][k][3] = -(2*
			jFaceAreaVector_unit_vector[2] * 
			jFaceAreaVector_unit_vector[0]) * 
			ConservedVariables[i][Nj-3][k][1] -
			(2*jFaceAreaVector_unit_vector[2] *
			jFaceAreaVector_unit_vector[1]) * 
			ConservedVariables[i][Nj-3][k][2] +
			(1-2*pow(jFaceAreaVector_unit_vector[2],2)) * 
			ConservedVariables[i][Nj-3][k][3] ;

			// using zero order extrapolation
			ConservedVariables[i][Nj-2][k][0] = 
			ConservedVariables[i][Nj-3][k][0] ;
			ConservedVariables[i][Nj-2][k][4] = 
			ConservedVariables[i][Nj-3][k][4] ;


			// using rightcell(j=Nj-4) filling j=Nj-1
			// here every term has been multiplied by -1 
			// because cell area vector is -A
			jFaceAreaVector_unit_vector[0] = - jFaceAreaVector[i][Nj-4][k][0] / 
			sqrt( pow(jFaceAreaVector[i][Nj-4][k][0] ,2) + 
			pow(jFaceAreaVector[i][Nj-4][k][1] ,2) + 
			pow(jFaceAreaVector[i][Nj-4][k][2] ,2) ) ;  
			
			jFaceAreaVector_unit_vector[1] = - jFaceAreaVector[i][Nj-4][k][1] / 
			sqrt( pow(jFaceAreaVector[i][Nj-4][k][0] ,2) + 
			pow(jFaceAreaVector[i][Nj-4][k][1] ,2) + 
			pow(jFaceAreaVector[i][Nj-4][k][2] ,2) ) ;
			
			jFaceAreaVector_unit_vector[2] = - jFaceAreaVector[i][Nj-4][k][2] / 
			sqrt( pow(jFaceAreaVector[i][Nj-4][k][0] ,2) + 
			pow(jFaceAreaVector[i][Nj-4][k][1] ,2) + 
			pow(jFaceAreaVector[i][Nj-4][k][2] ,2) ) ;
			

			ConservedVariables[i][Nj-1][k][1] = (1-2*
			pow(jFaceAreaVector_unit_vector[0],2)) * 
			ConservedVariables[i][Nj-4][k][1] -
			(2*jFaceAreaVector_unit_vector[0] * 
			jFaceAreaVector_unit_vector[1]) * 
			ConservedVariables[i][Nj-4][k][2] -
			(2*jFaceAreaVector_unit_vector[0] * 
			jFaceAreaVector_unit_vector[2]) * 
			ConservedVariables[i][Nj-4][k][3] ;

			ConservedVariables[i][Nj-1][k][2] = -(2*
			jFaceAreaVector_unit_vector[1] * 
			jFaceAreaVector_unit_vector[0]) * 
			ConservedVariables[i][Nj-4][k][1] +
			(1-2*pow(jFaceAreaVector_unit_vector[1],2)) * 
			ConservedVariables[i][Nj-4][k][2] -
			(2*jFaceAreaVector_unit_vector[1] * 
			jFaceAreaVector_unit_vector[2]) * 
			ConservedVariables[i][Nj-4][k][3] ;

			ConservedVariables[i][Nj-1][k][3] = -(2 *
			jFaceAreaVector_unit_vector[2]*
			jFaceAreaVector_unit_vector[0])* 
			ConservedVariables[i][Nj-4][k][1]-
			(2*jFaceAreaVector_unit_vector[2]* 
			jFaceAreaVector_unit_vector[1]) * 
			ConservedVariables[i][Nj-4][k][2]+
			(1-2*pow(jFaceAreaVector_unit_vector[2],2))
			*ConservedVariables[i][Nj-4][k][3] ;

			// using zero order extrapolation
			ConservedVariables[i][Nj-1][k][0] = 
			ConservedVariables[i][Nj-4][k][0] ;
			ConservedVariables[i][Nj-1][k][4] = 
			ConservedVariables[i][Nj-4][k][4] ;

		}
	}

	// updating the k ghost cells
	for (int i = 2; i < Ni-2; ++i)
	{
		for (int j = 2; j < Nj-2; ++j)
		{
			//using cell (k=2) filling k=1
			double kFaceAreaVectorUnitVector[3] ={0.0,0.0,1.0};
			kFaceAreaVectorUnitVector[0] = kFaceAreaVector[i][j][2][0] / 
			sqrt( pow(kFaceAreaVector[i][j][2][0] ,2) + 
			pow(kFaceAreaVector[i][j][2][1] ,2) + 
			pow(kFaceAreaVector[i][j][2][2] ,2) ) ;  
			
			kFaceAreaVectorUnitVector[1] = kFaceAreaVector[i][j][2][1] / 
			sqrt( pow(kFaceAreaVector[i][j][2][0] ,2) + 
			pow(kFaceAreaVector[i][j][2][1] ,2) + 
			pow(kFaceAreaVector[i][j][2][2] ,2) ) ;
			
			kFaceAreaVectorUnitVector[2] = kFaceAreaVector[i][j][2][2] / 
			sqrt( pow(kFaceAreaVector[i][j][2][0] ,2) + 
			pow(kFaceAreaVector[i][j][2][1] ,2) + 
			pow(kFaceAreaVector[i][j][2][2] ,2) ) ;


			ConservedVariables[i][j][1][1] = (1-2*
			pow(kFaceAreaVectorUnitVector[0],2)) * 
			ConservedVariables[i][j][2][1] -(2*kFaceAreaVectorUnitVector[0]*
			kFaceAreaVectorUnitVector[1])*ConservedVariables[i][j][2][2] -
			(2*kFaceAreaVectorUnitVector[0]*kFaceAreaVectorUnitVector[2])*
			ConservedVariables[i][j][2][3] ;

			ConservedVariables[i][j][1][2] = -(2*
			kFaceAreaVectorUnitVector[1]*kFaceAreaVectorUnitVector[0])*
			ConservedVariables[i][j][2][1] + 
			(1-2*pow(kFaceAreaVectorUnitVector[1],2)) *
			ConservedVariables[i][j][2][2] - (2*kFaceAreaVectorUnitVector[1] * 
			kFaceAreaVectorUnitVector[2])*ConservedVariables[i][j][2][3] ;

			ConservedVariables[i][j][1][3] = -(2*kFaceAreaVectorUnitVector[2] *
			kFaceAreaVectorUnitVector[0])*ConservedVariables[i][j][2][1] -
			(2*kFaceAreaVectorUnitVector[2]*kFaceAreaVectorUnitVector[1])*
			ConservedVariables[i][j][2][2] + 
			(1-2*pow(kFaceAreaVectorUnitVector[2],2))* 
			ConservedVariables[i][j][2][3] ;

			// using zero order extrapolation
			ConservedVariables[i][j][1][0] = ConservedVariables[i][j][2][0] ;
			ConservedVariables[i][j][1][4] = ConservedVariables[i][j][2][4] ;


			//using cell (k=3) filling k=0
			kFaceAreaVectorUnitVector[0] = kFaceAreaVector[i][j][3][0] / 
			sqrt( pow(kFaceAreaVector[i][j][3][0] ,2) + 
			pow(kFaceAreaVector[i][j][3][1] ,2) + 
			pow(kFaceAreaVector[i][j][3][2] ,2) ) ;  
			
			kFaceAreaVectorUnitVector[1] = kFaceAreaVector[i][j][3][1] / 
			sqrt( pow(kFaceAreaVector[i][j][3][0] ,2) + 
			pow(kFaceAreaVector[i][j][3][1] ,2) + 
			pow(kFaceAreaVector[i][j][3][2] ,2) ) ;
			
			kFaceAreaVectorUnitVector[2] = kFaceAreaVector[i][j][3][2] / 
			sqrt( pow(kFaceAreaVector[i][j][3][0] ,2) + 
			pow(kFaceAreaVector[i][j][3][1] ,2) + 
			pow(kFaceAreaVector[i][j][3][2] ,2) ) ;


			ConservedVariables[i][j][0][1] = (1-2*
			pow(kFaceAreaVectorUnitVector [0],2)) * 
			ConservedVariables[i][j][3][1] -
			(2*kFaceAreaVectorUnitVector[0]*
			kFaceAreaVectorUnitVector[1])*ConservedVariables[i][j][3][2] -
			(2*kFaceAreaVectorUnitVector[0]*kFaceAreaVectorUnitVector[2])*
			ConservedVariables[i][j][3][3] ;

			ConservedVariables[i][j][0][2] = -(2*kFaceAreaVectorUnitVector[1]*
			kFaceAreaVectorUnitVector[0])*ConservedVariables[i][j][3][1] +
			(1-2*pow(kFaceAreaVectorUnitVector [1],2))*
			ConservedVariables[i][j][3][2] -(2*kFaceAreaVectorUnitVector[1]*
			kFaceAreaVectorUnitVector[2])*ConservedVariables[i][j][3][3] ;

			ConservedVariables[i][j][0][3] = -(2*kFaceAreaVectorUnitVector[2]*
			kFaceAreaVectorUnitVector[0])*ConservedVariables[i][j][3][1] -
			(2*kFaceAreaVectorUnitVector [2]*kFaceAreaVectorUnitVector [1])*
			ConservedVariables[i][j][3][2] + 
			(1-2*pow(kFaceAreaVectorUnitVector [2],2)) * 
			ConservedVariables[i][j][3][3] ;

			// using zero order extrapolation
			ConservedVariables[i][j][0][0] = ConservedVariables[i][j][3][0] ;
			ConservedVariables[i][j][0][4] = ConservedVariables[i][j][3][4] ;

			// using cell (k=Nk-3) filling k=Nk-2

			kFaceAreaVectorUnitVector[0] = - kFaceAreaVector[i][j][Nk-3][0] / 
			sqrt( pow(kFaceAreaVector[i][j][Nk-3][0] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-3][1] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-3][2] ,2) ) ;  
			
			kFaceAreaVectorUnitVector[1] = - kFaceAreaVector[i][j][Nk-3][1] / 
			sqrt( pow(kFaceAreaVector[i][j][Nk-3][0] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-3][1] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-3][2] ,2) ) ;
			
			kFaceAreaVectorUnitVector[2] = - kFaceAreaVector[i][j][Nk-3][2] / 
			sqrt( pow(kFaceAreaVector[i][j][Nk-3][0] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-3][1] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-3][2] ,2) ) ;
			

			ConservedVariables[i][j][Nk-2][1] = (1-2*
			pow(kFaceAreaVectorUnitVector[0],2)) * 
			ConservedVariables[i][j][Nk-3][1] -
			(2*kFaceAreaVectorUnitVector[0]*kFaceAreaVectorUnitVector[1])*
			ConservedVariables[i][j][Nk-3][2] -
			(2*kFaceAreaVectorUnitVector[0]*kFaceAreaVectorUnitVector[2])*
			ConservedVariables[i][j][Nk-3][3] ;

			ConservedVariables[i][j][Nk-2][2] = -(2*
			kFaceAreaVectorUnitVector[1]*kFaceAreaVectorUnitVector[0])*
			ConservedVariables[i][j][Nk-3][1] +
			(1-2*pow(kFaceAreaVectorUnitVector[1],2))*
			ConservedVariables[i][j][Nk-3][2] -
			(2*kFaceAreaVectorUnitVector[1]*kFaceAreaVectorUnitVector[2])*
			ConservedVariables[i][j][Nk-3][3] ;

			ConservedVariables[i][j][Nk-2][3] = -(2*
			kFaceAreaVectorUnitVector[2]*kFaceAreaVectorUnitVector[0])*
			ConservedVariables[i][j][Nk-3][1] -
			(2*kFaceAreaVectorUnitVector[2]*kFaceAreaVectorUnitVector[1])*
			ConservedVariables[i][j][Nk-3][2] +
			(1-2*pow(kFaceAreaVectorUnitVector[2],2))*
			ConservedVariables[i][j][Nk-3][3] ;

			// using zero order extrapolation
			ConservedVariables[i][j][Nk-2][0] = 
			ConservedVariables[i][j][Nk-3][0] ;
			ConservedVariables[i][j][Nk-2][4] = 
			ConservedVariables[i][j][Nk-3][4] ;

			//using cell(k=Nk-4) filling k=Nk-1
			kFaceAreaVectorUnitVector[0] = - kFaceAreaVector[i][j][Nk-4][0] / 
			sqrt( pow(kFaceAreaVector[i][j][Nk-4][0] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-4][1] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-4][2] ,2) ) ;  
			
			kFaceAreaVectorUnitVector[1] = - kFaceAreaVector[i][j][Nk-4][1] / 
			sqrt( pow(kFaceAreaVector[i][j][Nk-4][0] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-4][1] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-4][2] ,2) ) ;
			
			kFaceAreaVectorUnitVector[2] = - kFaceAreaVector[i][j][Nk-4][2] / 
			sqrt( pow(kFaceAreaVector[i][j][Nk-4][0] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-4][1] ,2) + 
			pow(kFaceAreaVector[i][j][Nk-4][2] ,2) ) ;
			
			ConservedVariables[i][j][Nk-1][1] = (1-2*
			pow(kFaceAreaVectorUnitVector[0],2)) * 
			ConservedVariables[i][j][Nk-4][1] -
			(2*kFaceAreaVectorUnitVector[0]*kFaceAreaVectorUnitVector[1])*
			ConservedVariables[i][j][Nk-4][2]-(2*kFaceAreaVectorUnitVector[0]*
			kFaceAreaVectorUnitVector[2])*ConservedVariables[i][j][Nk-4][3] ;

			ConservedVariables[i][j][Nk-1][2] = -(2*
			kFaceAreaVectorUnitVector[1]*kFaceAreaVectorUnitVector[0])*
			ConservedVariables[i][j][Nk-4][1] +
			(1-2*pow(kFaceAreaVectorUnitVector[1],2))*
			ConservedVariables[i][j][Nk-4][2] -
			(2*kFaceAreaVectorUnitVector[1]*kFaceAreaVectorUnitVector[2])*
			ConservedVariables[i][j][Nk-4][3] ;

			ConservedVariables[i][j][Nk-1][3] = -(2*
			kFaceAreaVectorUnitVector[2]*kFaceAreaVectorUnitVector[0])*
			ConservedVariables[i][j][Nk-4][1] -
			(2*kFaceAreaVectorUnitVector[2]*kFaceAreaVectorUnitVector[1])*
			ConservedVariables[i][j][Nk-4][2] +
			(1-2*pow(kFaceAreaVectorUnitVector[2],2))*
			ConservedVariables[i][j][Nk-4][3] ; 

			// using zero order extrapolation
			ConservedVariables[i][j][Nk-1][0] = 
			ConservedVariables[i][j][Nk-4][0] ;
			ConservedVariables[i][j][Nk-1][4] = 
			ConservedVariables[i][j][Nk-4][4] ;
		}
	}
}

