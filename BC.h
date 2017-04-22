/*! \file BC.h
    \brief This header file implements all three boundary conditions. 
    - Inlet
    - Exit and 
    - Wall boundary

    \author Kuldeep Singh
    \date 2017 
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
	vector<vector<vector<vector<double> > > > & iFaceAreaVector,
	vector<vector<vector<vector<double> > > > & jFaceAreaVector,
	vector<vector<vector<vector<double> > > > & kFaceAreaVector,
	int Ni, int Nj, int Nk, string InletBC, double TemperatureFreestream,
	double PressureFreestream, double MachFreestream )
{
	// double TemperatureFreestream = 300.00;//5180.76 ; 
	// *\param TemperatureFreestream TemperatureFreestream at inlet   
	// double PressureFreestream = 1e5;//7927660.8; 
	/**\param PressureFreestream  Inletpressure at inlet */
	double DensityFreestream = PressureFreestream /
		(IdealGasConstant*TemperatureFreestream) ; 
	/**\param DensityFreestream DensityFreestream at inlet */
	// double MachFreestream = 0.40918; // For inlet area ratio 1.5585
	/**\param MachFreestream MachFreestream at inlet */

	/**\param VelocityFreestream VelocityFreestream at inlet */
	double VelocityFreestream = MachFreestream*sqrt(SpecificHeatRatio*
		IdealGasConstant*TemperatureFreestream);


	
	/* Inlet ghost cells are being updated using the stagnation quantities
	(\f$ P_0, T_0 \f$)*/
	if(InletBC=="ImposeStagnationParameters")
	{
		/**\param TemperatureStagnation Stagnation temperature at inlet */  
		double TemperatureStagnation = TemperatureFreestream*(1+
		(SpecificHeatRatio-1)*(MachFreestream*MachFreestream)/2);//5180.76 ; 
		
		/**\param PressureStagnation Stagnation pressure at inlet */
		double PressureStagnation = PressureFreestream*pow((1+(SpecificHeatRatio-1)*
		(MachFreestream*MachFreestream)/2),(SpecificHeatRatio/(SpecificHeatRatio-1)));//7927660.8; 
		
		/**\param DensityStagnation Stagnation density at inlet */
		double DensityStagnation = PressureStagnation /
			(IdealGasConstant*TemperatureStagnation) ; 
			

		/* Inlet ghost cells are being updated using the stagnation quantities
		(\f$ P_0, T_0 \f$) and flow direction */
		for (int j =2; j < Nj-2; ++j)
		{
			for (int k =2; k < Nk-2; ++k)
			{
				/**\param PressureInlet Static pressure at inlet */
				double PressureInlet = 
				(SpecificHeatRatio-1)*(ConservedVariables[2][j][k][4] - 0.5*(
					pow(ConservedVariables[2][j][k][1],2) +
					pow(ConservedVariables[2][j][k][2],2) + 
					pow(ConservedVariables[2][j][k][3],2)) /
					ConservedVariables[2][j][k][0]) ;

				/**\param MachInlet MachInlet number at inlet */
				double MachInlet=sqrt((2/(SpecificHeatRatio-1))*(pow((PressureStagnation/
				PressureInlet),((SpecificHeatRatio-1)/SpecificHeatRatio) ) -1 ) ) ;
				/**\param TemperatureInlet Static temperature at inlet */
				double TemperatureInlet = TemperatureStagnation/(1+
				((SpecificHeatRatio-1)*MachInlet*MachInlet)/2);
				
				/**\param VelocityInlet Flow velocity at inlet */
				double VelocityInlet = MachInlet * sqrt(SpecificHeatRatio*
				IdealGasConstant*TemperatureInlet) ;
				
				/**\param DensityInlet Flow density at inlet */
				double DensityInlet = DensityStagnation / 
				pow((PressureStagnation/PressureInlet),(1/SpecificHeatRatio)) ;

				ConservedVariables[0][j][k][0] = DensityInlet ;
				ConservedVariables[0][j][k][1] = DensityInlet*VelocityInlet ;
				ConservedVariables[0][j][k][2] = 0 ;
				ConservedVariables[0][j][k][3] = 0 ;
				ConservedVariables[0][j][k][4] = PressureInlet/(SpecificHeatRatio-1)
				+0.5*DensityInlet*VelocityInlet*VelocityInlet ;

				ConservedVariables[1][j][k][0] = DensityInlet ;
				ConservedVariables[1][j][k][1] = DensityInlet*VelocityInlet;
				ConservedVariables[1][j][k][2] = 0 ;
				ConservedVariables[1][j][k][3] = 0 ;
				ConservedVariables[1][j][k][4] = PressureInlet/(SpecificHeatRatio-1)
				+0.5*DensityInlet*VelocityInlet*VelocityInlet ;
			}
		}
	}

	/*Implementing the non-reflecting boundary condition at the inlet for 
	subsonic flow, where DeltaW1 and DeltaW2 are zero*/
	/* At inlet updating the i ghost cells(i=0, i=1),(this is true where flow is 
	subsonic)*/
	if(InletBC=="NRBC")
	{	
		/**\param Density  Density at local point */
		/**\param Pressure  Pressure at local point */
		/**\param Temperature  Temperature at local point */
		double Density ;
		double Velocity ;	
		double Pressure ;
		double SoundSpeed ;
		
		//change in the parameters
		double DeltaW ;  
		double DeltaDensity ;
		double DeltaVelocity ;	
		double DeltaPressure ;
		for (int j =2; j < Nj-2; ++j)
		{
			for (int k =2; k < Nk-2; ++k)
			{
				//Using second live cell(i=3) filling the first ghost cell (i=0)
				Density = ConservedVariables[3][j][k][0] ;
				Velocity = ConservedVariables[3][j][k][1]/ConservedVariables[3][j][k][0];

				Pressure = (ConservedVariables[3][j][k][4] - 0.5*Density*Velocity*
					Velocity)*(SpecificHeatRatio-1); 
				SoundSpeed = sqrt(SpecificHeatRatio*Pressure/Density);

				DeltaW = (Velocity-VelocityFreestream)-((Pressure-PressureFreestream)/(Density*SoundSpeed));

				DeltaPressure = -Density*SoundSpeed*DeltaW/2;
				DeltaDensity = DeltaPressure/(SoundSpeed*SoundSpeed);
				DeltaVelocity = -DeltaPressure/(Density*SoundSpeed);

				ConservedVariables[0][j][k][0] = Density+DeltaDensity;
				ConservedVariables[0][j][k][1] = ConservedVariables[0][j][k][0]*(Velocity+DeltaVelocity);
				ConservedVariables[0][j][k][2] = 0;
				ConservedVariables[0][j][k][3] = 0;
				ConservedVariables[0][j][k][4] = (Pressure+DeltaPressure)/(SpecificHeatRatio-1)
				+0.5*(Density+DeltaDensity)*(Velocity+DeltaVelocity)*(Velocity+DeltaVelocity);

				//Using second live cell(i=2) filling the second ghost cell (i=1)
				Density = ConservedVariables[2][j][k][0] ;
				Velocity = ConservedVariables[2][j][k][1]/ConservedVariables[2][j][k][0];

				Pressure = (ConservedVariables[2][j][k][4] - 0.5*Density*Velocity*
					Velocity)*(SpecificHeatRatio-1); 
				SoundSpeed = sqrt(SpecificHeatRatio*Pressure/Density);

				DeltaW = (Velocity-VelocityFreestream)-((Pressure-PressureFreestream)/(Density*SoundSpeed));

				DeltaPressure = -Density*SoundSpeed*DeltaW/2;
				DeltaDensity = DeltaPressure/(SoundSpeed*SoundSpeed);
				DeltaVelocity = -DeltaPressure/(Density*SoundSpeed);

				ConservedVariables[1][j][k][0] = Density+DeltaDensity;
				ConservedVariables[1][j][k][1] = ConservedVariables[1][j][k][0]*(Velocity+DeltaVelocity);
				ConservedVariables[1][j][k][2] = 0;
				ConservedVariables[1][j][k][3] = 0;
				ConservedVariables[1][j][k][4] = (Pressure+DeltaPressure)/(SpecificHeatRatio-1)
				+0.5*(Density+DeltaDensity)*(Velocity+DeltaVelocity)*(Velocity+DeltaVelocity);
			}
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
			double jFaceAreaVector_unit_vector[3] = {0.0,1.0,0.0} ;

			//using leftcell (j=2) filling j=1
			//calculating the normal vector
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

