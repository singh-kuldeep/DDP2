
#include "math.h"
#include <iostream>
#include "math.h"
#include <fstream> /* For file handling */
#include <string> /* For strings */
#include <vector> /* For vectors*/
#include <cstdlib> 

// Function reads the data from previous solution file
#include "Reader.h" 

/*! \file Initial_condition.h
    \brief This header file all conserved parameters inside the domain(Only the
    live cells are initialized). 
    Currently, there are three different ways to do this. Use the appropriate 
    if statement.  
    \author Kuldeep Singh
    \date 2017
*/


/** \brief Function find_throat() Finds the location of the throat*/
double find_throat(
	std::vector<std::vector<double> > UpperWallCoordinates, 
	int & throat_location)
{
	throat_location = 0;
	double throat_area = 1e101; /** \param throat_area Area of the throat*/
	for (int i = 0; i < UpperWallCoordinates.size(); ++i)
	{
		if(throat_area>UpperWallCoordinates[i][1])
		{
			throat_area = UpperWallCoordinates[i][1];
			throat_location = i;
		}
	}
	return throat_area;
}

double getPressure(double InletPressure, double InletMach, double Mach, double SpecificHeatRatio)
{	
	double InletRP = 1 + ((SpecificHeatRatio-1)*pow(InletMach,2)/2); 
	double RP = 1 + ((SpecificHeatRatio-1)*pow(Mach,2)/2); 
	return InletPressure*pow((InletRP/RP),(SpecificHeatRatio)/(SpecificHeatRatio-1));
}

double getDensity(double InletDensity, double InletMach, double Mach, double SpecificHeatRatio)
{	
	double InletRP = 1 + ((SpecificHeatRatio-1)*pow(InletMach,2)/2); 
	double RP = 1 + ((SpecificHeatRatio-1)*pow(Mach,2)/2); 
	return InletDensity*pow((InletRP/RP),1/(SpecificHeatRatio-1));
}

/** \brief Function getMachDivergingDuct() Finds the local Mach number for a 
given area ratio in the diverging duct*/
double getMachDivergingDuct(double areaRatio)
{
	double MachDiverging ;
	//Diverging section 
	double mach_tolaerence = 1; // 0.001;
	double MachDivergingLower = 1;
	double MachDivergingUpper = 100;

	MachDiverging = 0.5*(MachDivergingLower + MachDivergingUpper); 
	/**\param MachDiverging This will always be above 1 so initially let's 
	keep it 100*/

	double gamma = 1.4;
	// calculating the diverging Mach first
	while(mach_tolaerence > 0.001)
	{
		MachDiverging = 0.5*(MachDivergingUpper + MachDivergingLower);
		double area_difference = areaRatio - pow((gamma+1)/2,-((gamma+1)/
			(2*(gamma-1))))*pow((1+0.5*(gamma-1)*MachDiverging*MachDiverging),
		((gamma+1)/(2*(gamma-1))))/MachDiverging;

		if(area_difference<0)
		{
			MachDivergingUpper = MachDiverging;
		}
		else if(area_difference>=0)
		{
			MachDivergingLower = MachDiverging;
		}
		mach_tolaerence = MachDivergingUpper - MachDivergingLower;		
	}
	return MachDiverging;
}


/** \brief Function getMachConvergingDuct() Finds the local Mach number for a 
given area ratio in the Converging duct*/
double getMachConvergingDuct(double areaRatio)
{
	double MachConverging;
	//Converging section 
	double mach_tolaerence = 1; // 0.001;
	double MachConvergingLower = 0;
	double MachConvergingUpper = 1;

	MachConverging = 0.5*(MachConvergingLower+MachConvergingUpper);
	/**\param MachConverging This will always be below 1 so initially let's 
	keep it 100*/
	double gamma = 1.4;
	// calculating the diverging Mach first
	while(mach_tolaerence > 0.001)
	{
		MachConverging = 0.5*(MachConvergingUpper + MachConvergingLower);

		double area_difference = areaRatio - pow((gamma+1)/2,-((gamma+1)/
			(2*(gamma-1))))*pow((1+0.5*(gamma-1)*MachConverging*MachConverging),
		((gamma+1)/(2*(gamma-1))))/MachConverging;

		if(area_difference>0)
		{
			MachConvergingUpper = MachConverging;
		}
		else if(area_difference<=0)
		{
			MachConvergingLower = MachConverging;
		}
		mach_tolaerence = MachConvergingUpper - MachConvergingLower;		
	}
	
	return MachConverging;
}


double getAreaRatio(double Mach)
{
	double areaRatio ;
	double gamma = 1.4;
	areaRatio = pow((gamma+1)/2,-((gamma+1)/(2*(gamma-1))))*
	pow((1+0.5*(gamma-1)*Mach*Mach),((gamma+1)/(2*(gamma-1))))/Mach;
	return areaRatio;
}


/** \brief Function InitialCondition()
*\param [in] ConservedVariables This is the pointer to the 4D vector where all 
the conserved variables of previous time step are stored.
*\param [in] ConservedVariablesNew This is the pointer to the 4D vector where 
all the conserved variables of next/new time step are stored.
*\param [in] Ni Number of cells in in "i" direction.  
*\param [in] Nj Number of cells in in "j" direction.  
*\param [in] Nk Number of cells in in "k" direction.  
*\return void
*/
void initial_condition(
	vector<vector<vector<vector<double> > > > & ConservedVariables,
	vector<vector<vector<vector<double> > > > & ConservedVariablesNew,
	int Ni, int Nj, int Nk, double SpecificHeatRatio)
{
	// if premitive variables are given
	double InletDensity ; 
	double InletXVelocity ; 
	double InletYVelocity ; 
	double InletZVelocity ; 
	double InletPressure ;
	double InletTemperature;
	double InletSoundSpeed;

	// if total quantities are given at the inlet 
	double InletTotalTemperature;	
	double InletTotalPressure;
	double InletMach;
	
	string InitialCondition ;
	string GeometryOption ;

	ifstream infile("inputfile");
	string aline;

	// reading the input file
	while(!infile.eof())// file ended
	{
		getline(infile,aline); // reading line form file

		if (aline.find( "//" )!=0 && aline.empty()==false) 
		{
			if (aline.find("InitialCondition")!=string::npos)
			{
				InitialCondition = aline.substr(aline.find("=")+2);
			}
			else if (aline.find("GeometryOption")!=string::npos)
			{
				GeometryOption = aline.substr(aline.find("=")+2);
			}
			// else if(aline.find("InletDensity")!=string::npos)
			// {
			// 	InletDensity = stod (aline.substr(aline.find("=")+1));
			// }
			// else if(aline.find("InletXVelocity")!=string::npos)
			// {
			// 	InletXVelocity = stod (aline.substr(aline.find("=")+1));
			// }
			// else if(aline.find("InletYVelocity")!=string::npos)
			// {
			// 	InletYVelocity = stod (aline.substr(aline.find("=")+1));
			// }
			// else if(aline.find("InletZVelocity")!=string::npos)
			// {
			// 	InletZVelocity = stod (aline.substr(aline.find("=")+1));
			// }
			// else if(aline.find("InletPressure")!=string::npos)
			// {
			// 	InletPressure = stod (aline.substr(aline.find("=")+1));
			// }
			else if(aline.find("InletTotalTemperature")!=string::npos)
			{
				InletTotalTemperature = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletTotalPressure")!=string::npos)
			{
				InletTotalPressure = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletMach")!=string::npos)
			{
				InletMach = stod (aline.substr(aline.find("=")+1));
			}								
		}
	}
	/** 
	(1) InitialCondition = ZeroVelocity
	Zero velocity inside the domain and total quantities are given at the inlet 
	*/

	/**
	(2) InitialCondition = FreeStreamParameterEverywhare
	Free stream parameters are filled inside the domain*/

	/**
	(3) InitialCondition = NozzleBasedOn1DCalculation
	Parameters are varying inside the nozzle domain based on the 
	1D is-entropic relation.This is specifically designed for the nozzle.
	*/	
	
	cout << "InitialCondition are " << InitialCondition <<  endl;

	InletPressure = InletTotalPressure*
	pow((1+((SpecificHeatRatio-1)*pow(InletMach,2)/2)),(-SpecificHeatRatio/(SpecificHeatRatio-1)));
	
	InletTemperature = InletTotalTemperature/(1+((SpecificHeatRatio-1)*pow(InletMach,2)/2)); 
	
	InletDensity = InletPressure/(InletTemperature*RealGasConstant);

	InletSoundSpeed = sqrt(SpecificHeatRatio*InletPressure/InletDensity);
	InletXVelocity = InletSoundSpeed*InletMach;

	// cout << InletSoundSpeed << "," << InletPressure <<","<< InletTemperature << endl;
	if(InitialCondition =="ZeroVelocity")
	{
		for (int i =0; i < Ni; ++i)
		{
			for (int j =0; j < Nj; ++j)
			{
				for (int k =0; k < Nk; ++k)
				{
					ConservedVariables[i][j][k][0] = InletDensity;
					ConservedVariables[i][j][k][1] = 0; 
					ConservedVariables[i][j][k][2] = 0; 
					ConservedVariables[i][j][k][3] = 0; 
					ConservedVariables[i][j][k][4] = InletPressure*
					(SpecificHeatRatio-1) ;

					ConservedVariablesNew[i][j][k][0] = InletDensity;
					ConservedVariablesNew[i][j][k][1] = 0; 
					ConservedVariablesNew[i][j][k][2] = 0; 
					ConservedVariablesNew[i][j][k][3] = 0; 
					ConservedVariablesNew[i][j][k][4] = InletPressure*
					(SpecificHeatRatio-1) ;
				}
			}
		}
	}
	
	if (InitialCondition == "FreeStreamParameterEverywhare")
	{
		/*Sometime free stream conditions are filled in the domain*/
		/**
	     * Inlet conditions are user given data.
	     * one has to mention the free stream parameters at inlet (ex. static 
	     pressure (\f$ P_inf \f$), temperature(\f$ T \f$))
	     */
		for (int i =0; i < Ni; ++i)
		{
			for (int j =0; j < Nj; ++j)
			{
				for (int k =0; k < Nk; ++k)
				{
					ConservedVariables[i][j][k][0] = InletDensity;
					ConservedVariables[i][j][k][1] = InletDensity*InletXVelocity ;
					ConservedVariables[i][j][k][2] = 0 ;
					ConservedVariables[i][j][k][3] = 0 ;
					ConservedVariables[i][j][k][4] = (SpecificHeatRatio-1)*
					InletPressure+0.5*InletDensity*pow(InletXVelocity,2);

					ConservedVariablesNew[i][j][k][0] = InletDensity;
					ConservedVariablesNew[i][j][k][1] = InletDensity*InletXVelocity ;
					ConservedVariablesNew[i][j][k][2] = 0 ;
					ConservedVariablesNew[i][j][k][3] = 0 ;
					ConservedVariablesNew[i][j][k][4] = (SpecificHeatRatio-1)*
					InletPressure+0.5*InletDensity*pow(InletXVelocity,2);
				}
			}
		}
	}
	if (InitialCondition == "FreeStreamParameterAndZeroVelocity")
	{
		/*Sometime free stream conditions are filled in the domain*/
		/**
	     * Inlet conditions are user given data.
	     * one has to mention the free stream parameters at inlet (ex. static 
	     pressure (\f$ P_inf \f$), temperature(\f$ T \f$))
	     */
		for (int i =0; i < Ni; ++i)
		{
			for (int j =0; j < Nj; ++j)
			{
				for (int k =0; k < Nk; ++k)
				{
					ConservedVariables[i][j][k][0] = InletDensity;
					ConservedVariables[i][j][k][1] = 0 ;
					ConservedVariables[i][j][k][2] = 0 ;
					ConservedVariables[i][j][k][3] = 0 ;
					ConservedVariables[i][j][k][4] = (SpecificHeatRatio-1)*
					InletPressure+0.5*InletDensity*pow(InletXVelocity,2);

					ConservedVariablesNew[i][j][k][0] = InletDensity;
					ConservedVariablesNew[i][j][k][1] = 0 ;
					ConservedVariablesNew[i][j][k][2] = 0 ;
					ConservedVariablesNew[i][j][k][3] = 0 ;
					ConservedVariablesNew[i][j][k][4] = (SpecificHeatRatio-1)*
					InletPressure+0.5*InletDensity*pow(InletXVelocity,2);
				}
			}
		}
	}
	
	// Start from the previous solution
	if (InitialCondition == "StartFromPreviousSolution")
	{
		ReadConservedVariables(ConservedVariables, ConservedVariablesNew,
		Ni, Nj, Nk);
	}



	/* By calculating the mach number using the area ratio then filling the 
	Conserved Variables which are more close to the solution. This further 
	helps in fast convergance. */
	#if 1
	if(InitialCondition == "NozzleBasedOn1DCalculation") 
	{
		// Reading the nozzle geometry
		std::vector<std::vector<double> > UpperWallCoordinates;
		std::vector<std::vector<double> > LowerWallCoordinates;

		ifstream xup("./NozzleGeomatryGenrator/XCoordinatesUpperWall.csv");
		ifstream yup("./NozzleGeomatryGenrator/YCoordinatesUpperWall.csv");

		int pointCouter = 0;
	   	while(!xup.eof())
	   	{
		   string aline;
		   double xt,yt;

		   getline(xup,aline);
		   xt = atof( aline.c_str() );
		   
		   getline(yup,aline);
		   yt = atof( aline.c_str() );
		   
		   vector<double> temp;
		   temp.push_back(xt); 
		   temp.push_back(yt);
		   UpperWallCoordinates.push_back(temp);

		   temp[1] = 0.0; // change the y only and push it to the Down vector
		   LowerWallCoordinates.push_back(temp);
		   ++pointCouter;
	   	}
		// randomly extra zeros at the end so to remove them pop is used
		UpperWallCoordinates.pop_back();
		LowerWallCoordinates.pop_back();

		// closing the file
		xup.close();
		yup.close();

		/** \parm throat_location Location of the throat */
		int throat_location = 10; 

		/** \param throat_area Area at the throat*/
		double throat_area = find_throat(UpperWallCoordinates,throat_location);

	  	/**\param Density  Density at local point */
		/**\param Pressure  Pressure at local point */
		/**\param Temperature  Temperature at local point */
		/**\param VelocityMagnitude  Velocity magnitude at local point */
		/** \param Mach Mach number at a location */
		/** \param AreaRatio Area ratio at a location with throat area*/
		double AreaRatio, Mach, Pressure, Density, Temperature, VelocityMagnitude;

		/* Inlet velocity and density are taken form the input file
		1. Inlet mach is calculated from geomatery
		2. Then Inlet pressure is calculated p = \f$ (rho*V^2)/(M^2 * gamma)
		3. No density and pressure are calculated at every other point*/
		double InletAreaRatio = (UpperWallCoordinates[0][1]/throat_area);
		// double InletMach = getMachConvergingDuct(InletAreaRatio);

		// Subsonic Initial condition in the converging part
		for (int i = 0; i <=throat_location; ++i)
		{
			for (int j = 0; j < Nj; ++j)
			{
				for (int k = 0; k < Nk; ++k)
				{
					AreaRatio = UpperWallCoordinates[i][1]/throat_area;

					Mach = getMachConvergingDuct(AreaRatio);
					
					Pressure = InletTotalPressure/pow((1+((SpecificHeatRatio-1)*
					pow(Mach,2)/2)),(SpecificHeatRatio/(SpecificHeatRatio-1)));
					
					Temperature = InletTotalTemperature/(1+((SpecificHeatRatio-1)*pow(Mach,2)/2));

					Density = Pressure/(RealGasConstant*Temperature);
					/**\param VelocityMagnitude Velocity magnitude at any x(i) 
					location */
					VelocityMagnitude = Mach * 
					sqrt(SpecificHeatRatio*Pressure/Density);
					
					ConservedVariables[i][j][k][0] = Density;
					ConservedVariables[i][j][k][1] = Density*VelocityMagnitude ;
					ConservedVariables[i][j][k][2] = 0;
					ConservedVariables[i][j][k][3] = 0;
					ConservedVariables[i][j][k][4] = 
					Pressure/(SpecificHeatRatio-1)+ 0.5*Density*
					VelocityMagnitude*VelocityMagnitude ;

					ConservedVariablesNew[i][j][k][0] = Density;
					ConservedVariablesNew[i][j][k][1] = Density*VelocityMagnitude;
					ConservedVariablesNew[i][j][k][2] = 0;
					ConservedVariablesNew[i][j][k][3] = 0;
					ConservedVariablesNew[i][j][k][4] = 
					Pressure/(SpecificHeatRatio-1) + 0.5*Density*
					VelocityMagnitude*VelocityMagnitude ;
				}
			}
		}

		// Only the diverging section is initialized here
		// InletXVelocity = 0*InletXVelocity/MachFreestream;
		for (int i = throat_location+1; i < Ni; ++i) 
		{
			for (int j = 0; j < Nj; ++j)
			{
				for (int k = 0; k < Nk; ++k)
				{
					// cout << i << ", " << j <<", " << k << "  " << "Mach " << Mach << endl;
					// cout << "OKOKOKOK" << endl;
					AreaRatio = UpperWallCoordinates[i][1]/throat_area;

					Mach = getMachDivergingDuct(AreaRatio);
					
					Pressure = InletTotalPressure/pow((1+((SpecificHeatRatio-1)*
					pow(Mach,2)/2)),(SpecificHeatRatio/(SpecificHeatRatio-1)));
					
					Temperature = InletTotalTemperature/(1+((SpecificHeatRatio-1)*pow(Mach,2)/2));

					Density = Pressure/(RealGasConstant*Temperature);
					
					/**\param VelocityMagnitude Velocity magnitude at x(i) 
					location */
					VelocityMagnitude = Mach * 
					sqrt(SpecificHeatRatio*Pressure/Density);

					ConservedVariables[i][j][k][0] = Density;
					ConservedVariables[i][j][k][1] = Density*VelocityMagnitude ;
					ConservedVariables[i][j][k][2] = 0;
					ConservedVariables[i][j][k][3] = 0;
					ConservedVariables[i][j][k][4] = 
					Pressure/(SpecificHeatRatio-1) + 0.5*Density*
					VelocityMagnitude*VelocityMagnitude ;

					ConservedVariablesNew[i][j][k][0] = Density;
					ConservedVariablesNew[i][j][k][1] = Density*VelocityMagnitude;
					ConservedVariablesNew[i][j][k][2] = 0;
					ConservedVariablesNew[i][j][k][3] = 0;
					ConservedVariablesNew[i][j][k][4] = 
					Pressure/(SpecificHeatRatio-1)+ 0.5*Density*
					VelocityMagnitude*VelocityMagnitude ;
				}
			}
		}
	}
	#endif 
}
	