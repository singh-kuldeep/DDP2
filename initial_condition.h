#include "math.h"
#include <iostream>
#include "math.h"
#include <fstream> /* For file handling */
#include <string> /* For strings */
#include <vector> /* For vectors*/
#include <cstdlib> 

/*! \file Initial_condition.h
    \brief This header file all conserved parameters inside the domain(Only the
    live cells are initialized). 
    Currently, there are three different ways to do this. Use the appropriate 
    if statement.  
    \author Kuldeep Singh
    \date 2017
*/
#define SpecificHeatRatio 1.4 /*!< This is gas constant (Gamma). For air at 
room temperature it is almost equal to 1.4. If you are using some other gas at
some other temperature then change it*/
#define IdealGasConstant 287.14 
/*!< This is ideal gas constant \f$ R(J Kg^{-1}K^{-1}) = (c_p - c_v)\f$ */


/** \brief Function find_throat() Finds the location of the throat*/
double find_throat_area(
	std::vector<std::vector<double> > UpperWallCoordinates, 
	int & throat_location)
{
	throat_location = 0;
	double throat_area = 10000; /** \param throat_area Area of the throat*/
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
	int Ni, int Nj, int Nk)
{
	double InletDensity ; 
	double InletXVelocity ; 
	double InletYVelocity ; 
	double InletZVelocity ; 
	double InletStaticPressure ;
	
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
			else if(aline.find("InletDensity")!=string::npos)
			{
				InletDensity = stoi (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletXVelocity")!=string::npos)
			{
				InletXVelocity = stoi (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletYVelocity")!=string::npos)
			{
				InletYVelocity = stoi (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletZVelocity")!=string::npos)
			{
				InletZVelocity = stoi (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletStaticPressure")!=string::npos)
			{
				InletStaticPressure = stoi (aline.substr(aline.find("=")+1));
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
					ConservedVariables[i][j][k][4] = InletStaticPressure*
					(SpecificHeatRatio-1) ;

					ConservedVariablesNew[i][j][k][0] = InletDensity;
					ConservedVariablesNew[i][j][k][1] = 0; 
					ConservedVariablesNew[i][j][k][2] = 0; 
					ConservedVariablesNew[i][j][k][3] = 0; 
					ConservedVariablesNew[i][j][k][4] = InletStaticPressure*
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
					InletStaticPressure+0.5*InletDensity*pow(InletXVelocity,2);

					ConservedVariablesNew[i][j][k][0] = InletDensity;
					ConservedVariablesNew[i][j][k][1] = InletDensity*InletXVelocity ;
					ConservedVariablesNew[i][j][k][2] = 0 ;
					ConservedVariablesNew[i][j][k][3] = 0 ;
					ConservedVariablesNew[i][j][k][4] = (SpecificHeatRatio-1)*
					InletStaticPressure+0.5*InletDensity*pow(InletXVelocity,2);
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
					InletStaticPressure+0.5*InletDensity*pow(InletXVelocity,2);

					ConservedVariablesNew[i][j][k][0] = InletDensity;
					ConservedVariablesNew[i][j][k][1] = 0 ;
					ConservedVariablesNew[i][j][k][2] = 0 ;
					ConservedVariablesNew[i][j][k][3] = 0 ;
					ConservedVariablesNew[i][j][k][4] = (SpecificHeatRatio-1)*
					InletStaticPressure+0.5*InletDensity*pow(InletXVelocity,2);
				}
			}
		}
	}
	#if 0	
	// For the CD Nozzle (to make it converge faster)
	if(InitialCondition == "NozzleBasedOn1DCalculation") 
	{
		/**\param Density  Density at local point */
		/**\param Pressure  Pressure at local point */
		/**\param Temperature  Temperature at local point */
		/**\param VelocityMagnitude  Velocity magnitude at local point */
		double Density, Pressure, Temperature, VelocityMagnitude;

		/** \param Mach Mach number at a location */
		double Mach = 1.0; 
		
		/** \param localAreaRatio Area ratio at a location with throat area*/
		double localAreaRatio = 1.0;

		/**\param ratio_parameter \f$ 1 + ((Gamma-1)*M^2)/2 . Just to simplify 
		the	calculations */
		double ratio_parameter;	

		/* By calculating the area ratio and mach number then filling the 
		Conserved Variables which are more close to the solution. This further 
		helps in fast convergance. */

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
		double throat_area = find_throat_area(UpperWallCoordinates,throat_location);

		
		// Subsonic Initial condition in the converging part
		for (int i = 2; i <=throat_location+1; ++i)
		{
			for (int j = 0; j < Nj; ++j)
			{
				for (int k = 0; k < Nk; ++k)
				{
					localAreaRatio = UpperWallCoordinates[(i-2)][1]/throat_area;
					Mach = getMachConvergingDuct(localAreaRatio);
					ratio_parameter = 1 + (SpecificHeatRatio-1)*0.5*Mach*Mach;
					
					/**\param Temperature Temperature at x(i) location */
					Temperature = TemperatureStagnation/ratio_parameter;
					
					/**\param Pressure Static pressure at x(i) location */
					Pressure = PressureStagnation*pow(ratio_parameter,
						(-SpecificHeatRatio)/(SpecificHeatRatio-1));
					
					Density = Pressure/(IdealGasConstant*Temperature);

					/**\param VelocityMagnitude Velocity magnitude at x(i) 
					location */
					VelocityMagnitude = Mach * 
					sqrt(SpecificHeatRatio*IdealGasConstant*Temperature);
					
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
		for (int i = throat_location+2; i < Ni-2; ++i) 
		{
			for (int j =0; j < Nj; ++j)
			{
				for (int k =0; k < Nk; ++k)
				{
					localAreaRatio = UpperWallCoordinates[(i-2)][1]/throat_area;
					Mach = getMachDivergingDuct(localAreaRatio);
					// Mach = getMachConvergingDuct(localAreaRatio);
					ratio_parameter = 1 + (SpecificHeatRatio-1)*0.5*Mach*Mach;
					
					/**\param Temperature Temperature at x(i) location */
					Temperature = TemperatureStagnation/ratio_parameter;
					
					/**\param Pressure Static pressure at x(i) location */
					Pressure = PressureStagnation*pow(ratio_parameter,
						(-SpecificHeatRatio)/(SpecificHeatRatio-1));
					
					Density = Pressure/(IdealGasConstant*Temperature);

					/**\param VelocityMagnitude Velocity magnitude at x(i) 
					location */
					VelocityMagnitude = Mach * 
					sqrt(SpecificHeatRatio*IdealGasConstant*Temperature);
					
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
	#if 0 // for the testing purposes
	// storing the all conserved variables in one plane just after initialization
	ofstream kullu_2D_initial ;
	kullu_2D_initial.open("2D_parameters_B.csv");
	// kullu_2D_initial << "density" << "," << "density*u" << ","<< "density*v"
	// << "," << "density*w" << "," << "energy"  << endl ;
	for (int i = 2; i < Ni-2; ++i)
	{
		for (int j = 2; j < Nj-2; ++j)
		{
			kullu_2D_initial << ConservedVariables[i][j][Nk/2][0] << "," << 
			ConservedVariables[i][j][Nk/2][1] <<","<< 
			ConservedVariables[i][j][Nk/2][2] << "," <<
			ConservedVariables[i][j][Nk/2][3] << "," <<
			ConservedVariables[i][j][Nk/2][4] << endl ;
		}
	}
	// return 0;
	#endif
	#if 0
	/**@bug Every time simulation starts from first iteration. So, to save the
	simulation it is good to start from the last solution as the initial 
	condition*/
		cout << "Do you wants to start the simulation,where you left?(enter y)"
		<< endl << "Other wise press any key"<< endl;
		char oldstart;
		cin>> oldstart;
		if (oldstart=='y')
		{
			ifstream nozzleData ("./Results/outputfiles/restart.csv");
			string aline;

			// Initial condition are taken from the previously solved part
			for (int i = 0; i < Ni; ++i)
			{
				for (int j = 0; j < Nj; ++j)
				{
					for (int k = 0; k < Nk; ++k)
					{
						for (int l = 0; l < 5; ++l)
						{
							getline(nozzleData,aline);
							ConservedVariables[i][j][k][l]=atof(aline.c_str());
							ConservedVariablesNew[i][j][k][l] = 
							ConservedVariables[i][j][k][l];
						}
					}
				}
			}	
		}
	#endif
}
	