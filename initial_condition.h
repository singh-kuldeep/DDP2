#include <fstream>
#include "math.h"
#include <iostream>
#include "math.h"
#include <fstream> /* For file handling */
#include <string> /* For strings */
#include <vector> /* For vectors*/
#include <cstdlib> 

/*! \file initial_condition.h
    \brief This header file all conserved parameters inside the domain. 
    Cureently, there are three different ways to do this. Use the appropriate 
    switch statement.  
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
	std::vector<std::vector<double> > UpperCoordinates, 
	int & throat_location)
{
	throat_location = 0;
	double throat_area = 10000; /** \param throat_area Area of the throat*/
	for (int i = 0; i < UpperCoordinates.size(); ++i)
	{
		if(throat_area>UpperCoordinates[i][1])
		{
			throat_area = UpperCoordinates[i][1];
			throat_location = i;
		}
	}
	return throat_area;
}

/** \brief Function getMach() Finds the local Mach number for a given area 
ratio*/
double getMach(double areaRatio)
{
	double mach_tolaerence = 1 ; // 0.001;
	double Mach_down = 1.0001;
	double Mach_up = 100;
	double Mach = 0.5*(Mach_up + Mach_down);
	double gamma = 1.4;
	while(mach_tolaerence > 0.001)
	{
		Mach = 0.5*(Mach_up + Mach_down);
		double area_difference = areaRatio - pow((gamma+1)/2,-((gamma+1)/
			(2*(gamma-1))))*pow((1+0.5*(gamma-1)*Mach*Mach),((gamma+1)/
			(2*(gamma-1))))/Mach;

		if(area_difference<0)
		{
			Mach_up = Mach;
		}
		else if(area_difference>=0)
		{
			Mach_down = Mach;
		}
		mach_tolaerence = Mach_up - Mach_down;		
	}
	return Mach;
}

/** \brief Function getAreaRatio() Finds the area ratio for a given local Mach 
number*/
double getAreaRatio(double Mach)
{
	double areaRatio ;
	double gamma = 1.4;
	areaRatio = pow((gamma+1)/2,-((gamma+1)/(2*(gamma-1))))*pow((1+0.5*
		(gamma-1)*Mach*Mach),((gamma+1)/(2*(gamma-1))))/Mach;
	return areaRatio;
}


/** \brief Function initial_condition()
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
	int initial_condition = 1;// Zero velocity
	// int initial_condition = 2;// Nozzle varing mach number
	// int initial_condition = 3;// Start from the previous solution
	switch (initial_condition)
	{
		// Case (1) putting zero velocity every where
		case 1 :
			for (int i =0; i < Ni; ++i)
			{
				for (int j =0; j < Nj; ++j)
				{
					for (int k =0; k < Nk; ++k)
					{
						ConservedVariables[i][j][k][0] = 1.16;
						ConservedVariables[i][j][k][1] = 0 ;
						ConservedVariables[i][j][k][2] = 0 ;
						ConservedVariables[i][j][k][3] = 0 ;
						ConservedVariables[i][j][k][4] = 271767;

						ConservedVariablesNew[i][j][k][0] = 1.16;
						ConservedVariablesNew[i][j][k][1] = 0 ;
						ConservedVariablesNew[i][j][k][2] = 0 ;
						ConservedVariablesNew[i][j][k][3] = 0 ;
						ConservedVariablesNew[i][j][k][4] = 271767;
					}
				}
			}
		break;
		
		case 2 :
		/* By calculating the area ratio and mach number then filling the 
		Conserved Variables which are more close to the solution. This further 
		helps in fast convergance. */

			// Reading the nozzle geometry
			std::vector<std::vector<double> > UpperCoordinates;
			std::vector<std::vector<double> > DownCoordinates;

			ifstream xup("/home/kullu/Desktop/Acad/SEM10/DDP2/Code_DDP2/DDP2/NozzleGeomatryGenrator/XCoordinatesUpperWall.csv");
			ifstream yup("/home/kullu/Desktop/Acad/SEM10/DDP2/Code_DDP2/DDP2/NozzleGeomatryGenrator/YCoordinatesUpperWall.csv");

			int j = 0;
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
			   UpperCoordinates.push_back(temp);

			   temp[1] = 0.0; // change the y only and push it to the Down vector
			   DownCoordinates.push_back(temp);
			   ++j;
		   	}
		   // rendomaly extra zeros at the end so to remove them pop is used
		   UpperCoordinates.pop_back();
		   DownCoordinates.pop_back();

		   // closing the file
		   xup.close();
		   yup.close();

			/** \parm throat_location Location of the throat */
			int throat_location = 0; 

			/** \param throat_area Area at the throat*/
			double throat_area = find_throat_area(UpperCoordinates,throat_location);


			/** \param local_area_ratio Area ratio at a location with throat area*/
			double local_area_ratio = 1.0;
			
			/** \param Mach Mach number at a location */
			double Mach = 2.0; 

			/**
		     * Inlet conditions are user given data.
		     * one has to mention the stagnation parameters at inlet (ex. stagnation 
		     pressure (\f$ P_0 \f$), temperature(\f$ T_0 \f$))
		     */

			/**\param TemperatureStagnation Stagnation temperature at inlet */  
			double TemperatureStagnation = 5180.76 ; 
			
			/**\param PressureStagnation Stagnation pressure at inlet */
			double PressureStagnation = 7927660.8; 
			
			/**\param DensityStagnation Stagnation density at inlet */
			double DensityStagnation = PressureStagnation /
				(IdealGasConstant*TemperatureStagnation) ; 

			double Density, Pressure, Temperature, VelocityMagnitude;

			/**\param ratio_parameter \f$ 1 + ((Gamma-1)*M^2)/2 . Just to simplify the 
			calulations */
			double ratio_parameter;	

			/**\param theta Geometry rotation angle about the k (or z axis which is 
			passing through origin) direction */
			double theta = 3.14159265 * 0 / 180 ;
			
			/**\param FlowAngle Flow angle at upper wall*/
			double FlowAngle = 3.14159265 * 10 / 180 ; // upper wall is at 10 degree

			// Only the diverging section is initialised here
			for (int i =throat_location+1; i < Ni; ++i) 
			{
				for (int j =0; j < Nj; ++j)
				{
					for (int k =0; k < Nk; ++k)
					{
						FlowAngle = atan((UpperCoordinates[i][1]-
							UpperCoordinates[i-1][1])/
						(UpperCoordinates[i][0]-UpperCoordinates[i-1][0]));
						
						Mach = getMach(UpperCoordinates[i][1]/throat_area);

						ratio_parameter = 1 + (SpecificHeatRatio-1)*0.5*Mach*Mach;
						
						/**\param Temperature Temperature at x(i) location */
						Temperature = TemperatureStagnation/ratio_parameter;
						
						/**\param Pressure Static pressure at x(i) location */
						Pressure = PressureStagnation*
						pow(ratio_parameter,(-SpecificHeatRatio)/(SpecificHeatRatio-1));
						
						Density = Pressure/(IdealGasConstant*Temperature);

						/**\param VelocityMagnitude Velocity magnitude at x(i) location */
						VelocityMagnitude = Mach * 
						sqrt(Density*IdealGasConstant*Temperature);
						
						
						ConservedVariables[i][j][k][0] = Density;
						ConservedVariables[i][j][k][1] = cos(FlowAngle*(j-2)/(Nj-5))*
						cos(theta)*Density*VelocityMagnitude ;
						ConservedVariables[i][j][k][2] = -sin(FlowAngle*(j-2)/(Nj-5))*
						cos(theta)*Density*VelocityMagnitude ;
						ConservedVariables[i][j][k][3] = -sin(theta)*Density*
						VelocityMagnitude ;
						ConservedVariables[i][j][k][4] = Pressure/(SpecificHeatRatio-1)
						+0.5*Density*VelocityMagnitude*VelocityMagnitude ;

						ConservedVariablesNew[i][j][k][0] = Density ;
						ConservedVariablesNew[i][j][k][1] = cos(FlowAngle*(j-2)/(Nj-5))*
						cos(theta)*Density*VelocityMagnitude;
						ConservedVariablesNew[i][j][k][2] = -sin(FlowAngle*(j-2)/(Nj-5))*
						cos(theta)*Density*VelocityMagnitude ;
						ConservedVariablesNew[i][j][k][3] = -sin(theta)*Density*
						VelocityMagnitude;
						ConservedVariablesNew[i][j][k][4] = Pressure/(SpecificHeatRatio-1)
						+0.5*Density*VelocityMagnitude*VelocityMagnitude ;
					}
				}
			}

			// Subsonic Initial condition in the converging part
			for (int i = 0; i <=throat_location; ++i)
			{
				for (int j = 0; j < Nj; ++j)
				{
					for (int k = 0; k < Nk; ++k)
					{
						#if 1
						ConservedVariables[i][j][k][0] = DensityStagnation;
						ConservedVariables[i][j][k][1] = 0; 
						ConservedVariables[i][j][k][2] = 0; 
						ConservedVariables[i][j][k][3] = 0; 
						ConservedVariables[i][j][k][4] = PressureStagnation/(SpecificHeatRatio-1) ;

						ConservedVariablesNew[i][j][k][0] = DensityStagnation;
						ConservedVariablesNew[i][j][k][1] = 0; 
						ConservedVariablesNew[i][j][k][2] = 0; 
						ConservedVariablesNew[i][j][k][3] = 0; 
						ConservedVariablesNew[i][j][k][4] = PressureStagnation/(SpecificHeatRatio-1) ;
						#endif
						#if 0
						ConservedVariables[i][j][k][0] = 1.16;
						ConservedVariables[i][j][k][1] = 0 ;
						ConservedVariables[i][j][k][2] = 0 ;
						ConservedVariables[i][j][k][3] = 0 ;
						ConservedVariables[i][j][k][4] = 271767;

						ConservedVariablesNew[i][j][k][0] = 1.16;
						ConservedVariablesNew[i][j][k][1] = 0 ;
						ConservedVariablesNew[i][j][k][2] = 0 ;
						ConservedVariablesNew[i][j][k][3] = 0 ;
						ConservedVariablesNew[i][j][k][4] = 271767;
						#endif
					}
				}
			}
		break;
		#if 0
		case 3 :
		/**@bug Every time simulation starts from first iteration. So, to save the
		simulation it is good to start from the last solution as the initial 
		condition*/
			cout << "Do you wants to start the simulation, where you left?(enter y)" <<
			endl << "Other wise press any key"<< endl;
			char oldstart;
			cin>> oldstart;
			if (oldstart=='y')
			{
				ifstream nozzleData ("restart.csv");
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
								ConservedVariables[i][j][k][l] = atof(aline.c_str());
								ConservedVariablesNew[i][j][k][l] = 
								ConservedVariables[i][j][k][l];
							}
						}
					}
				}	
			}
		break;
		#endif
	}
	
}