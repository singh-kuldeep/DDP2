#include <fstream>
#include "math.h"
#include <iostream>
#include "math.h"
#include <fstream> /* For file handling */
#include <string> /* For strings */
#include <vector> /* For vectors*/
#include <cstdlib>
using namespace std;

// will return the Mach numbers in the diverging section 
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

// will return the Mach numbers in the diverging section 
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

int main()
{

	double MachConverging = 0.5;
	double MachDiverging = 5;

	std::vector<std::vector<double> > UpperWallCoordinates;
	std::vector<std::vector<double> > LowerWallCoordinates;

	ifstream xup("./../NozzleGeomatryGenrator/XCoordinatesUpperWall.csv");
	ifstream yup("./../NozzleGeomatryGenrator/YCoordinatesUpperWall.csv");

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
	double area ; 
	cout << "size " << UpperWallCoordinates.size() << endl;
	for (int i =0; i<throat_location; ++i)
	{	
		area = UpperWallCoordinates[i][1]/throat_area;
		// cout << getMachConvergingDuct(area) << endl;
	}

	for(int i=throat_location; i < UpperWallCoordinates.size(); ++i)
	{
		area = UpperWallCoordinates[i][1]/throat_area;
		// cout << getMachDivergingDuct(area) << endl;
	}

	// cout << "throat_area :  " << throat_area << endl;
	// cout << "Inlet mach :  " << getMachConvergingDuct(1.4786) << endl; // Me=4
	cout << "Inlet mach :  " << getMachConvergingDuct(2.2603) << endl; // Me=2
	// cout << "atan(1)  " << atan((1)/(1))*180/acos(-1) << endl;
	return 0;
}