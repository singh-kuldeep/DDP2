#include <fstream>
#include "math.h"
#include <iostream>
#include "math.h"
#include <fstream> /* For file handling */
#include <string> /* For strings */
#include <vector> /* For vectors*/
#include <cstdlib>
using namespace std;

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

	for (int i = 0; i < UpperWallCoordinates.size(); ++i)
	{
		cout << UpperWallCoordinates[i][1]/throat_area << endl;
	}
	cout << "throat_area :  " << throat_area << endl;
	cout << "throat_location :  " << throat_location << endl;
	return 0;

}
