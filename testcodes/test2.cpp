
#include <iostream>
#include "math.h"
#include <fstream> /* For file handling */
#include <string> /* For strings */
#include <vector> /* For vectors*/
#include <cstdlib> /* For converting string into numerical value */

using namespace std;
int main()
{	
	std::vector<std::vector<double> > UpperCoordinates;
	std::vector<std::vector<double> > DownCoordinates;

	ifstream xup("/home/kullu/Desktop/Acad/SEM10/DDP2/Code_DDP2/DDP2/NozzleGeomatryGenrator/XCoordinatesUpperWall.csv");
	ifstream yup("/home/kullu/Desktop/Acad/SEM10/DDP2/Code_DDP2/DDP2/NozzleGeomatryGenrator/YCoordinatesUpperWall.csv");

	int j = 0;
	   
   while(!xup.eof()){

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
   

 //   // writeing coordinates into the file 
	// ofstream Coordinates ;
	// Coordinates.open("CoordinatesNozzle.csv");
	// for (int i = 0; i < UpperCoordinates.size(); ++i)
	// {
	// 	Coordinates << UpperCoordinates[i][0] << "," << UpperCoordinates[i][1] << "," << DownCoordinates[i][1] << endl ; 
	// }

#if 0
   	std::vector<std::vector<double> > UpperCoordinatesNew;
	std::vector<std::vector<double> > DownCoordinatesNew;

	UpperCoordinatesNew.push_back(UpperCoordinates[0]); 
	// starting point is same
	DownCoordinatesNew.push_back(DownCoordinates[0]);
   	

   double dx;
   double dy;
   double x = UpperCoordinates[0][0];
   double y = UpperCoordinates[0][1];
   int N = 25 ; /**@param N Total cells in j direction*/
   /**@param N+1 Total grid points in j direction after including the 
   boundary points*/
   int i = 0 ;

   while(UpperCoordinates[UpperCoordinates.size()-1][0]>x)
   {
	   std::vector<double> xyup;
	   std::vector<double> xydown;
	   // cout<<UpperCoordinatesNew[i][0]<< "   " <<UpperCoordinatesNew[i][1]<<
	   //"   "<<DownCoordinatesNew[i][0]<<"   "<<DownCoordinatesNew[i][1]<<endl;

	   dy = (UpperCoordinatesNew[i][1] - DownCoordinatesNew[i][1])/N ;
	   dx = dy ;
	   
	   x = x + dx ;

	   if(UpperCoordinates[UpperCoordinates.size()-1][0]>x){
		   y = findY(x, UpperCoordinates);
		   xyup.push_back(x);
		   xyup.push_back(y);
		   UpperCoordinatesNew.push_back(xyup);
		
		   xydown.push_back(x);
		   xydown.push_back(0.0);
			   
 		   DownCoordinatesNew.push_back(xydown);
		   i++;
	   }
   }

#endif
   return 0;
}