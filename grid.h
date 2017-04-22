/*! \file  grid_ideal_nozzle.h
    \brief This header file functions find the grid points, cell area vectors 
    and the cell volumes for the ideal nozzle. In which nozzle has ans uniform
    flow at the exit. Because in the cancellation of expansion fan has been
    done by compression waves.  
    \author Kuldeep Singh
    \date 2017
    \warning For different geometries change this file accordingly.
*/

#include <iostream>
#include "math.h"
#include <fstream> /* For file handling */
#include <string> /* For strings */
#include <vector> /* For vectors*/
#include <cstdlib> /* For converting string into numerical value */

using namespace std;
// using std::cout;
// using std::endl;
// using std::ifstream;

/** \brief Find the cell side in z direction by taking average of all dx for deltaz.
*\param [in] DownCoordinates (x,y)coordinates of the down wall of the nozzle.
*\return double
*/
double finddeltaz(std::vector<std::vector<double> > DownCoordinates)
{	
	int size = DownCoordinates.size();
	double deltaz = 0;
	for (int i = 1; i < size; ++i)
	{
		deltaz = deltaz + DownCoordinates[i][0] - DownCoordinates[i-1][0];  
	}
	return deltaz/(size-1);
} 

/** \brief Calculates the distance between the two points in 3D space.
*\param [in] p1 First point.
*\param [in] p2 Second point.
*\return Distance between the two points
*/
double distance(std::vector<double> p1, std::vector<double> p2)
{
	return sqrt(pow((p2[0]-p1[0]),2) + pow((p2[1]-p1[1]),2) + 
		pow((p2[2]-p1[2]),2));
}

/** \brief This function will take boundary cell grid points and will calculate 
the ghost cell grid points by taking the mirror image about the boundary.
*\param [in] &x Pointer to x coordinate after taking mirror image
*\param [in] &y Pointer to y coordinate after taking mirror image
*\param [in] l x coordinate of the point which is to mirrored
*\param [in] m y coordinate of the point which is to mirrored
*\param [in] (x1,y1) Starting point of the line about which mirror is taken
*\param [in] (x2,y2) End point of the line about which mirror is taken
*\return void
*/
void takeMirror(
	double &x,double&y,
	double x1,double y1,
	double x2,double y2,
	double l,double m)
{
	double slop = (y2-y1)/(x2-x1);
	y = (slop*(2*l-2*x1)+m*slop*slop+2*y1-m)/(1+slop*slop);
	x =  (m-y)*slop+l;
} 


double find_y(double x, std::vector<std::vector<double> > UpperCoordinates){
	int i=0;
	while(x>UpperCoordinates[i][0]){
		i++;
	}
	return UpperCoordinates[i-1][1] + (UpperCoordinates[i][1] - 
		UpperCoordinates[i-1][1])*(x-UpperCoordinates[i-1][0])/
	(UpperCoordinates[i][0] - UpperCoordinates[i-1][0]);
}
/** \brief This function calculates the cell area and the cell volumes of all
*cells including the ghost cells.
*\param [in] iFaceAreaVectorIn Input pointer to "i" faces area vector 
*\param [in] jFaceAreaVectorIn Input pointer to "j" faces area vector 
*\param [in] kFaceAreaVectorIn Input pointer to "k" faces area vector 
*\param [in] CellVolumeIn Input pointer to cell volumes 
*\param [in] dsIn Input pointer to minimum distance    
*\param UpperCoordinates Upper wall coordinates (x,y) of the nozzle geometry
*\param DownCoordinates Down wall coordinates (x,y) of the nozzle geometry
/**\param [in] Ni Number of cells(Including ghost cells) in "i" direction.*/ 
// extra 4 is added for ghost cell
/**\param [in] Nj Number of cells(Including ghost cells) in "j" direction.*/ 
/**\param [in] Nk Number of cells(Including ghost cells) in "k" direction.*/
/* Here Nk = 6 because this is 2D-simulation so no need to take large 
number of cells in z direction
*\return void
*/

// Function defines the area vector and cell volumes 
void grid(vector<vector<vector<vector<double> > > > & iFaceAreaVectorIn,
		  vector<vector<vector<vector<double> > > > & jFaceAreaVectorIn,
		  vector<vector<vector<vector<double> > > > & kFaceAreaVectorIn,
		  vector<vector<vector<double> > > & CellVolumeIn,
		  vector<vector<vector<double> > > & dsIn,
		  int & Ni, int & Nj, int & Nk, string GeometryOption)
{
	cout << "Generating grid for " << GeometryOption << endl ;
	
	// Creating a 4D vector object for grid points
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> Dim4;

	/**\brief Only declaration is being done here initialization will be done in
	inside the individual geometry option*/
	Dim4 Coordinate, iFaceAreaVector, jFaceAreaVector, kFaceAreaVector;
	Dim3 CellVolume, ds;
	/**\param Coordinate 4D vector which stores the all coordinates of all cells
	(Including ghost) inside the domain*/
	/**\param iFaceAreaVector 4D vector which stores the all "i" face area 
	vectors of all cells(Including ghost) inside the domain*/
	/**\param jFaceAreaVector 4D vector which stores the all "j" face area 
	vectors of all cells(Including ghost) inside the domain*/
	/**\param kFaceAreaVector 4D vector which stores the all "k" face area 
	vectors of all cells(Including ghost) inside the domain*/
	/**\param CellVolume 3D vector which stores the cell volume of all cells
	(Including ghost) inside the domain*/

	/**@param N Total cells in all direction somewhat depends on N*/
	int N = 10; 
	Ni =3*N+4 ; 
	Nj =N+4 ;  
	Nk = 2+4 ; 
	/**\warning Do not reduce Nk below 6(atleast 2 live cells)*/
	
	// Default delta
	float delta = 1 ; 
	float deltax = delta ;
	float deltay = delta ;
	float deltaz = delta ;

	// Here only live cell coordinates will be defined
	if(GeometryOption =="StraightDuct")// straight duct
	{
		cout << "Generating grid for " << GeometryOption << endl ;
		N = 21 ;
		Ni = 3*N+4 ;
		Nj = N+4 ;  
		Nk = 2+4 ; 

		// Resizing the vectors
		Coordinate.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		

		for (int i =2; i <= Ni+1-2; ++i) 
		{
			for (int  j=2;  j <= Nj+1-2; j++)
			{
				for (int  k=2;  k <= Nk+1-2; k++)
				{
					Coordinate[i][j][k][0] = i*deltax ;   
					Coordinate[i][j][k][1] = j*deltay ;
					Coordinate[i][j][k][2] = k*deltaz ;
				}
			}	
		}
	}
	else if(GeometryOption == "BumpInsidetheStraightSuct") // Bump inside the straight duct
	{
		cout << "Generating grid for " << GeometryOption << endl ;
		
		/**\warning To increase the grid density, change the "N"*/
		/**\warning Do not change the Ni and Nj otherwise you will have to 
		change the code for grid as well, written inside the for loop below*/
		N = 21 ;
		Ni = 3*N+4 ;
		Nj = N+4 ;  
		Nk = 2+4 ; 

		// Resizing the vectors
		Coordinate.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		

		double wedgangle = 10*acos(-1)/180; // In radian
		double xtantheta = 0 ;
		for (int i = 2; i <= Ni+1-2; ++i) // 3N+4+1 grid points in x
		{
			for (int  j=2;  j <= Nj+1-2; j++)
			{
				for (int  k=2;  k <= Nk+1-2; k++)
				{
					Coordinate[i][j][k][0] = i*deltax ;   
					Coordinate[i][j][k][2] = k*deltaz ;
					if (i < N+2)
					{
						Coordinate[i][j][k][1] = (j-2)*deltay ;
					}
					else if( i >= N+2 && i <= floor(3*N/2) + 2 )
					{
						xtantheta = (i-N-1)*deltax*tan(wedgangle);

						Coordinate[i][j][k][1] = xtantheta + (j-2)*((Nj-4)*
							deltay - xtantheta)/(Nj-4) ; 
					}
					else 
					{
						Coordinate[i][j][k][1] = Coordinate[3*N+4-i][j][k][1] ;
					}
				}
			}	
		}
	}
	else if(GeometryOption =="IdelNozzleDesignedUsingMOC") // Nozzle 
	{
		cout << "Generating grid for " << GeometryOption << endl ;

		std::vector<std::vector<double> > UpperCoordinates;
		std::vector<std::vector<double> > DownCoordinates;
	   	ifstream nozzleData("./NozzleGeomatryGenrator/CoordinatesUpperWall.csv");
		int j = 0;
		
		// Reading the nozzle upper wall coordinates from the file
		while(!nozzleData.eof())
		{
		   string aline,xst,yst;
		   int comma_pos;
		   double xt,yt;
		   getline(nozzleData,aline);
		   comma_pos = aline.find(',',0);
		   xst = aline.substr(0,comma_pos);
		   yst = aline.substr(comma_pos+1,aline.length() - comma_pos - 1);
		   xt = atof( xst.c_str() );
		   yt = atof( yst.c_str() );
		   
		   vector<double> temp;
		   temp.push_back(xt); 
		   temp.push_back(yt);
		   UpperCoordinates.push_back(temp);

		   temp[1] = 0.0; // change the y only and push it to the Down vector
		   DownCoordinates.push_back(temp);
		   ++j;
		}
		// randomly extra zeros at the end so to remove them pop is used
		UpperCoordinates.pop_back();
		DownCoordinates.pop_back();

		nozzleData.close();
	   
		std::vector<std::vector<double> > UpperCoordinatesNew;
		std::vector<std::vector<double> > DownCoordinatesNew;

		// starting point is same.So,
		UpperCoordinatesNew.push_back(UpperCoordinates[0]); 
		DownCoordinatesNew.push_back(DownCoordinates[0]);
			
		double dx;
		double dy;
		double x = UpperCoordinates[0][0];
		double y = UpperCoordinates[0][1];

		// total N+1 points after including the boundary points because N cells
		int N = 20 ; 
		int i = 0 ;

		/*This while loop improves the grid quality, by making the deltax and 
		deltay same at every location*/
		while(UpperCoordinates[UpperCoordinates.size()-1][0]>x)
		{
		   std::vector<double> xyup;
		   std::vector<double> xydown;
		   /* cout << UpperCoordinatesNew[i][0] << "   " << 
		   UpperCoordinatesNew[i][1] << "   " << DownCoordinatesNew[i][0] 
		   << "   " << DownCoordinatesNew[i][1] << endl;*/

		   dy = (UpperCoordinatesNew[i][1] - DownCoordinatesNew[i][1])/N ;
		   dx = dy ;
		   
		   x = x + dx ;

		   if(UpperCoordinates[UpperCoordinates.size()-1][0]>x)
		   {
			   y = find_y(x, UpperCoordinates);
			   xyup.push_back(x);
			   xyup.push_back(y);
			   UpperCoordinatesNew.push_back(xyup);
			
			   xydown.push_back(x);
			   xydown.push_back(0.0);
				   
			   DownCoordinatesNew.push_back(xydown);
			   i++;
		    }
		}

	   	// N = 20 ;
	   	Ni = UpperCoordinatesNew.size()-1+4; // Total cells in X-dir  
		Nj = N+4 ;  // Total cell in Y-dir 
		Nk = 2+4 ; 	/* Because this is 2D-simulation so no need to take large 
		number of grids in z direction */

		// Resizing the vectors
		Coordinate.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));

		deltaz = finddeltaz(DownCoordinatesNew);
		
		// First defining the grid points/coordinates
		for (int i =2; i < Ni+1-2; i++) 
		// Will extend remaining x dir'n  after 
		{
			for (int  j=2;  j < Nj+1-2; j++)  
			//Will extend remaining y dir'n  after
			{
				for (int  k=2;  k < Nk+1-2; k++)					
				{
					Coordinate[i][j][k][0] = DownCoordinatesNew[i+1-3][0] ;   
					Coordinate[i][j][k][1] = (j-2)*(
						UpperCoordinatesNew[i+1-3][1]- 
						DownCoordinatesNew[i+1-3][1])/N ;
					Coordinate[i][j][k][2] = (k-2)*deltaz;
				}
			}	
		}

		#if 0 // Trying to check ghost cells are fine or not
		// these for loop are for ghost and live cells together 
		for (int i = 2; i  <= Ni-2; ++i)
		{
			for (int  j = 2;  j <= Nj-2; ++j)
			{
				for (int  k = 0;  k <= 1; ++k)
				{
					Coordinate[i][j][k][0] = Coordinate[i][j][2][0];
					Coordinate[i][j][k][1] = Coordinate[i][j][2][1];
					Coordinate[i][j][k][2] = (k-2)*deltaz; 
				}	
				for (int k = Nk-1; k <= Nk; ++k)
				{
					Coordinate[i][j][k][0] = Coordinate[i][j][Nk-1][0];
					Coordinate[i][j][k][1] = Coordinate[i][j][Nk-1][1];
					Coordinate[i][j][k][2] = (k-2)*deltaz; 
				}
			}
		}			
		
		for (int j = 2; j  <= Nj-2; ++j)
		{
			for (int  k = 2;  k <= Nk-2; ++k)
			{
				for (int  i = 1;  i <= 0; --i)
				{
					Coordinate[i][j][k][0] = Coordinate[i+1][j][k][0]-
					(Coordinate[3][j][k][0]-Coordinate[2][j][k][0]);
					Coordinate[i][j][k][1] = Coordinate[i+1][j][k][1];
					Coordinate[i][j][k][2] = Coordinate[i+1][j][k][2]; 
				}	
				for (int i = Ni-1; i <= Ni; ++i)
				{
					Coordinate[i][j][k][0] = Coordinate[i-1][j][k][0]+
					(Coordinate[Ni-2][j][k][0]-Coordinate[Ni-3][j][k][0]);
					Coordinate[i][j][k][1] = Coordinate[i-1][j][k][1];
					Coordinate[i][j][k][2] = Coordinate[i-1][j][k][2]; 
				}
			}
		}					
		#endif

		/* writing the updated coordinate in the file, which will be used in 
		the initial condition implementation*/

		/** @brief Structure of grid out put file ("grids_Nozzle_2D.csv") 
		* - First line of the grid file will contain grid points
		*(excluding ghost cells) in x and y direction 
		* - This will exclude the ghost, only live cells or actual geometry 
		*points
		*/
		ofstream NewWallX ;
		ofstream NewWallY ;

		NewWallX.open("./NozzleGeomatryGenrator/XCoordinatesUpperWall.csv");
		NewWallY.open("./NozzleGeomatryGenrator/YCoordinatesUpperWall.csv");
		
		//taking the lower left corner for the plotting  
		for (int i = 0; i < UpperCoordinatesNew.size(); ++i)
		{
		
			NewWallX << UpperCoordinatesNew[i][0] << endl;
			NewWallY << UpperCoordinatesNew[i][1] << endl; 
		}   
	}
	

	/* here comes the live cells area vectors
	Once the grid coordinates are defined then defining the area vectors and 
	cell volumes will be same for all geometries*/
	// Resizing the vectors
	iFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
	jFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
	kFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
	CellVolume.resize(Ni+1,Dim2(Nj+1,Dim1(Nk+1)));
	ds.resize(Ni+1,Dim2(Nj+1,Dim1(Nk+1)));

	//live cells area vector
	for (int i = 2; i  <= Ni+1-2; ++i)
	{
		for (int  j = 2;  j <= Nj+1-2; ++j)
		{
			for (int  k = 2;  k <= Nk+1-2; ++k)
	
			{
				iFaceAreaVector[i][j][k][0] = (Coordinate[i][j+1][k][1]-
				Coordinate[i][j][k][1])*deltaz ;
				iFaceAreaVector[i][j][k][1] = 0 ;
				iFaceAreaVector[i][j][k][2] = 0 ;

				jFaceAreaVector[i][j][k][0] = -deltaz*(Coordinate[i+1][j][k][1]-
				Coordinate[i][j][k][1]) ;
				jFaceAreaVector[i][j][k][1] =  deltaz*(Coordinate[i+1][j][k][0]-
				Coordinate[i][j][k][0]) ;
				jFaceAreaVector[i][j][k][2] = 0 ;

				kFaceAreaVector[i][j][k][0] = 0 ; 
				kFaceAreaVector[i][j][k][1] = 0 ;
				kFaceAreaVector[i][j][k][2] = 0.5*(Coordinate[i+1][j][k][0] - 
				Coordinate[i][j][k][0])*( (Coordinate[i][j+1][k][1] -
				Coordinate[i][j][k][1]) + (Coordinate[i+1][j+1][k][1] - 
				Coordinate[i+1][j][k][1]) ); 
			}
		}	
	}

	// live cell volumes 
	for (int i = 2; i  < Ni-2; ++i)
	{
		for (int  j= 2;  j < Nj-2; ++j)
		{
			for (int  k= 2;  k < Nk-2; ++k)
			{
				CellVolume[i][j][k] = 0.5 * (Coordinate[i+1][j][k][0] - 
				Coordinate[i][j][k][0]) * ( (Coordinate[i][j+1][k][1] - 
				Coordinate[i][j][k][1]) + (Coordinate[i+1][j+1][k][1] - 
				Coordinate[i+1][j][k][1]) ) * deltaz ; 
			}
		}	
	}

	 
	// All (i, j, k direction) ghost cells area vectors and volumes
	/*This method of defining the area vector and volume for the ghost cells is 
	only valid when 
	1) i/x direction wall are flat
	2) variation in the j/y wall can be there 
	3) k/z direction wall are flat */
	
	// here comes the i_ghost cells area vectors and volumes
	for(int i=0; i<2; ++i)
	{
		for (int j = 2; j < Nj+1-2; ++j)
		{
			for (int k = 2; k < Nk+1-2; ++k)
			{
				iFaceAreaVector[i][j][k][0] = iFaceAreaVector[4-i][j][k][0];
				iFaceAreaVector[i][j][k][1] = 0 ;
				iFaceAreaVector[i][j][k][2] = 0 ;

				iFaceAreaVector[Ni-1+i][j][k][0] = 
				iFaceAreaVector[Ni-3-i][j][k][0];
				iFaceAreaVector[Ni-1+i][j][k][1] = 0 ;
				iFaceAreaVector[Ni-1+i][j][k][2] = 0 ;

				CellVolume[i][j][k] = CellVolume[3-i][j][k];
				CellVolume[Ni-2+i][j][k] = CellVolume[Ni-3-i][j][k];
			}
		}
	}	

	// here comes the j_ghost cells area vectors
	double x0,y0; 
	/**\param (x0,y0) Live cell coordinates which needs 
	to be mirrored to get the ghost cell coordinates*/
	double x1,y1;
	 /**\param (x1,y1) Next live cell coordinates which needs to be 
	 mirrored to get the ghost cell coordinates*/
	double l0,m0,l1,m1; 
	/**\param (l0,m0),(l1,m1) Line about which reflection needs to be taken */
	double rx0,ry0; /**\param (rx0,ry0) Ghost cell grid point*/ 
	double rx1,ry1; /**\param (rx1,ry1) Ghost cell next grid point*/
	
	for (int i = 2; i < Ni+1-2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			for (int k = 2; k < Nk+1-2; ++k)
			{
				// Bottom ghost cells 
				// M1 point 
				l0 = Coordinate[i][2][k][0];
				m0 = Coordinate[i][2][k][1];

				// M2 point  
				l1 = Coordinate[i+1][2][k][0];
				m1 = Coordinate[i+1][2][k][1];
				//Mirror has been taken about the line M1 M2  


				x0 = Coordinate[i][4-j][k][0];
				y0 = Coordinate[i][4-j][k][1];

				x1 = Coordinate[i+1][4-j][k][0];
				y1 = Coordinate[i+1][4-j][k][1];

				takeMirror(rx0,ry0,l0,m0,l1,m1,x0,y0);
				takeMirror(rx1,ry1,l0,m0,l1,m1,x1,y1);

				jFaceAreaVector[i][j][k][0] = -deltaz*(ry1-ry0) ;
				jFaceAreaVector[i][j][k][1] =  deltaz*(rx1-rx0) ;
				jFaceAreaVector[i][j][k][2] = 0;
				
				// // Assigning the mirrored points to the coordinates but this is 
				// //for the testing purposes only
				// Coordinate[i][j][k][0] = rx0;
				// Coordinate[i][j][k][1] = ry0;
				// Coordinate[i][j][k][2] = Coordinate[i][4-j][k][2];

				// Coordinate[i+1][j][k][0] = rx1;
				// Coordinate[i+1][j][k][1] = ry1;
				// Coordinate[i+1][j][k][2] = Coordinate[i+1][4-j][k][2];
				
				// Top ghost cells
				l0 = Coordinate[i][Nj-2][k][0];
				m0 = Coordinate[i][Nj-2][k][1];

				l1 = Coordinate[i+1][Nj-2][k][0];
				m1 = Coordinate[i+1][Nj-2][k][1];

				x0 = Coordinate[i][Nj-3-j][k][0];
				y0 = Coordinate[i][Nj-3-j][k][1];

				x1 = Coordinate[i+1][Nj-3-j][k][0];
				y1 = Coordinate[i+1][Nj-3-j][k][1];

				takeMirror(rx0,ry0,l0,m0,l1,m1,x0,y0);
				takeMirror(rx1,ry1,l0,m0,l1,m1,x1,y1);

				jFaceAreaVector[i][Nj-2+1+j][k][0] = -deltaz*(ry1-ry0) ;
				jFaceAreaVector[i][Nj-2+1+j][k][1] =  deltaz*(rx1-rx0) ;
				jFaceAreaVector[i][Nj-2+1+j][k][2] = 0;

				CellVolume[i][j][k] = CellVolume[i][3-j][k];
				CellVolume[i][Nj-2+j][k] = CellVolume[i][Nj-3-j][k];

				// // Assigning the mirrored points to the coordinates but this is 
				// //for the testing purposes only
				// Coordinate[i][j][k][0] = rx0;
				// Coordinate[i][j][k][1] = ry0;
				// Coordinate[i][j][k][2] = Coordinate[i][Nj-3-j][k][2];

				// Coordinate[i+1][j][k][0] = rx1;
				// Coordinate[i+1][j][k][1] = ry1;
				// Coordinate[i+1][j][k][2] = Coordinate[i+1][Nj-3-j][k][2];
			}
		}
	}

	// here comes the k_ghost cells area vectors
	for(int i=2; i<Ni+1-2; ++i)
	{
		for (int j = 2; j < Nj+1-2; ++j)
		{
			for (int k = 0; k < 2; ++k)
			{
				kFaceAreaVector[i][j][k][0] = 0;
				kFaceAreaVector[i][j][k][1] = 0 ;
				kFaceAreaVector[i][j][k][2] = kFaceAreaVector[i][j][4-k][2] ;

				kFaceAreaVector[i][j][k][0] = 0;
				kFaceAreaVector[i][j][k][1] = 0 ;
				kFaceAreaVector[i][j][Nk-1+k][2] = 
				kFaceAreaVector[i][j][Nk-3-k][2] ;

				CellVolume[i][j][k] = CellVolume[i][j][3-k];
				CellVolume[i][j][Nk-2+k] = CellVolume[i][j][Nk-3-k];


			}
		}
	}	
	
	
	#if 0
	/**@bug Yet to calculate the ds value properly*/
	for (int i = 1; i  < Ni-2; ++i)
	{
		for (int  j= 1;  j < Nj-2; ++j)
		{
			for (int  k= 1;  k < Nk-2; ++k)
			{
				//ds[i][j][k] = min(distance(&Coordinate[i][j][k],
				//&Coordinate[i+1][j][k]),distance(&Coordinate[i][j+1][k],
				//&Coordinate[i+1][j+1][k]),
				//distance(&Coordinate[i][j][k+1],&Coordinate[i+1][j][k+1]),
				//distance(&Coordinate[i][j+1][k+1],&Coordinate[i+1][j+1][k+1]),
				//distance(&Coordinate[i][j][k],&Coordinate[i][j+1][k]),
				//distance(&Coordinate[i+1][j][k],&Coordinate[i+1][j+1][k]),
				//distance(&Coordinate[i][j][k+1],&Coordinate[i][j+1][k+1]),
				//distance(&Coordinate[i+1][j][k+1],&Coordinate[i+1][j+1][k+1]),
				//distance(&Coordinate[i][j][k],&Coordinate[i][j][k+1]),
				//distance(&Coordinate[i+1][j][k],&Coordinate[i+1][j][k+1]),
				//distance(&Coordinate[i][j+1][k],&Coordinate[i][j+1][k+1]),
				//distance(&Coordinate[i+1][j+1][k],&Coordinate[i+1][j+1][k+1]));

				ds[i][j][k] = Coordinate[i+1][Nj/2][Nk/2][1] - 
				Coordinate[i+1][Nj/2 - 1][Nk/2][1];	
			}
		}	
	}
	
	// writeing ds into the file 
	ofstream kullu_ds ;
	kullu_ds.open("Nozzle_ds.csv");
	for (int i = 2; i < Ni-2; ++i)
	{
		kullu_ds << ds[i][Nj/2][Nk/2] << endl ; 
	}
	#endif 
	

	/** @brief Structure of grid out put file ("grids_Nozzle_2D.csv") 
	* - First line of the grid file will contain grid points
	*(excluding ghost cells) in x and y direction 
	* - This will exclude the ghost, only live cells or actual geomatry points
	*/
	ofstream kullu_grid ;
	kullu_grid.open("./Results/outputfiles/grids_2D.csv");
	kullu_grid << Ni-4 << "," << Nj-4 << endl ; 
	//taking the lower left corner for the plotting  
	for (int i = 2; i <= Ni-3; ++i)
	{
		for (int j = 2; j <= Nj-3; ++j)
		{
			// kullu_grid << 0.5*(Coordinate[i][j][4][0]+
			// Coordinate[i+1][j][4][0]) << "," 
			// << 0.5*(Coordinate[i][j][4][1]+Coordinate[i][j+1][4][1]) << endl;
			kullu_grid <<  Coordinate[i][j][Nk/2][0] << ","<< 
			Coordinate[i][j][Nk/2][1] << endl;
		}
	}   

	#if 0 /* plotting the full domain including the ghost cell just to 
	check the correctness of the ghost cell generation*/
	#if 1
	/*
	Because the corner points(Corner ghost cells) do not play any 
	role(in the simulation) so they 
	are not generated above, but to plot the complete domain, these corner 
	points are added and defined below  
	*/
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			for (int k = 0; k < Nk+1; ++k)
			{
				Coordinate[i][j][k][0] = Coordinate[i][2][k][0];
				Coordinate[i][j][k][1] = Coordinate[2][j][k][1];
				Coordinate[i][j][k][2] = Coordinate[2][2][k][2]; 

				Coordinate[Ni-i][j][k][0] = Coordinate[Ni-i][2][k][0];
				Coordinate[Ni-i][j][k][1] = Coordinate[Ni-2][j][k][1];
				Coordinate[Ni-i][j][k][2] = Coordinate[Ni-2][2][k][2];

				Coordinate[i][Nj-j][k][0] = Coordinate[i][Nj-2][k][0];
				Coordinate[i][Nj-j][k][1] = Coordinate[2][Nj-j][k][1]; 
				Coordinate[i][Nj-j][k][2] = Coordinate[2][Nj-2][k][2];

				Coordinate[Ni-i][Nj-j][k][0] = Coordinate[Ni-i][Nj-2][k][0];
				Coordinate[Ni-i][Nj-j][k][1] = Coordinate[Ni-2][Nj-j][k][1]; 
				Coordinate[Ni-i][Nj-j][k][2] = Coordinate[Ni-2][Nj-2][k][2];
			}
		}
	}

	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < Nj+1; ++j)
		{
			for (int k = 0; k < 2; ++k)
			{
				Coordinate[i][j][k][0] = Coordinate[i][j][2][0]; 
				Coordinate[i][j][k][1] = Coordinate[2][j][2][1];  
				Coordinate[i][j][k][2] = Coordinate[2][j][k][2];

				Coordinate[Ni-i][j][k][0] = Coordinate[Ni-i][j][2][0];
				Coordinate[Ni-i][j][k][1] = Coordinate[Ni-2][j][2][1]; 
				Coordinate[Ni-i][j][k][2] = Coordinate[Ni-2][j][k][2];

				Coordinate[i][j][Nk-k][0] = Coordinate[i][j][Nk-2][0]; 
				Coordinate[i][j][Nk-k][1] = Coordinate[2][j][Nk-2][1];  
				Coordinate[i][j][Nk-k][2] = Coordinate[2][j][Nk-k][2];

				Coordinate[Ni-i][j][Nk-k][0] = Coordinate[Ni-i][j][Nk-2][0];
				Coordinate[Ni-i][j][Nk-k][1] = Coordinate[Nk-2][j][Nk-2][1]; 
				Coordinate[Ni-i][j][Nk-k][2] = Coordinate[Ni-2][j][Nk-2][2];
			}
		}
	}
	for (int i = 0; i < Ni+1; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			for (int k = 0; k < 2; ++k)
			{
				Coordinate[i][j][k][0] = Coordinate[i][2][2][0]; 
				Coordinate[i][j][k][1] = Coordinate[i][j][2][1];  
				Coordinate[i][j][k][2] = Coordinate[i][2][k][2];

				Coordinate[i][Nj-j][k][0] = Coordinate[i][Nj-2][2][0];
				Coordinate[i][Nj-j][k][1] = Coordinate[i][Nj-j][2][1]; 
				Coordinate[i][Nj-j][k][2] = Coordinate[i][Nj-2][k][2];

				Coordinate[i][j][Nk-k][0] = Coordinate[i][2][Nk-2][0]; 
				Coordinate[i][j][Nk-k][1] = Coordinate[i][j][Nk-2][1];  
				Coordinate[i][j][Nk-k][2] = Coordinate[i][2][Nk-k][2];

				Coordinate[i][Nj-j][Nk-k][0] = Coordinate[i][Nj-2][Nk-2][0];
				Coordinate[i][Nj-j][Nk-k][1] = Coordinate[i][Nj-2][Nk-k][1]; 
				Coordinate[i][Nj-j][Nk-k][2] = Coordinate[i][Nj-2][Nk-k][2];
			}
		}
	}
	#endif
	ofstream kullu_grid_with_ghost ;
	kullu_grid_with_ghost.open("./Results/outputfiles/grids_2D_with_ghost.csv");
	kullu_grid_with_ghost << Ni+1 << "," << Nj+1 << endl ; 
	//taking the lower left corner for the plotting  
	for (int i = 0; i <= Ni; ++i)
	{
		for (int j = 0; j <= Nj; ++j)
		{
			// kullu_grid_with_ghost << 0.5*(Coordinate[i][j][4][0]+
			// Coordinate[i+1][j][4][0]) << "," 
			// << 0.5*(Coordinate[i][j][4][1]+Coordinate[i][j+1][4][1]) << endl;
			kullu_grid_with_ghost <<  Coordinate[i][j][Nk/2][0] << ","<< 
			Coordinate[i][j][Nk/2][1] << endl;
		}
	}
	#endif

	// assigning the vector pointers 
	iFaceAreaVectorIn = iFaceAreaVector;
	jFaceAreaVectorIn = jFaceAreaVector;
	kFaceAreaVectorIn = kFaceAreaVector;
	CellVolumeIn = CellVolume;
	dsIn = ds;
	cout << "Grid genration succsesful :)  " << endl ;
}