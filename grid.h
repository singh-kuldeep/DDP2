/*! \file  grid_ideal_nozzle.h
    \brief This header file functions find the live cell vertices, live area 
    vectors and live cell volumes for the multiple geomatries based on the input 
    given.
    \author Kuldeep Singh
    \date 2017
    \warning For different geometries chnage the option in the input file.
*/

#include <iostream>
#include "math.h"
#include <fstream> /* For file handling */
#include <string> /* For strings */
#include <vector> /* For vectors*/
#include <cstdlib> /* For converting string into numerical value */

using namespace std;

/** \brief This function calculates the cell area and the cell volumes of all
*cells including the ghost cells.
*\param [in] iFaceAreaVectorIn Input pointer to "i" faces area vector 
*\param [in] jFaceAreaVectorIn Input pointer to "j" faces area vector 
*\param [in] kFaceAreaVectorIn Input pointer to "k" faces area vector 
*\param [in] CellVolumeIn Input pointer to cell volumes 
*\param [in] dsIn Input pointer to minimum distance    
*\param UpperCoordinates Upper wall coordinates (x,y) of the nozzle geometry
*\param DownCoordinates Down wall coordinates (x,y) of the nozzle geometry
/**\param [in] Ni Number of live cells in "i" direction.*/ 
/**\param [in] Nj Number of live cells in "j" direction.*/ 
/**\param [in] Nk Number of live cells in "k" direction.*/
/* Here Nk = 2 because this is 2D-simulation so no need to take large 
number of cells in z direction
*\return void
*/

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

double find_y(double x, std::vector<std::vector<double> > UpperCoordinates){
	int i=0;
	while(x>UpperCoordinates[i][0]){
		i++;
	}
	return UpperCoordinates[i-1][1] + (UpperCoordinates[i][1] - 
		UpperCoordinates[i-1][1])*(x-UpperCoordinates[i-1][0])/
	(UpperCoordinates[i][0] - UpperCoordinates[i-1][0]);
}

// Function defines the area vector and cell volumes 
void grid(vector<vector<vector<vector<double> > > > & Coordinate,
		  vector<vector<vector<vector<double> > > > & iFaceAreaVector,
		  vector<vector<vector<vector<double> > > > & jFaceAreaVector,
		  vector<vector<vector<vector<double> > > > & kFaceAreaVector,
		  vector<vector<vector<double> > > & CellVolume,
		  vector<vector<vector<double> > > & ds,
		  int & Ni, int & Nj, int & Nk, string GeometryOption)
{	
	// Creating a 4D vector object for grid points
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> Dim4;

	/**\brief Only declaration is being done here initialization will be done in
	inside the individual geometry option*/
	// Dim4 iFaceAreaVector, jFaceAreaVector, kFaceAreaVector;
	// Dim3 CellVolume, ds;

	/**\param Coordinate 4D vector, stores the all coordinates of live cells*/
	/**\param iFaceAreaVector 4D vector, stores the all "i" face area  
	vectors of live cells*/
	/**\param jFaceAreaVector 4D vector which stores the all "j" face area 
	vectors of live cells*/
	/**\param kFaceAreaVector 4D vector which stores the all "k" face area 
	vectors of live cells*/
	/**\param CellVolume 3D vector, stores the cell volume live cells */

	/**@param N Total cells in all direction somewhat depends on N*/
	int N = 100; 
	Ni = N ; 
	Nj = 10 ;  
	Nk = 1 ; 
	/**\warning Do not reduce Nk below 6(atleast 2 live cells)*/
	
	// Default delta
	float delta = 1 ; 
	float deltax = delta ;
	float deltay = delta ;
	float deltaz = delta ;

	// Here only live cell coordinates will be defined
	if(GeometryOption =="StraightDuct")// straight duct
	{
		// Resizing the vectors
		Coordinate.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		
		for (int i = 0; i < Ni+1; ++i) // 3N+4+1 grid points in x
		{
			for (int  j=0;  j < Nj+1; j++)
			{
				for (int  k=0;  k < Nk+1; k++)
				{
					Coordinate[i][j][k][0] = i*deltax ;   
					Coordinate[i][j][k][1] = j*deltay ;
					Coordinate[i][j][k][2] = k*deltaz ;
				}
			}	
		}

		/* here comes the live cells area vectors*/
		// Resizing the vectors
		iFaceAreaVector.resize(Ni+1,Dim3(Nj,Dim2(Nk,Dim1(3))));
		jFaceAreaVector.resize(Ni,Dim3(Nj+1,Dim2(Nk,Dim1(3))));
		kFaceAreaVector.resize(Ni,Dim3(Nj,Dim2(Nk+1,Dim1(3))));
		CellVolume.resize(Ni,Dim2(Nj,Dim1(Nk)));
		ds.resize(Ni,Dim2(Nj,Dim1(Nk)));
		//i 
		for (int i = 0; i<Ni+1; ++i)
		{
			for (int  j = 0;  j <Nj; ++j)
			{
				for (int  k = 0; k < Nk; ++k)
		
				{
					iFaceAreaVector[i][j][k][0] = deltay*deltaz ;
					iFaceAreaVector[i][j][k][1] = 0 ;
					iFaceAreaVector[i][j][k][2] = 0 ;
				}
			}	
		}
		cout << "iFaceAreaVector OK\n";
		//j
		for (int i = 0; i<Ni; ++i)
		{
			for (int  j = 0;  j <Nj+1; ++j)
			{
				for (int  k = 0; k < Nk; ++k)
		
				{
					jFaceAreaVector[i][j][k][0] = 0;
					jFaceAreaVector[i][j][k][1] = deltaz*deltax;
					jFaceAreaVector[i][j][k][2] = 0 ;					
				}
			}	
		}
		cout << "jFaceAreaVector OK\n";

		//k
		for (int i = 0; i<Ni; ++i)
		{
			for (int  j = 0;  j <Nj; ++j)
			{
				for (int  k = 0; k < Nk+1; ++k)
		
				{
					kFaceAreaVector[i][j][k][0] = 0 ; 
					kFaceAreaVector[i][j][k][1] = 0 ;
					kFaceAreaVector[i][j][k][2] = deltax*deltay; 					
				}
			}	
		}
		cout << "kFaceAreaVector OK\n";

		cout << "Generating grid for " << GeometryOption << endl ;

		// live cell volumes 
		for (int i = 0; i < Ni; ++i)
		{
			for (int  j= 0; j < Nj; ++j)
			{
				for (int  k= 0; k < Nk; ++k)
				{
					CellVolume[i][j][k] = deltax*deltay*deltaz; 
					ds[i][j][k] = delta;
				}
			}	
		}

		// writeing ds into the file 
		ofstream kullu_ds ;
		kullu_ds.open("./Results/outputfiles/ds.csv");
		for (int i = 0; i < Ni; ++i)
		{
			kullu_ds << ds[i][Nj/2][Nk/2] << endl ; 
		}
		
		cout << "grid has been genrated for " << GeometryOption << endl ;
		
		/** @brief Structure of grid out put file ("grids_Cell_Center_xy_plane.csv") 
		* - First line : Number of cell center in i direction and j direction 
		* - Coordinates of the cell centers 
		*/
		ofstream CellCenter ;
		CellCenter.open("./Results/outputfiles/CellCenter_ij.csv");
		CellCenter << Ni << "," << Nj << endl ; 
		//taking the lower left corner for the plotting  
		for (int i = 0; i < Ni; ++i)
		{
			for (int j = 0; j < Nj; ++j)
			{
				CellCenter << (Coordinate[i][j][0][0] + Coordinate[i+1][j][0][0] + 
							Coordinate[i][j+1][0][0] + Coordinate[i+1][j+1][0][0])/4 
				<< ","<< (Coordinate[i][j][0][1] + Coordinate[i+1][j][0][1] + 
							Coordinate[i][j+1][0][1] + Coordinate[i+1][j+1][0][1])/4
				<< endl;
			}
		}   
	}
	cout << "Grid genration succsesful for" << GeometryOption << ":)" << endl ;

	// Bump inside the straight duct
	// Here only live cell coordinates will be defined
	if(GeometryOption =="BumpInsidetheStraightSuct")// straight duct
	{
		/**\warning To increase the grid density, change the "N"*/
		/**\warning Do not change the Ni and Nj otherwise you will have to 
		change the code for grid as well, written inside the for loop below*/
		N = 21 ;
		Ni = 3*N ;
		Nj = N ;  
		Nk = 1 ;

		// Resizing the vectors
		Coordinate.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		
		for (int i = 0; i < Ni+1; ++i) // 3N+4+1 grid points in x
		{
			for (int  j=0;  j < Nj+1; j++)
			{
				for (int  k=0;  k < Nk+1; k++)
				{
					Coordinate[i][j][k][0] = i*deltax ;   
					Coordinate[i][j][k][2] = k*deltaz ;
					if (i < N)
					{
						Coordinate[i][j][k][1] = j*deltay ;
					}
					else if( i >= N && i <= floor(3*N/2))
					{
						double wedgangle = 10*acos(-1)/180; // In radian
						double xtantheta = (i-N)*deltax*tan(wedgangle);

						Coordinate[i][j][k][1] = xtantheta + j*(Nj*deltay - xtantheta)/Nj ; 
					}
					else 
					{
						Coordinate[i][j][k][1] = Coordinate[Ni-i][j][k][1] ;
					}
				}
			}	
		}

		/* here comes the live cells area vectors*/
		// Resizing the vectors
		iFaceAreaVector.resize(Ni+1,Dim3(Nj,Dim2(Nk,Dim1(3))));
		jFaceAreaVector.resize(Ni,Dim3(Nj+1,Dim2(Nk,Dim1(3))));
		kFaceAreaVector.resize(Ni,Dim3(Nj,Dim2(Nk+1,Dim1(3))));
		CellVolume.resize(Ni,Dim2(Nj,Dim1(Nk)));
		ds.resize(Ni,Dim2(Nj,Dim1(Nk)));
		//i 
		for (int i = 0; i<Ni+1; ++i)
		{
			for (int  j = 0;  j <Nj; ++j)
			{
				for (int  k = 0; k < Nk; ++k)
		
				{
					iFaceAreaVector[i][j][k][0] = (Coordinate[i][j+1][k][1]-
					Coordinate[i][j][k][1])*deltaz;
					iFaceAreaVector[i][j][k][1] = 0 ;
					iFaceAreaVector[i][j][k][2] = 0 ;
				}
			}	
		}
		cout << "iFaceAreaVector OK\n";
		//j
		for (int i = 0; i<Ni; ++i)
		{
			for (int  j = 0;  j <Nj+1; ++j)
			{
				for (int  k = 0; k < Nk; ++k)
		
				{
					jFaceAreaVector[i][j][k][0] = -deltaz*(Coordinate[i+1][j][k][1]-
					Coordinate[i][j][k][1]) ;
					jFaceAreaVector[i][j][k][1] =  deltaz*(Coordinate[i+1][j][k][0]-
					Coordinate[i][j][k][0]) ;
					jFaceAreaVector[i][j][k][2] = 0 ;					
				}
			}	
		}
		cout << "jFaceAreaVector OK\n";

		//k
		for (int i = 0; i<Ni; ++i)
		{
			for (int  j = 0;  j <Nj; ++j)
			{
				for (int  k = 0; k < Nk+1; ++k)
		
				{
					kFaceAreaVector[i][j][k][0] = 0 ; 
					kFaceAreaVector[i][j][k][1] = 0 ;
					kFaceAreaVector[i][j][k][2] = 0.5*(Coordinate[i+1][j][k][0] - 
					Coordinate[i][j][k][0])*( (Coordinate[i][j+1][k][1] -
					Coordinate[i][j][k][1]) + (Coordinate[i+1][j+1][k][1] - 
					Coordinate[i+1][j][k][1]) );  					
				}
			}	
		}
		cout << "kFaceAreaVector OK\n";

		cout << "Generating grid for " << GeometryOption << endl ;

		// live cell volumes 
		for (int i = 0; i < Ni; ++i)
		{
			for (int  j= 0; j < Nj; ++j)
			{
				for (int  k= 0; k < Nk; ++k)
				{
					CellVolume[i][j][k] = 0.5 * (Coordinate[i+1][j][k][0] - 
					Coordinate[i][j][k][0]) * ( (Coordinate[i][j+1][k][1] - 
					Coordinate[i][j][k][1]) + (Coordinate[i+1][j+1][k][1] - 
					Coordinate[i+1][j][k][1]) ) * deltaz ; 
					double dyleft = (Coordinate[i][j+1][k][1]-Coordinate[i][j][k][1]);
					double dyright = (Coordinate[i+1][j+1][k][1]-Coordinate[i+1][j][k][1]);
					ds[i][j][k] = fmin(dyleft,dyright);
				}
			}	
		}

		// writeing ds into the file 
		ofstream kullu_ds ;
		kullu_ds.open("./Results/outputfiles/ds.csv");
		for (int i = 0; i < Ni; ++i)
		{
			kullu_ds << ds[i][Nj/2][Nk/2] << endl ; 
		}
		
		cout << "grid has been genrated for " << GeometryOption << endl ;
		
		/** @brief Structure of grid out put file ("grids_Cell_Center_xy_plane.csv") 
		* - First line : Number of cell center in i direction and j direction 
		* - Coordinates of the cell centers 
		*/
		ofstream CellCenter ;
		CellCenter.open("./Results/outputfiles/CellCenter_ij.csv");
		CellCenter << Ni << "," << Nj << endl ; 
		//taking the lower left corner for the plotting  
		for (int i = 0; i < Ni; ++i)
		{
			for (int j = 0; j < Nj; ++j)
			{
				CellCenter << (Coordinate[i][j][0][0] + Coordinate[i+1][j][0][0] + 
							Coordinate[i][j+1][0][0] + Coordinate[i+1][j+1][0][0])/4 
				<< ","<< (Coordinate[i][j][0][1] + Coordinate[i+1][j][0][1] + 
							Coordinate[i][j+1][0][1] + Coordinate[i+1][j+1][0][1])/4
				<< endl;
			}
		}   
	}
	cout << "Grid genration succsesful for" << GeometryOption << ":)" << endl ;


	// Nozzle
	if(GeometryOption =="IdelNozzleDesignedUsingMOC")// straight duct
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
		int N = 10 ; 
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
	   	Ni = UpperCoordinatesNew.size()-1; // Total cells in X-dir  
		Nj = N ;  // Total cell in Y-dir 
		Nk = 2 ; 	/* Because this is 2D-simulation so no need to take large 
		number of grids in z direction */

		// Resizing the vectors
		Coordinate.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));

		deltaz = finddeltaz(DownCoordinatesNew);
		
		// First defining the grid points/coordinates
		for (int i =0; i < Ni+1; i++) 
		// Will extend remaining x dir'n  after 
		{
			for (int  j=0;  j < Nj+1; j++)  
			//Will extend remaining y dir'n  after
			{
				for (int  k=0;  k < Nk+1; k++)					
				{
					Coordinate[i][j][k][0] = DownCoordinatesNew[i][0] ;   
					Coordinate[i][j][k][1] = j*
					(UpperCoordinatesNew[i][1]-DownCoordinatesNew[i][1])/N ;
					Coordinate[i][j][k][2] = k*deltaz;
				}
			}	
		}

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

		/* here comes the live cells area vectors*/
		// Resizing the vectors
		iFaceAreaVector.resize(Ni+1,Dim3(Nj,Dim2(Nk,Dim1(3))));
		jFaceAreaVector.resize(Ni,Dim3(Nj+1,Dim2(Nk,Dim1(3))));
		kFaceAreaVector.resize(Ni,Dim3(Nj,Dim2(Nk+1,Dim1(3))));
		CellVolume.resize(Ni,Dim2(Nj,Dim1(Nk)));
		ds.resize(Ni,Dim2(Nj,Dim1(Nk)));
		//i 
		for (int i = 0; i<Ni+1; ++i)
		{
			for (int  j = 0;  j <Nj; ++j)
			{
				for (int  k = 0; k < Nk; ++k)
		
				{
					iFaceAreaVector[i][j][k][0] = (Coordinate[i][j+1][k][1]-
					Coordinate[i][j][k][1])*deltaz;
					iFaceAreaVector[i][j][k][1] = 0 ;
					iFaceAreaVector[i][j][k][2] = 0 ;
				}
			}	
		}
		cout << "iFaceAreaVector OK\n";
		//j
		for (int i = 0; i<Ni; ++i)
		{
			for (int  j = 0;  j <Nj+1; ++j)
			{
				for (int  k = 0; k < Nk; ++k)
		
				{
					jFaceAreaVector[i][j][k][0] = -deltaz*(Coordinate[i+1][j][k][1]-
					Coordinate[i][j][k][1]) ;
					jFaceAreaVector[i][j][k][1] =  deltaz*(Coordinate[i+1][j][k][0]-
					Coordinate[i][j][k][0]) ;
					jFaceAreaVector[i][j][k][2] = 0 ;					
				}
			}	
		}
		cout << "jFaceAreaVector OK\n";

		//k
		for (int i = 0; i<Ni; ++i)
		{
			for (int  j = 0;  j <Nj; ++j)
			{
				for (int  k = 0; k < Nk+1; ++k)
		
				{
					kFaceAreaVector[i][j][k][0] = 0 ; 
					kFaceAreaVector[i][j][k][1] = 0 ;
					kFaceAreaVector[i][j][k][2] = 0.5*(Coordinate[i+1][j][k][0] - 
					Coordinate[i][j][k][0])*( (Coordinate[i][j+1][k][1] -
					Coordinate[i][j][k][1]) + (Coordinate[i+1][j+1][k][1] - 
					Coordinate[i+1][j][k][1]) );  					
				}
			}	
		}
		cout << "kFaceAreaVector OK\n";

		cout << "Generating grid for " << GeometryOption << endl ;

		// live cell volumes 
		for (int i = 0; i < Ni; ++i)
		{
			for (int  j= 0; j < Nj; ++j)
			{
				for (int  k= 0; k < Nk; ++k)
				{
					CellVolume[i][j][k] = 0.5 * (Coordinate[i+1][j][k][0] - 
					Coordinate[i][j][k][0]) * ( (Coordinate[i][j+1][k][1] - 
					Coordinate[i][j][k][1]) + (Coordinate[i+1][j+1][k][1] - 
					Coordinate[i+1][j][k][1]) ) * deltaz ; 
					double dyleft = (Coordinate[i][j+1][k][1]-Coordinate[i][j][k][1]);
					double dyright = (Coordinate[i+1][j+1][k][1]-Coordinate[i+1][j][k][1]);
					ds[i][j][k] = fmin(dyleft,dyright);
				}
			}	
		}

		// writeing ds into the file 
		ofstream kullu_ds ;
		kullu_ds.open("./Results/outputfiles/ds.csv");
		for (int i = 0; i < Ni; ++i)
		{
			kullu_ds << ds[i][Nj/2][Nk/2] << endl ; 
		}
		
		cout << "grid has been genrated for " << GeometryOption << endl ;
		
		/** @brief Structure of grid out put file ("grids_Cell_Center_xy_plane.csv") 
		* - First line : Number of cell center in i direction and j direction 
		* - Coordinates of the cell centers 
		*/
		ofstream CellCenter ;
		CellCenter.open("./Results/outputfiles/CellCenter_ij.csv");
		CellCenter << Ni << "," << Nj << endl ; 
		//taking the lower left corner for the plotting  
		for (int i = 0; i < Ni; ++i)
		{
			for (int j = 0; j < Nj; ++j)
			{
				CellCenter << (Coordinate[i][j][0][0] + Coordinate[i+1][j][0][0] + 
							Coordinate[i][j+1][0][0] + Coordinate[i+1][j+1][0][0])/4 
				<< ","<< (Coordinate[i][j][0][1] + Coordinate[i+1][j][0][1] + 
							Coordinate[i][j+1][0][1] + Coordinate[i+1][j+1][0][1])/4
				<< endl;
			}
		}   
	}
	cout << "Grid genration succsesful for" << GeometryOption << ":)" << endl ;
}