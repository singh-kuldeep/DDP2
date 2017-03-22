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
/* Here Nk = 5 because this is 2D-simulation so no need to take large 
number of cells in z direction
*\return void
*/

// Function defines the area vector and cell volumes 
void grid(vector<vector<vector<vector<double> > > > & iFaceAreaVectorIn,
		  vector<vector<vector<vector<double> > > > & jFaceAreaVectorIn,
		  vector<vector<vector<vector<double> > > > & kFaceAreaVectorIn,
		  vector<vector<vector<double> > > & CellVolumeIn,
		  vector<vector<vector<double> > > & dsIn,
		  int & Ni, int & Nj, int & Nk, int GeometryOption)
{
	/**@param N+1 Total "grid points" in j direction after including the 
	boundary points*/
	int N; /**@param N Total cells in j direction*/
	// Ni =3*N+4 ; //3*N+4 ;
	// Nj =N+4 ;//N+4 ;  
	// Nk = 2+4 ;
	
	// Creating a 4D vector object for grid points
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> Dim4;

	/**\brief Only declaration is being done here initialization will be done in
	inside the individual geometry option*/
	Dim4 Coordinate, iFaceAreaVector, jFaceAreaVector, kFaceAreaVector;
	Dim3 CellVolume, ds;
	
	#if 0
	/**\param Coordinate 4D vector which stores the all coordinates of all cells
	(Including ghost) inside the domain*/
	Dim4 Coordinate(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3)))); 
	
	/**\param iFaceAreaVector 4D vector which stores the all "i" face area 
	vectors of all cells(Including ghost) inside the domain*/
	Dim4 iFaceAreaVector(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3)))); 
	/**\param jFaceAreaVector 4D vector which stores the all "j" face area 
	vectors of all cells(Including ghost) inside the domain*/
	Dim4 jFaceAreaVector(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3)))); 
	/**\param kFaceAreaVector 4D vector which stores the all "k" face area 
	vectors of all cells(Including ghost) inside the domain*/
	Dim4 kFaceAreaVector(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3)))); 

	/**\param CellVolume 3D vector which stores the cell volume of all cells
	(Including ghost) inside the domain*/
	Dim3 CellVolume(Ni,Dim2(Nj,Dim1(Nk)));

	Dim3 ds(Ni,Dim2(Nj,Dim1(Nk)));
	#endif

	float delta = 0.1 ; 
	float deltax = delta ;
	float deltay = delta ;
	float deltaz = delta ;

	// Here only live cell coordinates will be defined
	if(GeometryOption ==1)// straight duct
	{
		N = 20 ;
		Ni = 3*N+4 ;
		Nj = N+4 ;  
		Nk = 2+4 ; 

		// Resizing the vectors
		Coordinate.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		iFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		jFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		kFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		CellVolume.resize(Ni+1,Dim2(Nj+1,Dim1(Nk+1)));
		ds.resize(Ni+1,Dim2(Nj+1,Dim1(Nk+1)));

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
	else if(GeometryOption == 2) // Bump inside the straight duct
	{
		/**\warning Do not change the Ni and Nj otherwise you will have to 
		change the code for grid as well written inside the for loop below*/
		N = 21 ;
		Ni = 3*N+4 ;
		Nj = N+4 ;  
		Nk = 2+4 ; 

		// Resizing the vectors
		Coordinate.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		iFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		jFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		kFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		CellVolume.resize(Ni+1,Dim2(Nj+1,Dim1(Nk+1)));
		ds.resize(Ni+1,Dim2(Nj+1,Dim1(Nk+1)));

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

						Coordinate[i][j][k][1] = xtantheta + (j-2)*((Nj-4)*deltay - xtantheta)/(Nj-4) ; 
						// Coordinate[i][j][k][1] = (j-2) * (deltay - 0.2*(i-N-2)/(Ni-2))  ; 
					}
					else 
					{
						Coordinate[i][j][k][1] = Coordinate[3*N+4-i][j][k][1] ;
					}
				}
			}	
		}
	}
	else if(GeometryOption == 3 || GeometryOption == 4) // Nozzle 
	{
		std::vector<std::vector<double> > UpperCoordinates;
		std::vector<std::vector<double> > DownCoordinates;

		ifstream xup("./NozzleGeomatryGenrator/XCoordinatesUpperWall.csv");
		ifstream yup("./NozzleGeomatryGenrator/YCoordinatesUpperWall.csv");

		int pointCount = 0;
		   
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
		   ++pointCount;
		}
		// rendomaly extra zeros at the end so to remove them pop is used
		UpperCoordinates.pop_back();
		DownCoordinates.pop_back();

		// closing the file
		xup.close();
		yup.close();


		#if 1
		// Grids for hypersonic nozzle which has uniform flow at the exit

		// int N; /**@param N Total cells in j direction*/
		/**@param N+1 Total "grid points" in j direction after including the 
		boundary points*/
		// N = 25;
		N = UpperCoordinates.size();
		// extra 4 is added for ghost cell
		Ni = UpperCoordinates.size()-1+4; 
		/**\param [in] Ni Number of cells(Including ghost cells) in "i" direction.*/ 
		Nj = N+4 ;  
		/**\param [in] Nj Number of cells(Including ghost cells) in "j" direction.*/ 
		Nk = 2+4 ; 
		/**\param [in] Nk Number of cells(Including ghost cells) in "k" direction.*/
		/* Here Nk = 2+4(Don't reduce it further) because this is 2D-simulation so 
		no need to take large number of cells in z direction */ 

		// Resizing the vectors
		Coordinate.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		iFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		jFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		kFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		CellVolume.resize(Ni+1,Dim2(Nj+1,Dim1(Nk+1)));
		ds.resize(Ni+1,Dim2(Nj+1,Dim1(Nk+1)));
		
		
		deltaz = finddeltaz(DownCoordinates);
		// First defining the grid points
		for (int i =2; i < Ni+1-2; i++) 
		// Will extend remaining x dir'n  after 
		{
			for (int  j=2;  j < Nj+1-2; j++)  
			//Will extend remaining y dir'n  after
			{
				for (int  k=2;  k < Nk+1-2; k++)
				// for (int  k=2;  k < Nk+1-2; k++)
				{
					Coordinate[i][j][k][0] = DownCoordinates[i+1-3][0] ;   
					Coordinate[i][j][k][1] = (j-2)*(UpperCoordinates[i+1-3][1]- 
						DownCoordinates[i+1-3][1])/N ;
					Coordinate[i][j][k][2] = (k-2)*deltaz;
				}
			}	
		}
		#endif
	}
	else if(GeometryOption == 4) // some other geometry
	{
		N = 11 ;
		Ni = 3*N+4 ;
		Nj = N+4 ;  
		Nk = 2+4 ; 

		// Resizing the vectors
		Coordinate.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		iFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		jFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		kFaceAreaVector.resize(Ni+1,Dim3(Nj+1,Dim2(Nk+1,Dim1(3))));
		CellVolume.resize(Ni+1,Dim2(Nj+1,Dim1(Nk+1)));
		ds.resize(Ni+1,Dim2(Nj+1,Dim1(Nk+1)));

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

	#if 1
	// here comes the live cells area vectors
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
				// Initial/Inlet j cells
				l0 = Coordinate[i][2][k][0];
				m0 = Coordinate[i][2][k][1];

				l1 = Coordinate[i+1][2][k][0];
				m1 = Coordinate[i+1][2][k][1];

				x0 = Coordinate[i][4-j][k][0];
				y0 = Coordinate[i][4-j][k][1];

				x1 = Coordinate[i+1][4-j][k][0];
				y1 = Coordinate[i+1][4-j][k][1];

				takeMirror(rx0,ry0,l0,m0,l1,m1,x0,y0);
				takeMirror(rx1,ry1,l0,m0,l1,m1,x1,y1);

				jFaceAreaVector[i][j][k][0] = -deltaz*(ry1-ry0) ;
				jFaceAreaVector[i][j][k][1] =  deltaz*(rx1-rx0) ;
				jFaceAreaVector[i][j][k][2] = 0;

				// Final/Exit j cells
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
	
#endif

	/** @brief Structure of grid out put file ("grids_Nozzle_2D.csv") 
	* - First line of the grid file will contain grid points
	*(excluding ghost cells) in x and y direction 
	* - This will exclude the ghost, only live cells or actual geomatry points
	*/
	ofstream kullu_grid ;
	kullu_grid.open("grids_2D.csv");
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


	// assigning the vector pointers 
	iFaceAreaVectorIn = iFaceAreaVector;
	jFaceAreaVectorIn = jFaceAreaVector;
	kFaceAreaVectorIn = kFaceAreaVector;
	CellVolumeIn = CellVolume;
	dsIn = ds;

}