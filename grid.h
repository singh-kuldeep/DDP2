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
	Nj = 1 ;  
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

		N = 20 ;
		Ni = 1.5*N ;
		Nj = N ;  
		Nk = 20 ; 

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
				CellCenter << (Coordinate[i][j][(Nk+1)/2][0] + Coordinate[i+1][j][(Nk+1)/2][0] + 
							Coordinate[i][j+1][(Nk+1)/2][0] + Coordinate[i+1][j+1][(Nk+1)/2][0])/4 
				<< ","<< (Coordinate[i][j][(Nk+1)/2][1] + Coordinate[i+1][j][(Nk+1)/2][1] + 
							Coordinate[i][j+1][(Nk+1)/2][1] + Coordinate[i+1][j+1][(Nk+1)/2][1])/4
				<< endl;
			}
		}   
	}
	cout << "Grid genration succsesful :)  " << endl ;
}