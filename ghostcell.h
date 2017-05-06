/*! \file  ghostcell.h
    \brief This header file functions find the ghost cells area 
    vectors and ghost cell volumes
    \author Kuldeep Singh
    \date 2017
*/
#include <iostream>
#include "math.h"
#include <fstream> /* For file handling */
#include <string> /* For strings */
#include <vector> /* For vectors*/
#include <cstdlib> /* For converting string into numerical value */

using namespace std;
void ghostcell(
	vector<vector<vector<vector<double> > > > Coordinates,
	vector<vector<vector<vector<double> > > > iFaceAreaVector,
	vector<vector<vector<vector<double> > > > jFaceAreaVector,
	vector<vector<vector<vector<double> > > > kFaceAreaVector,
	vector<vector<vector<double> > > CellVolume,

	vector<vector<vector<double> > > & i0GhostCellVolume,
	vector<vector<vector<double> > > & j0GhostCellVolume,
	vector<vector<vector<double> > > & k0GhostCellVolume,

	vector<vector<vector<double> > > & iNiGhostCellVolume,
	vector<vector<vector<double> > > & jNjGhostCellVolume,
	vector<vector<vector<double> > > & kNkGhostCellVolume,
	int Ni, int Nj, int Nk)
{
	for (int j = 0; j < Nj; ++j)
	{
		for (int k = 0; k < Nk; ++k)
		{
			i0GhostCellVolume[0][j][k] = CellVolume[0][j][k];
			iNiGhostCellVolume[0][j][k] = CellVolume[Ni-1][j][k];
		}
	}

	for (int i = 0; i < Ni; ++i)
	{
		for (int k = 0; k < Nk; ++k)
		{
			j0GhostCellVolume[i][0][k] = CellVolume[i][0][k];
			jNjGhostCellVolume[i][0][k] = CellVolume[i][Nj-1][k];
		}
	}

	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			k0GhostCellVolume[i][j][0] = CellVolume[i][j][0];
			kNkGhostCellVolume[i][j][0] = CellVolume[i][j][Nk-1];
		}
	}
}
	#if 0
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
				l0 = Coordinates[i][2][k][0];
				m0 = Coordinates[i][2][k][1];

				// M2 point  
				l1 = Coordinates[i+1][2][k][0];
				m1 = Coordinates[i+1][2][k][1];
				//Mirror has been taken about the line M1 M2  


				x0 = Coordinates[i][4-j][k][0];
				y0 = Coordinates[i][4-j][k][1];

				x1 = Coordinates[i+1][4-j][k][0];
				y1 = Coordinates[i+1][4-j][k][1];

				takeMirror(rx0,ry0,l0,m0,l1,m1,x0,y0);
				takeMirror(rx1,ry1,l0,m0,l1,m1,x1,y1);

				jFaceAreaVector[i][j][k][0] = -deltaz*(ry1-ry0) ;
				jFaceAreaVector[i][j][k][1] =  deltaz*(rx1-rx0) ;
				jFaceAreaVector[i][j][k][2] = 0;
				
				// // Assigning the mirrored points to the coordinates but this is 
				// //for the testing purposes only
				// Coordinates[i][j][k][0] = rx0;
				// Coordinates[i][j][k][1] = ry0;
				// Coordinates[i][j][k][2] = Coordinates[i][4-j][k][2];

				// Coordinates[i+1][j][k][0] = rx1;
				// Coordinates[i+1][j][k][1] = ry1;
				// Coordinates[i+1][j][k][2] = Coordinates[i+1][4-j][k][2];
				
				// Top ghost cells
				l0 = Coordinates[i][Nj-2][k][0];
				m0 = Coordinates[i][Nj-2][k][1];

				l1 = Coordinates[i+1][Nj-2][k][0];
				m1 = Coordinates[i+1][Nj-2][k][1];

				x0 = Coordinates[i][Nj-3-j][k][0];
				y0 = Coordinates[i][Nj-3-j][k][1];

				x1 = Coordinates[i+1][Nj-3-j][k][0];
				y1 = Coordinates[i+1][Nj-3-j][k][1];

				takeMirror(rx0,ry0,l0,m0,l1,m1,x0,y0);
				takeMirror(rx1,ry1,l0,m0,l1,m1,x1,y1);

				jFaceAreaVector[i][Nj-2+1+j][k][0] = -deltaz*(ry1-ry0) ;
				jFaceAreaVector[i][Nj-2+1+j][k][1] =  deltaz*(rx1-rx0) ;
				jFaceAreaVector[i][Nj-2+1+j][k][2] = 0;

				CellVolume[i][j][k] = CellVolume[i][3-j][k];
				CellVolume[i][Nj-2+j][k] = CellVolume[i][Nj-3-j][k];

				// // Assigning the mirrored points to the coordinates but this is 
				// //for the testing purposes only
				// Coordinates[i][j][k][0] = rx0;
				// Coordinates[i][j][k][1] = ry0;
				// Coordinates[i][j][k][2] = Coordinates[i][Nj-3-j][k][2];

				// Coordinates[i+1][j][k][0] = rx1;
				// Coordinates[i+1][j][k][1] = ry1;
				// Coordinates[i+1][j][k][2] = Coordinates[i+1][Nj-3-j][k][2];
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
	#endif