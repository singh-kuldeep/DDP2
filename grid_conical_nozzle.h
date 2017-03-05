#include "math.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;
// using std::cout;
// using std::endl;
// using std::ifstream;


double distance(std::vector<double> v1, std::vector<double> v2)
{
	double dl = sqrt(pow((v2[0]-v1[0]),2) + pow((v2[1]-v1[1]),2) + pow((v2[2]-v1[2]),2));
	cout << dl << endl;
	return dl;
}


// this function will take boundary cell points and will calculate the ghost cell grid ponits by taking the mirror image
void mirror(double &x,double&y, double a,double b,double c,double d,double l,double m)
{
	double slop = (d-b)/(c-a);
	y = (slop*(2*l-2*a)+m*slop*slop+2*b-m)/(1+slop*slop);
	x =  (m-y)*slop+l;
} 
// this class calculates the area and volues in the domain

// Function defines the area vector and cell volumes 
// int main()
void grid( vector<vector<vector<vector<double> > > > & x_face_area_in,
		  vector<vector<vector<vector<double> > > > & y_face_area_in,
		   vector<vector<vector<vector<double> > > > & z_face_area_in,
		  vector<vector<vector<double> > > & cell_volume_in,
		  vector<vector<vector<double> > > & delta_s_in,
		   int & Nx, int & Ny, int & Nz)
{
	
#if 1
// Grids for hypersonic nozzle (this is the second case which will test the scheme)
	
   	int N = 15; // Keep is odd for symmatry 

	double angle = 60.0;
   	double corner_2_y = 50.0;
   	double corner_1_y = 1.0;
   	double corner_1_x = 0.0;
   	double corner_2_x = corner_1_x + (corner_2_y- corner_1_y);
	double x = corner_1_x;
	double h = 2*corner_1_y; // height at x=corner_1_x
	double dy = h/N;
	double dx = dy;
	double dz = 10*dx;
	
	// extra 4 is added for ghost cell
	Ny = N+4 ;  // Total cell in Y-dir 
	Nx = N+4; // Total cells in X-dir  
	Nz = 1+4 ; // Because this is 2D-simulation so no need to take large number of grids in z direction 

	// Creating a 4D vector object for grid points
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> matrix4D;

	// this store previous values of variables (density , three momentum, energy)
	matrix4D grid_point(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	// matrix4D grid_point_rot(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	matrix4D x_face_area(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	matrix4D y_face_area(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	matrix4D z_face_area(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 

	Dim3 cell_volume(Nx,Dim2(Ny,Dim1(Nz)));
	Dim3 delta_s(Nx,Dim2(Ny,Dim1(Nz)));

	// First defining the grid points
	for (int i =2; i < Nx+1-2; i++) // Will extend remaining x dir'n  after 
	{
		for (int  j=2;  j < Ny+1-2; j++)  //Will extend remaining y dir'n  after
		{
			for (int  k=0;  k < Nz+1; k++)
			{
				grid_point[i][j][k][0] = x;   
				grid_point[i][j][k][1] = dy * (j - Ny/2);   
				grid_point[i][j][k][2] = (k-2)*dz;
			}
		}
		x = x+dy;
		// h = 2*corner_2_y;
		// h = 2*(corner_1_y+((corner_2_y- corner_1_y)/(corner_2_x - corner_1_x))*(x- corner_1_x));
		dy = dy*1.3;
		// dx = dy;
	}

#endif	

#if 1
	// here comes the live cells area vectors
	for (int i = 2; i  < Nx; ++i)
	{
		for (int  j = 2;  j < Ny; ++j)
		{
			for (int  k = 2;  k < Nz; ++k)
			{
				x_face_area[i][j][k][0] = (grid_point[i][j+1][k][1]-grid_point[i][j][k][1])*dz ;
				x_face_area[i][j][k][1] = 0 ;
				x_face_area[i][j][k][2] = 0 ;

				y_face_area[i][j][k][0] = -dz*(grid_point[i+1][j][k][1]-grid_point[i][j][k][1]) ;
				y_face_area[i][j][k][1] =  dz*(grid_point[i+1][j][k][0]-grid_point[i][j][k][0]) ;
				y_face_area[i][j][k][2] = 0 ;

				z_face_area[i][j][k][0] = 0 ; 
				z_face_area[i][j][k][1] = 0 ;
				z_face_area[i][j][k][2] = 0.5*(grid_point[i+1][j][k][0] - grid_point[i][j][k][0])*
				( (grid_point[i][j+1][k][1]-grid_point[i][j][k][1]) + 
					(grid_point[i+1][j+1][k][1]-grid_point[i+1][j][k][1]) ); 
			}
		}	
	}

	// live cell volumes 
	for (int i = 2; i  < Nx-2; ++i)
	{
		for (int  j= 2;  j < Ny-2; ++j)
		{
			for (int  k= 2;  k < Nz-2; ++k)
			{
				cell_volume[i][j][k] = 0.5 * (grid_point[i+1][j][k][0] - grid_point[i][j][k][0]) * 
				( (grid_point[i][j+1][k][1] - grid_point[i][j][k][1])  + 
				(grid_point[i+1][j+1][k][1]-grid_point[i+1][j][k][1]) ) * dz ; 
			}
		}	
	}

	// here comes the x_ghost cells area vectors and volumes
	for(int i=0; i<2; ++i)
	{
		for (int j = 2; j < Ny+1-2; ++j)
		{
			for (int k = 2; k < Nz+1-2; ++k)
			{
				x_face_area[i][j][k][0] = x_face_area[4-i][j][k][0];
				x_face_area[i][j][k][1] = 0 ;
				x_face_area[i][j][k][2] = 0 ;

				x_face_area[Nx-1+i][j][k][0] = x_face_area[Nx-3-i][j][k][0];
				x_face_area[Nx-1+i][j][k][1] = 0 ;
				x_face_area[Nx-1+i][j][k][2] = 0 ;

				cell_volume[i][j][k] = cell_volume[3-i][j][k];
				cell_volume[Nx-2+i][j][k] = cell_volume[Nx-3-i][j][k];
			}
		}
	}	

	double x0,y0,x1,y1; // y_grid before reflection or to be reflected
	double l0,m0,l1,m1; // line about which reflection needs to be done
	double rx0,ry0,rx1,ry1; // y_grid after reflection 
	// here comes the y_ghost cells area vectors
	for (int i = 2; i < Nx+1-2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			for (int k = 2; k < Nz+1-2; ++k)
			{
				l0 = grid_point[i][2][k][0];
				m0 = grid_point[i][2][k][1];

				l1 = grid_point[i+1][2][k][0];
				m1 = grid_point[i+1][2][k][1];

				x0 = grid_point[i][4-j][k][0];
				y0 = grid_point[i][4-j][k][1];

				x1 = grid_point[i+1][4-j][k][0];
				y1 = grid_point[i+1][4-j][k][1];

				mirror(rx0,ry0,l0,m0,l1,m1,x0,y0);
				mirror(rx1,ry1,l0,m0,l1,m1,x1,y1);

				y_face_area[i][j][k][0] = -dz*(ry1-ry0) ;
				y_face_area[i][j][k][1] =  dz*(rx1-rx0) ;
				y_face_area[i][j][k][2] = 0;

				l0 = grid_point[i][Ny-2][k][0];
				m0 = grid_point[i][Ny-2][k][1];

				l1 = grid_point[i+1][Ny-2][k][0];
				m1 = grid_point[i+1][Ny-2][k][1];

				x0 = grid_point[i][Ny-3-j][k][0];
				y0 = grid_point[i][Ny-3-j][k][1];

				x1 = grid_point[i+1][Ny-3-j][k][0];
				y1 = grid_point[i+1][Ny-3-j][k][1];

				mirror(rx0,ry0,l0,m0,l1,m1,x0,y0);
				mirror(rx1,ry1,l0,m0,l1,m1,x1,y1);

				y_face_area[i][Ny-2+1+j][k][0] = -dz*(ry1-ry0) ;
				y_face_area[i][Ny-2+1+j][k][1] =  dz*(rx1-rx0) ;
				y_face_area[i][Ny-2+1+j][k][2] = 0;

				cell_volume[i][j][k] = cell_volume[i][3-j][k];
				cell_volume[i][Ny-2+j][k] = cell_volume[i][Ny-3-j][k];
			}
		}
	}


	// here comes the z_ghost cells area vectors
	for(int i=2; i<Nx+1-2; ++i)
	{
		for (int j = 2; j < Ny+1-2; ++j)
		{
			for (int k = 0; k < 2; ++k)
			{
				z_face_area[i][j][k][0] = 0;
				z_face_area[i][j][k][1] = 0 ;
				z_face_area[i][j][k][2] = z_face_area[i][j][4-k][2] ;

				z_face_area[i][j][k][0] = 0;
				z_face_area[i][j][k][1] = 0 ;
				z_face_area[i][j][Nz-1+k][2] = z_face_area[i][j][Nz-3-k][2] ;

				cell_volume[i][j][k] = cell_volume[i][j][3-k];
				cell_volume[i][j][Nz-2+k] = cell_volume[i][j][Nz-3-k];
			}
		}
	}	
	

	
	#if 1
	for (int i = 1; i  < Nx-2; ++i)
	{
		for (int  j= 1;  j < Ny-2; ++j)
		{
			for (int  k= 1;  k < Nz-2; ++k)
			{
				// delta_s[i][j][k] = min(distance(&grid_point[i][j][k],&grid_point[i+1][j][k]),distance(&grid_point[i][j+1][k],&grid_point[i+1][j+1][k]),
				// 	distance(&grid_point[i][j][k+1],&grid_point[i+1][j][k+1]),distance(&grid_point[i][j+1][k+1],&grid_point[i+1][j+1][k+1]),
				// 	distance(&grid_point[i][j][k],&grid_point[i][j+1][k]),distance(&grid_point[i+1][j][k],&grid_point[i+1][j+1][k]),
				// 	distance(&grid_point[i][j][k+1],&grid_point[i][j+1][k+1]),distance(&grid_point[i+1][j][k+1],&grid_point[i+1][j+1][k+1]),
				// 	distance(&grid_point[i][j][k],&grid_point[i][j][k+1]),distance(&grid_point[i+1][j][k],&grid_point[i+1][j][k+1]),
				// 	distance(&grid_point[i][j+1][k],&grid_point[i][j+1][k+1]),distance(&grid_point[i+1][j+1][k],&grid_point[i+1][j+1][k+1]));

				delta_s[i][j][k] = grid_point[i+1][Ny/2][Nz/2][1] - grid_point[i+1][Ny/2 - 1][Nz/2][1];	
			}
		}	
	}
	
	/////////////////////////////////////////////////////////////////////////////
	// writeing delta_s into the file 
	/////////////////////////////////////////////////////////////////////////////	
		ofstream kullu_ds ;
		kullu_ds.open("Nozzle_ds.csv");
		for (int i = 2; i < Nx-2; ++i)
		{
			// for (int j = 2; j < Ny-2; ++j)
			// {
				kullu_ds << delta_s[i][Ny/2][Nz/2] << endl ; 
			// }
		}
	#endif 
	
#endif
/////////////////////////////////////////////////////////////////////////////
// STUCTURE OF GRID FILE	
// 1. First line of the grid file will contain grid points(excluding ghost cells) in x and y direction 
// 2. This will exclude the ghost, only live cells or actual geomatry points
// 3.  	
/////////////////////////////////////////////////////////////////////////////	
	ofstream kullu_grid ;
	kullu_grid.open("grids_Nozzle_2D.csv");
	kullu_grid << Nx-4 << "," << Ny-4 << endl ; 
	for (int i = 2; i < Nx-2; ++i)
	{
		for (int j = 2; j < Ny-2; ++j)
		{
			// kullu_grid << 0.5*(grid_point[i][j][4][0]+grid_point[i+1][j][4][0]) << "," 
			// << 0.5*(grid_point[i][j][4][1]+grid_point[i][j+1][4][1]) << endl ;
			kullu_grid <<  grid_point[i][j][4][0] << ","<< grid_point[i][j][4][1] << endl;
		}
	}   


// assigning the vector 
	x_face_area_in = x_face_area;
	y_face_area_in = y_face_area;
	z_face_area_in = z_face_area;
	cell_volume_in = cell_volume;
	delta_s_in = delta_s;

}