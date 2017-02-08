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
double find_y(double x, std::vector<std::vector<double> > UpCods){
	int i=0;
	while(x>UpCods[i][0]){
		i++;
	}
	return UpCods[i-1][1] + (UpCods[i][1] - UpCods[i-1][1])*(x-UpCods[i-1][0])/(UpCods[i][0] - UpCods[i-1][0]);
}
// taking average of all dx for dz
double find_dz(std::vector<std::vector<double> > DownCodsNew)
{	
	int size = DownCodsNew.size();
	double dz = 0;
	for (int i = 1; i < size; ++i)
	{
		dz = dz + DownCodsNew[i][0] - DownCodsNew[i-1][0];  
	}
	return dz/(size-1);
} 

// this class calculates the area and volues in the domain

// Function defines the area vector and cell volumes 
// int main()
void grid( vector<vector<vector<vector<double> > > > & x_face_area_in,
		  vector<vector<vector<vector<double> > > > & y_face_area_in,
		   vector<vector<vector<vector<double> > > > & z_face_area_in,
		  vector<vector<vector<double> > > & cell_volume_in,
		   int & Nx, int & Ny, int & Nz)
{
	std::vector<std::vector<double> > UpCods;
	std::vector<std::vector<double> > DownCods;
   	ifstream nozzleData ("NozzleGeomData.csv");
	int j = 0;
	   
   while(!nozzleData.eof()){

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
	   UpCods.push_back(temp);

	   temp[1] = 0.0; // change the y only and push it to the Down vector
	   DownCods.push_back(temp);
	   ++j;
   }
   // rendomaly extra zeros at the end so to remove them pop is used
   UpCods.pop_back();
   DownCods.pop_back();

   nozzleData.close();
   
#if 1
   	std::vector<std::vector<double> > UpCodsNew;
	std::vector<std::vector<double> > DownCodsNew;

	UpCodsNew.push_back(UpCods[0]); // starting point is same
	DownCodsNew.push_back(DownCods[0]);
   	

   double dx;
   double dy;
   double x = UpCods[0][0];
   double y = UpCods[0][1];
   int N = 20 ; // total N+1 points after including the boundary points because N cells
   int i = 0 ;

   while(UpCods[UpCods.size()-1][0]>x)
   {
	   std::vector<double> xyup;
	   std::vector<double> xydown;
	   // cout << UpCodsNew[i][0] << "   " << UpCodsNew[i][1] << "   " << DownCodsNew[i][0] << "   " << DownCodsNew[i][1] << endl;

	   dy = (UpCodsNew[i][1] - DownCodsNew[i][1])/N ;
	   dx = dy ;
	   
	   x = x + dx ;

	   if(UpCods[UpCods.size()-1][0]>x){
		   y = find_y(x, UpCods);
		   xyup.push_back(x);
		   xyup.push_back(y);
		   UpCodsNew.push_back(xyup);
		
		   xydown.push_back(x);
		   xydown.push_back(0.0);
			   
 		   DownCodsNew.push_back(xydown);
		   i++;
	   }
   }

#endif


#if 1
// Grids for hypersonic nozzle (this is the second case which will test the scheme)
	
	// extra 4 is added for ghost cell
	Nx = UpCodsNew.size()-1+4; // Total cells in X-dir  
	Ny = N+4 ;  // Total cell in Y-dir 
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

	double dz = find_dz(DownCodsNew);
	// First defining the grid points
	for (int i =2; i < Nx+1-2; i++) // Will extend remaining x dir'n  after 
	{
		for (int  j=2;  j < Ny+1-2; j++)  //Will extend remaining y dir'n  after
		{
			for (int  k=0;  k < Nz+1; k++)
			{
				grid_point[i][j][k][0] = DownCodsNew[i+1-3][0] ;   
				grid_point[i][j][k][1] = (j-2)*(UpCodsNew[i+1-3][1]- DownCodsNew[i+1-3][1])/N ;   
				grid_point[i][j][k][2] = k*dz;
			}
		}	
	}

	// y right ghost cell
	for (int i =2; i < Nx+1-2; ++i) 
	{
		for (int j = 0; j < 2; ++j)
		{
			for (int  k=0;  k < Nz+1; k++)
			{
				grid_point[i][1-j][k][0] = grid_point[i][2-j][k][0] ;   
				grid_point[i][1-j][k][1] = grid_point[i][2-j][k][1] - dz ;
				grid_point[i][1-j][k][2] = grid_point[i][2-j][k][2] ;
			}
		}
	}
	//y left ghost cell
	for (int i =2; i < Nx+1-2; ++i) 
	{
		for (int j = Ny+1-2; j < Ny+1; ++j)
		{
			for (int  k=0;  k < Nz+1; k++)
			{
				grid_point[i][j][k][0] = grid_point[i][j-1][k][0] ;   
				grid_point[i][j][k][1] = grid_point[i][j-1][k][1] + dz ;
				grid_point[i][j][k][2] = grid_point[i][j-1][k][2] ;
			}
		}
	}

	// x right ghost cell
	for (int i =0; i < 2; ++i) 
	{
		for (int j = 0; j < Ny+1; ++j)
		{
			for (int  k=0;  k < Nz+1; k++)
			{
				grid_point[1-i][j][k][0] = grid_point[2-i][j][k][0] - dz ;   
				grid_point[1-i][j][k][1] = grid_point[2-i][j][k][1] ;
				grid_point[1-i][j][k][2] = grid_point[2-i][j][k][2] ;
			}
		}
	}
	//x left ghost cell
	for (int i = Nx+1-2; i < Nx+1; ++i) 
	{
		for (int j = 0; j < Ny+1; ++j)
		{
			for (int  k=0;  k < Nz+1; k++)
			{
				grid_point[i][j][k][0] = grid_point[i-1][j][k][0] + dz ;   
				grid_point[i][j][k][1] = grid_point[i-1][j][k][1] ;
				grid_point[i][j][k][2] = grid_point[i-1][j][k][2] ;
			}
		}
	}
	
#endif	
#if 1
	// here comes the area vectors
	for (int i = 0; i  < Nx; ++i)
	{
		for (int  j = 0;  j < Ny; ++j)
		{
			for (int  k = 0;  k < Nz; ++k)
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

	// cell volume 
	for (int i = 0; i  < Nx; ++i)
	{
		for (int  j= 0;  j < Ny; ++j)
		{
			for (int  k= 0;  k < Nz; ++k)
			{
				cell_volume[i][j][k] = 0.5 * (grid_point[i+1][j][k][0] - grid_point[i][j][k][0]) * 
				( (grid_point[i][j+1][k][1] - grid_point[i][j][k][1])  + 
				(grid_point[i+1][j+1][k][1]-grid_point[i+1][j][k][1]) ) * dz ; 
			}
		}	
	}

#endif

ofstream kullu_grid ;
	kullu_grid.open("grids_Nozzle_2D.csv");
	for (int i = 0; i < Nx+1; ++i)
	{
		for (int j = 0; j < Ny+1; ++j)
		{
			kullu_grid << grid_point[i][j][4][0] << "," << grid_point[i][j][4][1] << endl ; 
		}
	}   
   // return 0;
// assigning the vector 
	x_face_area_in = x_face_area;
	y_face_area_in = y_face_area;
	z_face_area_in = z_face_area;
	cell_volume_in = cell_volume;

}