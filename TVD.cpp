// assuming that only these three things are given 
// (1)Grid points
// (2)Boundary condition
// (3)Initial condition are given
#include "iostream"
#include <vector>
#include <fstream>
#include "math.h"
#include "time.h"
#include "netfluxinterface.h"
#include "BC.h"
// #include "grid_straight_duct.h"
// #include "grid_bump.h"
#include "grid_nozzle.h"

#define gamma 1.4
#define gasconstant 287.14
#define heatcapacityconstantvolume 717.5

using namespace std ;
// BC function implements the boundary condition 
// x_face_area is not required Because flow in x direction and 2D flow
void BC(vector<vector<vector<vector<double> > > > & variablesvector, vector<vector<vector<vector<double> > > > &
 y_face_area, vector<vector<vector<vector<double> > > > & z_face_area, int Nx, int Ny, int Nz) ;

// grid function genrates the area vector and volume for the all cells in the domain
void grid(vector<vector<vector<vector<double> > > > & x_face_area , 
	vector<vector<vector<vector<double> > > > & y_face_area,
	vector<vector<vector<vector<double> > > > & z_face_area, 
	vector<vector<vector<double> > >& cell_volume,
	int & Nx, int & Ny, int & Nz) ;

// double dt(int CFL, ); // time step at every iteration 

int main ()
{
	time_t start, end ;
	time(&start); // noteing the starting time

	double deltat = 0.0000015; // this is for CFL = 0.2
	double TIME = 1000000*deltat;
	int totaltimesteps = floor(TIME/deltat) ;

	double lenght = 26 ; // keep it even
	double delta = 1.0 ; // this basically defines the grid size

	int N = floor(lenght/delta);
	// extra 4 is added for ghost cell
	int Nx = (3*N + 4);
	int Ny = N+4;
	int Nz = 1+4; // Because this is 2D-simulation so no need to take large number of grids in z direction 

	// Creating a 4D vector object
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> matrix4D;

	
	// this store grid information in the domain
	// matrix4D x_face_area.resize(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	// matrix4D y_face_area(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	// matrix4D z_face_area(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	// Dim3 cell_volume(Nx,Dim2(Ny,Dim1(Nz)));

	std::vector<std::vector<std::vector<std::vector<double> > > > x_face_area ;
	std::vector<std::vector<std::vector<std::vector<double> > > > y_face_area ;
	std::vector<std::vector<std::vector<std::vector<double> > > > z_face_area ;
	std::vector<std::vector<std::vector<double> > > cell_volume ;
	// Now the grid function has defined the area vector and volume od cells 
	// grid function will change the grid points as well 
	grid(x_face_area,y_face_area,z_face_area,cell_volume,Nx,Ny,Nz);

	cout << Nx<< "  " << Ny << "  " << Nz << endl;

	// this store previous values of variables (density , three momentum, energy)
	matrix4D variablesvector(Nx,Dim3(Ny,Dim2(Nz,Dim1(5)))); 
	// this store new values of variables(density , three momentum, energy) 
	matrix4D variablevectornew(Nx,Dim3(Ny,Dim2(Nz,Dim1(5)))); 
	// Initial condition(these are just rendom values )
	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			for (int k = 0; k < Nz; ++k)
			{
				variablesvector[i][j][k][0] = 1.16;
				variablesvector[i][j][k][1] = 0 ;
				variablesvector[i][j][k][2] = 0 ;
				variablesvector[i][j][k][3] = 0 ;
				variablesvector[i][j][k][4] = 271767;

				variablevectornew[i][j][k][0] = 1.16;
				variablevectornew[i][j][k][1] = 0 ;
				variablevectornew[i][j][k][2] = 0 ;
				variablevectornew[i][j][k][3] = 0 ;
				variablevectornew[i][j][k][4] = 271767;
			}
		}
	}

	// time progression
	// this file is opend to store the mass residual at each time step
	ofstream kullu_mass ;
	kullu_mass.open("Residual_Nozzle.dat");
	// kullu_mass <<  "t(secs)" << "," << "density_res"  << "," << "x_momentum_res" << "," <<
	// 	 "y_momentum_res" <<","<< "z_momentum_res" << "," << "energy_res" << endl ;

	for (int t = 0; t < totaltimesteps; ++t)
	{
		// cout << "timestep  = " << t << endl ; 
		// Befor every time step we need to have proper value in the ghost cells 
		// So applying the Boundary condition
		// BC takes  care of Inlet, Exit, y-wall and Z-wall boundary condition
		BC(variablesvector,y_face_area,z_face_area,Nx,Ny,Nz); 
		
		// next timestep variablesvectors calcuation
		for (int i = 1; i < Nx-2; ++i) // Or (int i = 2; i < (Nx+1)-2; ++i) Total Nx+1 interface, 2 used in BC implementation
		{
			for (int j = 1; j < Ny-2; ++j)
			{
				for (int k = 1; k < Nz-2; ++k)
				{
					//x right interface volume
					double xcellinterfacevolumeright = 0.5*(cell_volume[i][j][k] + cell_volume[i+1][j][k]);

					//y right interface volume
					double ycellinterfacevolumeright = 0.5*(cell_volume[i][j][k] + cell_volume[i][j+1][k]);

					//z right interface volume
					double zcellinterfacevolumeright = 0.5*( cell_volume[i][j][k] + cell_volume[i][j][k+1]);

					// net flux using the class netfluxinterface
					netfluxinterface xrightface(variablesvector[i-1][j][k],variablesvector[i][j][k],variablesvector[i+1][j][k],
							variablesvector[i+2][j][k] ,x_face_area[i-1][j][k] ,x_face_area[i][j][k] ,x_face_area[i+1][j][k],
						cell_volume[i-1][j][k], cell_volume[i][j][k], cell_volume[i+1][j][k], cell_volume[i+2][j][k] ,deltat) ;

					netfluxinterface yrightface(variablesvector[i][j-1][k],variablesvector[i][j][k],variablesvector[i][j+1][k],
							variablesvector[i][j+2][k] ,y_face_area[i][j-1][k] ,y_face_area[i][j][k] ,y_face_area[i][j+1][k],
						cell_volume[i][j-1][k], cell_volume[i][j][k], cell_volume[i][j+1][k], cell_volume[i][j+2][k] ,deltat) ;

					netfluxinterface zrightface(variablesvector[i][j][k-1],variablesvector[i][j][k],variablesvector[i][j][k+1],
							variablesvector[i][j][k+2] ,z_face_area[i][j][k-1] ,z_face_area[i][j][k] ,z_face_area[i][j][k+1],
						cell_volume[i][j][k-1], cell_volume[i][j][k], cell_volume[i][j][k+1], cell_volume[i][j][k+2] ,deltat) ;

					// updating the variablevectornew using flux at the right interfaces
					for (int l = 0; l < 5; ++l)
					{
						variablevectornew[i][j][k][l] -=(deltat/xcellinterfacevolumeright)*(xrightface.netflux[l]);
						variablevectornew[i+1][j][k][l] +=(deltat/xcellinterfacevolumeright)*(xrightface.netflux[l]);

						variablevectornew[i][j][k][l] -=(deltat/ycellinterfacevolumeright)*(yrightface.netflux[l]);
						variablevectornew[i][j+1][k][l] +=(deltat/ycellinterfacevolumeright)*(yrightface.netflux[l]);

						variablevectornew[i][j][k][l] -=(deltat/zcellinterfacevolumeright)*(zrightface.netflux[l]);
						variablevectornew[i][j][k+1][l] +=(deltat/zcellinterfacevolumeright)*(zrightface.netflux[l]);
					}
				}
			}
		}

		// Residual calculation after each timestep and writing that into the mass_residual file
		double density_res = 0.0 ;
		double x_momentum_res = 0.0 ;
		double y_momentum_res = 0.0 ;
		double z_momentum_res = 0.0 ;
		double energy_res = 0.0 ;

		int number_ponit = 0 ; 
		for (int x = 2; x < Nx-2; ++x)
			{
				for (int y = 2; y < Ny-2; ++y)
				{
					number_ponit   += 1 ; 
					density_res    += pow((variablevectornew[x][y][2][0] - variablesvector[x][y][2][0]),2);
					x_momentum_res += pow((variablevectornew[x][y][2][1] - variablesvector[x][y][2][1]),2);     
					y_momentum_res += pow((variablevectornew[x][y][2][2] - variablesvector[x][y][2][2]),2);     
					z_momentum_res += pow((variablevectornew[x][y][2][3] - variablesvector[x][y][2][3]),2);     
					energy_res     += pow((variablevectornew[x][y][2][4] - variablesvector[x][y][2][4]),2);     
				}
			}
		// cout << "number_ponit" << number_ponit << endl ;
		kullu_mass << t*deltat << "," << sqrt(density_res/1875)  << "," << sqrt(x_momentum_res/1875) << "," <<
		 sqrt(y_momentum_res/1875) <<","<< sqrt(z_momentum_res/1875) << "," << sqrt(energy_res/1875) << endl ;
		
		// cout << "timestep     "<<  t << "    Residual    "<< density_res   << endl ;
		cout <<  t << "  -  "<< density_res << endl ;

		// before going to the new timestep update variablesvector by variablevectornew
		for (int i = 2; i < Nx-2; ++i)
		{
			for (int j = 2; j < Ny-2; ++j)
			{
				for (int k = 2; k < Nz-2; ++k)
				{
					for (int l = 0; l < 5; ++l)
					{
						variablesvector[i][j][k][l] = variablevectornew[i][j][k][l] ;
					}				
				}
			}
		}


		if (t%50 == 0)
		{
			// storing the velocity in one plane
			ofstream kullu_2D ;
			kullu_2D.open("2D_parameters_Nozzle.dat");
			// kullu_2D << "density" << "," << "density*u" << ","<< "density*v" << "," << "density*w" << "," << "energy"  << endl ;
			for (int i = 2; i < Nx-2; ++i)
			{
				for (int j = 2; j < Ny-2; ++j)
				{
					kullu_2D << variablesvector[i][j][2][0] << "," << variablesvector[i][j][2][1] <<","<< 
					variablesvector[i][j][2][2] << "," << variablesvector[i][j][2][3] << "," << variablesvector[i][j][2][4] << endl ;
				}
			}
		}
	} 
	// time progression ends here

	time(&end) ;
	double diff = difftime (end,start);
	cout << "Time taken by the solver in secs = " << diff << endl ;
	return 0;
}

