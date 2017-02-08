#include "iostream"
#include "math.h"
#include <vector>

#define gamma 1.4
#define gasconstant 287.14
#define heatcapacityconstantvolume 717.5

// vector<double> normal(vector<double> & v)
void normal(vector<double> & v)
{
	vector<double> n(3);
	double mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

	v[0] = v[0]/mag;
	v[1] = v[1]/mag;
	v[2] = v[2]/mag;
	// return n;
}
// Function Implements the boundary condition
void BC(vector<vector<vector<vector<double> > > > & variablesvector , vector<vector<vector<vector<double> > > > & y_face_area,
	vector<vector<vector<vector<double> > > > & z_face_area, int Nx, int Ny, int Nz)
{
	// inilet conditions (user given data)
	// one has to maintion only inlet Mach number totalpressure and the totaltemperature
	double totaltemperature = 518.76 ;
	double totalpressure = 792766.8;
	double totaldensity = totalpressure / (gasconstant*totaltemperature) ;
	double theta = 3.14159265 * 0 / 180 ;// rotation
	 	
	// double theta = 3.14159265 * 25 / 180 ;// rotation 	

	// Inlet ghost cell updating usinig the totel quantities(T_0, P_0) and flow direction  
	for (int j =2; j < Ny-2; ++j)
	{
		for (int k =2; k < Nz-2; ++k)
		{
			double inletpressure = (gamma-1)*(variablesvector[2][j][k][4] - 0.5*(pow(variablesvector[2][j][k][1],2) +
				pow(variablesvector[2][j][k][2],2)+pow(variablesvector[2][j][k][3],2))/variablesvector[2][j][k][0]) ;

			double Mach = sqrt( (2/(gamma-1)) * ( pow( (totalpressure/inletpressure),((gamma-1)/gamma) ) -1 ) ) ;

			double inlettemperature = totaltemperature / (1+((gamma-1)*Mach*Mach)/2);

			double inletvelocity = Mach * sqrt(gamma*gasconstant*inlettemperature) ;

			double inletdensity = totaldensity / pow((totalpressure/inletpressure),(1/gamma)) ;

			variablesvector[0][j][k][0] = inletdensity ;
			variablesvector[0][j][k][1] = cos(theta)*inletdensity*inletvelocity ;
			variablesvector[0][j][k][2] = 0 ;
			variablesvector[0][j][k][3] = -sin(theta)*inletdensity*inletvelocity ;
			variablesvector[0][j][k][4] = inletpressure/(gamma-1) + 0.5*inletdensity*inletvelocity*inletvelocity ;

			variablesvector[1][j][k][0] = inletdensity ;
			variablesvector[1][j][k][1] = cos(theta)*inletdensity*inletvelocity ;
			variablesvector[1][j][k][2] = 0 ;
			variablesvector[1][j][k][3] = -sin(theta)*inletdensity*inletvelocity ;
			variablesvector[1][j][k][4] = inletpressure/(gamma-1) + 0.5*inletdensity*inletvelocity*inletvelocity ;
		}
	}

	// at exit updating the x ghost cell (this is true whene flow is supersonic)
	for (int j =2; j < Ny-2; ++j)
	{
		for (int k =2; k < Nz-2; ++k)
		{
			for (int l = 0; l < 5; ++l)
			{
				variablesvector[Nx-2][j][k][l] = variablesvector[Nx-3][j][k][l] ;
				variablesvector[Nx-1][j][k][l] = variablesvector[Nx-4][j][k][l] ;
			}
		}
	}

	// updating the wall ghost cell(Y - wall)
	for (int i = 2; i < Nx-2; ++i)
	{
		for (int k = 2; k < Nz-2; ++k)
		{
			double y_face_area_unit_vector[3] = {0.0,2.0,0.0} ;
			// double mag = sqrt( pow(y_face_area[i][2][k][0] ,2) + pow(y_face_area[i][2][k][1] ,2) + pow(y_face_area[i][2][k][2] ,2) ) ;
			//using leftcell(j=2) filling j=1
			y_face_area_unit_vector[0] = y_face_area[i][2][k][0] / sqrt( pow(y_face_area[i][2][k][0] ,2) + pow(y_face_area[i][2][k][1] ,2) + pow(y_face_area[i][2][k][2] ,2) ) ;  
			y_face_area_unit_vector[1] = y_face_area[i][2][k][1] / sqrt( pow(y_face_area[i][2][k][0] ,2) + pow(y_face_area[i][2][k][1] ,2) + pow(y_face_area[i][2][k][2] ,2) ) ;
			y_face_area_unit_vector[2] = y_face_area[i][2][k][2] / sqrt( pow(y_face_area[i][2][k][0] ,2) + pow(y_face_area[i][2][k][1] ,2) + pow(y_face_area[i][2][k][2] ,2) ) ;

			variablesvector[i][1][k][1] = (1-2*pow(y_face_area_unit_vector[0],2)) * variablesvector[i][2][k][1] -
			(2*y_face_area_unit_vector[0]*y_face_area_unit_vector[1])*variablesvector[i][2][k][2] -
			(2*y_face_area_unit_vector[0]*y_face_area_unit_vector[2])*variablesvector[i][2][k][3] ;

			variablesvector[i][1][k][2] = -(2*y_face_area_unit_vector[1]*y_face_area_unit_vector[0])*variablesvector[i][2][k][1] +
			(1-2*pow(y_face_area_unit_vector[1],2))*variablesvector[i][2][k][2] -
			(2*y_face_area_unit_vector[1]*y_face_area_unit_vector[2])*variablesvector[i][2][k][3] ;

			variablesvector[i][1][k][3] = -(2*y_face_area_unit_vector[2]*y_face_area_unit_vector[0])*variablesvector[i][2][k][1] -
			(2*y_face_area_unit_vector[2]*y_face_area_unit_vector[1])*variablesvector[i][2][k][2] +
			(1-2*pow(y_face_area_unit_vector[2],2))*variablesvector[i][2][k][3] ;

			// using zero order extrapolation
			variablesvector[i][1][k][0] = variablesvector[i][1][k][0] ;
			variablesvector[i][1][k][4] = variablesvector[i][2][k][4] ;

			// cout << "variablesvector[i][1][k][0] =  " << variablesvector[i][1][k][0] << endl ;
			// cout <<  "v    "<< variablesvector[i][1][k][4] << endl ;

			//using rightcell(j=3) filling j=0

			y_face_area_unit_vector[0] = y_face_area[i][3][k][0] / sqrt( pow(y_face_area[i][3][k][0] ,2) + pow(y_face_area[i][3][k][1] ,2) + pow(y_face_area[i][3][k][2] ,2) ) ;  
			y_face_area_unit_vector[1] = y_face_area[i][3][k][1] / sqrt( pow(y_face_area[i][3][k][0] ,2) + pow(y_face_area[i][3][k][1] ,2) + pow(y_face_area[i][3][k][2] ,2) ) ;
			y_face_area_unit_vector[2] = y_face_area[i][3][k][2] / sqrt( pow(y_face_area[i][3][k][0] ,2) + pow(y_face_area[i][3][k][1] ,2) + pow(y_face_area[i][3][k][2] ,2) ) ;

			variablesvector[i][0][k][1] = (1-2*pow(y_face_area_unit_vector[0],2)) * variablesvector[i][3][k][1] -
			(2*y_face_area_unit_vector[0]*y_face_area_unit_vector[1])*variablesvector[i][3][k][2] -
			(2*y_face_area_unit_vector[0]*y_face_area_unit_vector[2])*variablesvector[i][3][k][3] ;

			variablesvector[i][0][k][2] = -(2*y_face_area_unit_vector[1]*y_face_area_unit_vector[0])*variablesvector[i][3][k][1] +
			(1-2*pow(y_face_area_unit_vector[1],2))*variablesvector[i][3][k][2] -
			(2*y_face_area_unit_vector[1]*y_face_area_unit_vector[2])*variablesvector[i][3][k][3] ;

			variablesvector[i][0][k][3] = -(2*y_face_area_unit_vector[2]*y_face_area_unit_vector[0])*variablesvector[i][3][k][1] -
			(2*y_face_area_unit_vector[2]*y_face_area_unit_vector[1])*variablesvector[i][3][k][2] +
			(1-2*pow(y_face_area_unit_vector[2],2))*variablesvector[i][3][k][3] ;

			// using zero order extrapolation
			variablesvector[i][0][k][0] = variablesvector[i][3][k][0] ;
			variablesvector[i][0][k][4] = variablesvector[i][3][k][4] ;



			//using rightcell(j=Ny-3) filling j=Ny-2
			// here every term has been multiplied by -1 because cell area vector is -A
			y_face_area_unit_vector[0] = - y_face_area[i][Ny-3][k][0] / sqrt( pow(y_face_area[i][Ny-3][k][0] ,2) + pow(y_face_area[i][Ny-3][k][1] ,2) + pow(y_face_area[i][Ny-3][k][2] ,2) ) ;  
			y_face_area_unit_vector[1] = - y_face_area[i][Ny-3][k][1] / sqrt( pow(y_face_area[i][Ny-3][k][0] ,2) + pow(y_face_area[i][Ny-3][k][1] ,2) + pow(y_face_area[i][Ny-3][k][2] ,2) ) ;
			y_face_area_unit_vector[2] = - y_face_area[i][Ny-3][k][2] / sqrt( pow(y_face_area[i][Ny-3][k][0] ,2) + pow(y_face_area[i][Ny-3][k][1] ,2) + pow(y_face_area[i][Ny-3][k][2] ,2) ) ;
			
			variablesvector[i][Ny-2][k][1] = (1-2*pow(y_face_area_unit_vector[0],2)) * variablesvector[i][Ny-3][k][1] -
			(2*y_face_area_unit_vector[0]*y_face_area_unit_vector[1])*variablesvector[i][Ny-3][k][2] -
			(2*y_face_area_unit_vector[0]*y_face_area_unit_vector[2])*variablesvector[i][Ny-3][k][3] ;

			variablesvector[i][Ny-2][k][2] = -(2*y_face_area_unit_vector[1]*y_face_area_unit_vector[0])*variablesvector[i][Ny-3][k][1] +
			(1-2*pow(y_face_area_unit_vector[1],2))*variablesvector[i][Ny-3][k][2] -
			(2*y_face_area_unit_vector[1]*y_face_area_unit_vector[2])*variablesvector[i][Ny-3][k][3] ;

			variablesvector[i][Ny-2][k][3] = -(2*y_face_area_unit_vector[2]*y_face_area_unit_vector[0])*variablesvector[i][Ny-3][k][1] -
			(2*y_face_area_unit_vector[2]*y_face_area_unit_vector[1])*variablesvector[i][Ny-3][k][2] +
			(1-2*pow(y_face_area_unit_vector[2],2))*variablesvector[i][Ny-3][k][3] ;

			// using zero order extrapolation
			variablesvector[i][Ny-2][k][0] = variablesvector[i][Ny-3][k][0] ;
			variablesvector[i][Ny-2][k][4] = variablesvector[i][Ny-3][k][4] ;

			// using rightcell(j=Ny-4) filling j=Ny-1
			// here every term has been multiplied by -1 because cell area vector is -A
			y_face_area_unit_vector[0] = - y_face_area[i][Ny-4][k][0] / sqrt( pow(y_face_area[i][Ny-4][k][0] ,2) + pow(y_face_area[i][Ny-4][k][1] ,2) + pow(y_face_area[i][Ny-4][k][2] ,2) ) ;  
			y_face_area_unit_vector[1] = - y_face_area[i][Ny-4][k][1] / sqrt( pow(y_face_area[i][Ny-4][k][0] ,2) + pow(y_face_area[i][Ny-4][k][1] ,2) + pow(y_face_area[i][Ny-4][k][2] ,2) ) ;
			y_face_area_unit_vector[2] = - y_face_area[i][Ny-4][k][2] / sqrt( pow(y_face_area[i][Ny-4][k][0] ,2) + pow(y_face_area[i][Ny-4][k][1] ,2) + pow(y_face_area[i][Ny-4][k][2] ,2) ) ;
			
			variablesvector[i][Ny-1][k][1] = (1-2*pow(y_face_area_unit_vector[0],2)) * variablesvector[i][Ny-4][k][1] -
			(2*y_face_area_unit_vector[0]*y_face_area_unit_vector[1])*variablesvector[i][Ny-4][k][2] -
			(2*y_face_area_unit_vector[0]*y_face_area_unit_vector[2])*variablesvector[i][Ny-4][k][3] ;

			variablesvector[i][Ny-1][k][2] = -(2*y_face_area_unit_vector[1]*y_face_area_unit_vector[0])*variablesvector[i][Ny-4][k][1] +
			(1-2*pow(y_face_area_unit_vector[1],2))*variablesvector[i][Ny-4][k][2] -
			(2*y_face_area_unit_vector[1]*y_face_area_unit_vector[2])*variablesvector[i][Ny-4][k][3] ;

			variablesvector[i][Ny-1][k][3] = -(2*y_face_area_unit_vector[2]*y_face_area_unit_vector[0])*variablesvector[i][Ny-4][k][1] -
			(2*y_face_area_unit_vector[2]*y_face_area_unit_vector[1])*variablesvector[i][Ny-4][k][2] +
			(1-2*pow(y_face_area_unit_vector[2],2))*variablesvector[i][Ny-4][k][3] ;

			// using zero order extrapolation
			variablesvector[i][Ny-1][k][0] = variablesvector[i][Ny-4][k][0] ;
			variablesvector[i][Ny-1][k][4] = variablesvector[i][Ny-4][k][4] ;

		}
	}

	// // updating the z ghost cells
	for (int i = 2; i < Nx-2; ++i)
	{
		for (int j = 2; j < Ny-2; ++j)
		{
			//using cell (k=2) filling k=1
			double z_face_area_unit_vector[3] ={0.0,0.0,1.0};
			z_face_area_unit_vector[0] = z_face_area[i][j][2][0] / sqrt( pow(z_face_area[i][j][2][0] ,2) + pow(z_face_area[i][j][2][1] ,2) + pow(z_face_area[i][j][2][2] ,2) ) ;  
			z_face_area_unit_vector[1] = z_face_area[i][j][2][1] / sqrt( pow(z_face_area[i][j][2][0] ,2) + pow(z_face_area[i][j][2][1] ,2) + pow(z_face_area[i][j][2][2] ,2) ) ;
			z_face_area_unit_vector[2] = z_face_area[i][j][2][2] / sqrt( pow(z_face_area[i][j][2][0] ,2) + pow(z_face_area[i][j][2][1] ,2) + pow(z_face_area[i][j][2][2] ,2) ) ;

			variablesvector[i][j][1][1] = (1-2*pow(z_face_area_unit_vector[0],2)) * variablesvector[i][j][2][1] -
			(2*z_face_area_unit_vector[0]*z_face_area_unit_vector[1])*variablesvector[i][j][2][2] -
			(2*z_face_area_unit_vector[0]*z_face_area_unit_vector[2])*variablesvector[i][j][2][3] ;

			variablesvector[i][j][1][2] = -(2*z_face_area_unit_vector[1]*z_face_area_unit_vector[0])*variablesvector[i][j][2][1] +
			(1-2*pow(z_face_area_unit_vector[1],2))*variablesvector[i][j][2][2] -
			(2*z_face_area_unit_vector[1]*z_face_area_unit_vector[2])*variablesvector[i][j][2][3] ;

			variablesvector[i][j][1][3] = -(2*z_face_area_unit_vector[2]*z_face_area_unit_vector[0])*variablesvector[i][j][2][1] -
			(2*z_face_area_unit_vector[2]*z_face_area_unit_vector[1])*variablesvector[i][j][2][2] +
			(1-2*pow(z_face_area_unit_vector[2],2))*variablesvector[i][j][2][3] ;

			// using zero order extrapolation
			variablesvector[i][j][1][0] = variablesvector[i][j][2][0] ;
			variablesvector[i][j][1][4] = variablesvector[i][j][2][4] ;

			//using cell (k=3) filling k=0

			z_face_area_unit_vector[0] = z_face_area[i][j][3][0] / sqrt( pow(z_face_area[i][j][3][0] ,2) + pow(z_face_area[i][j][3][1] ,2) + pow(z_face_area[i][j][3][2] ,2) ) ;  
			z_face_area_unit_vector[1] = z_face_area[i][j][3][1] / sqrt( pow(z_face_area[i][j][3][0] ,2) + pow(z_face_area[i][j][3][1] ,2) + pow(z_face_area[i][j][3][2] ,2) ) ;
			z_face_area_unit_vector[2] = z_face_area[i][j][3][2] / sqrt( pow(z_face_area[i][j][3][0] ,2) + pow(z_face_area[i][j][3][1] ,2) + pow(z_face_area[i][j][3][2] ,2) ) ;

			variablesvector[i][j][0][1] = (1-2*pow(z_face_area_unit_vector [0],2)) * variablesvector[i][j][3][1] -
			(2*z_face_area_unit_vector[0]*z_face_area_unit_vector[1])*variablesvector[i][j][3][2] -
			(2*z_face_area_unit_vector[0]*z_face_area_unit_vector[2])*variablesvector[i][j][3][3] ;

			variablesvector[i][j][0][2] = -(2*z_face_area_unit_vector[1]*z_face_area_unit_vector[0])*variablesvector[i][j][3][1] +
			(1-2*pow(z_face_area_unit_vector [1],2))*variablesvector[i][j][3][2] -
			(2*z_face_area_unit_vector[1]*z_face_area_unit_vector[2])*variablesvector[i][j][3][3] ;

			variablesvector[i][j][0][3] = -(2*z_face_area_unit_vector [2]*z_face_area_unit_vector [0])*variablesvector[i][j][3][1] -
			(2*z_face_area_unit_vector [2]*z_face_area_unit_vector [1])*variablesvector[i][j][3][2] +
			(1-2*pow(z_face_area_unit_vector [2],2))*variablesvector[i][j][3][3] ;

			// using zero order extrapolation
			variablesvector[i][j][0][0] = variablesvector[i][j][3][0] ;
			variablesvector[i][j][0][4] = variablesvector[i][j][3][4] ;

			// using cell (k=Nz-3) filling k=Nz-2

			z_face_area_unit_vector[0] = - z_face_area[i][j][Nz-3][0] / sqrt( pow(z_face_area[i][j][Nz-3][0] ,2) + pow(z_face_area[i][j][Nz-3][1] ,2) + pow(z_face_area[i][j][Nz-3][2] ,2) ) ;  
			z_face_area_unit_vector[1] = - z_face_area[i][j][Nz-3][1] / sqrt( pow(z_face_area[i][j][Nz-3][0] ,2) + pow(z_face_area[i][j][Nz-3][1] ,2) + pow(z_face_area[i][j][Nz-3][2] ,2) ) ;
			z_face_area_unit_vector[2] = - z_face_area[i][j][Nz-3][2] / sqrt( pow(z_face_area[i][j][Nz-3][0] ,2) + pow(z_face_area[i][j][Nz-3][1] ,2) + pow(z_face_area[i][j][Nz-3][2] ,2) ) ;
			
			variablesvector[i][j][Nz-2][1] = (1-2*pow(z_face_area_unit_vector[0],2)) * variablesvector[i][j][Nz-3][1] -
			(2*z_face_area_unit_vector[0]*z_face_area_unit_vector[1])*variablesvector[i][j][Nz-3][2] -
			(2*z_face_area_unit_vector[0]*z_face_area_unit_vector[2])*variablesvector[i][j][Nz-3][3] ;

			variablesvector[i][j][Nz-2][2] = -(2*z_face_area_unit_vector[1]*z_face_area_unit_vector[0])*variablesvector[i][j][Nz-3][1] +
			(1-2*pow(z_face_area_unit_vector[1],2))*variablesvector[i][j][Nz-3][2] -
			(2*z_face_area_unit_vector[1]*z_face_area_unit_vector[2])*variablesvector[i][j][Nz-3][3] ;

			variablesvector[i][j][Nz-2][3] = -(2*z_face_area_unit_vector[2]*z_face_area_unit_vector[0])*variablesvector[i][j][Nz-3][1] -
			(2*z_face_area_unit_vector[2]*z_face_area_unit_vector[1])*variablesvector[i][j][Nz-3][2] +
			(1-2*pow(z_face_area_unit_vector[2],2))*variablesvector[i][j][Nz-3][3] ;
			// using zero order extrapolation
			variablesvector[i][j][Nz-2][0] = variablesvector[i][j][Nz-3][0] ;
			variablesvector[i][j][Nz-2][4] = variablesvector[i][j][Nz-3][4] ;

			//using cell(k=Nz-4) filling k=Nz-1
			z_face_area_unit_vector[0] = - z_face_area[i][j][Nz-4][0] / sqrt( pow(z_face_area[i][j][Nz-4][0] ,2) + pow(z_face_area[i][j][Nz-4][1] ,2) + pow(z_face_area[i][j][Nz-4][2] ,2) ) ;  
			z_face_area_unit_vector[1] = - z_face_area[i][j][Nz-4][1] / sqrt( pow(z_face_area[i][j][Nz-4][0] ,2) + pow(z_face_area[i][j][Nz-4][1] ,2) + pow(z_face_area[i][j][Nz-4][2] ,2) ) ;
			z_face_area_unit_vector[2] = - z_face_area[i][j][Nz-4][2] / sqrt( pow(z_face_area[i][j][Nz-4][0] ,2) + pow(z_face_area[i][j][Nz-4][1] ,2) + pow(z_face_area[i][j][Nz-4][2] ,2) ) ;
			

			variablesvector[i][j][Nz-1][1] = (1-2*pow(z_face_area_unit_vector[0],2)) * variablesvector[i][j][Nz-4][1] -
			(2*z_face_area_unit_vector[0]*z_face_area_unit_vector[1])*variablesvector[i][j][Nz-4][2] -
			(2*z_face_area_unit_vector[0]*z_face_area_unit_vector[2])*variablesvector[i][j][Nz-4][3] ;

			variablesvector[i][j][Nz-1][2] = -(2*z_face_area_unit_vector[1]*z_face_area_unit_vector[0])*variablesvector[i][j][Nz-4][1] +
			(1-2*pow(z_face_area_unit_vector[1],2))*variablesvector[i][j][Nz-4][2] -
			(2*z_face_area_unit_vector[1]*z_face_area_unit_vector[2])*variablesvector[i][j][Nz-4][3] ;

			variablesvector[i][j][Nz-1][3] = -(2*z_face_area_unit_vector[2]*z_face_area_unit_vector[0])*variablesvector[i][j][Nz-4][1] -
			(2*z_face_area_unit_vector[2]*z_face_area_unit_vector[1])*variablesvector[i][j][Nz-4][2] +
			(1-2*pow(z_face_area_unit_vector[2],2))*variablesvector[i][j][Nz-4][3] ; 
			// using zero order extrapolation
			variablesvector[i][j][Nz-1][0] = variablesvector[i][j][Nz-4][0] ;
			variablesvector[i][j][Nz-1][4] = variablesvector[i][j][Nz-4][4] ;
			// cout<<" x vel " << variablesvector[i][j][Nz-1][1] <<" y vel " << variablesvector[i][j][Nz-1][2] <<
			// "    z vel   "<< variablesvector[i][j][Nz-1][3] << endl ; 
		}
	}
}

