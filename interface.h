// checked
#include "math.h"
#include "iostream"

#define gamma 1.4
#define gasconstant 287.14
#define heatcapacityconstantvolume 717.5

using namespace std ;

class interface
{
	public:
	double densityinterface ;
	double xvelocityinterface ;
	double yvelocityinterface ;
	double zvelocityinterface ;
	double enthalpyinterface ;

	double vectorjumpinterface[5] ;
	double eigenvalue[5] ;


	double eigenvectormatrix[5][5] ;
	double eigenvectormatrixinvers[5][5] ; 

	double alphavectorinterface[5] ; // eigenvectormatrixinvers[5][5]*vectorjumpinterface

	double muvectorinterface[5] ; // this is some paramere = delta t * eigenvalue
	double Zvectorinterface[5] ; // this is dame as muvectorinterface

	double pshivectorinterface[5] ;
	double gvectorinterface[5];


	interface(vector<double>& vectorleft,vector<double>& vectorright,vector<double>& areavectorinterface,
				double volumeleft, double volumeright, double deltat) 
	{
		double volumeinterface = 0.5*(volumeleft + volumeright) ;

		double pressureleft = (gamma -1)*( vectorleft[4] - 0.5*(pow(vectorleft[1],2)+pow(vectorleft[2],2)+
			pow(vectorleft[3],2))/vectorleft[0] ) ;
		double pressureright = (gamma -1)*( vectorright[4] - 0.5*(pow(vectorright[1],2)+pow(vectorright[2],2)+
			pow(vectorright[3],2))/vectorright[0] ) ;
		
		double enthalpyleft = (vectorleft[4] + pressureleft)/vectorleft[0];
		double enthalpyright = (vectorright[4] + pressureright)/vectorright[0];

		densityinterface = sqrt(vectorleft[0]*vectorright[0]) ; 
		xvelocityinterface = (vectorleft[1] + vectorright[1]*sqrt(vectorleft[0]/vectorright[0]))/
		(vectorleft[0] +sqrt(vectorright[0]*vectorleft[0])) ;
		yvelocityinterface = (vectorleft[2] + vectorright[2]*sqrt(vectorleft[0]/vectorright[0]))/
		(vectorleft[0] +sqrt(vectorright[0]*vectorleft[0])) ;
		zvelocityinterface = (vectorleft[3] + vectorright[3]*sqrt(vectorleft[0]/vectorright[0]))/
		(vectorleft[0] +sqrt(vectorright[0]*vectorleft[0])) ;
		enthalpyinterface = (enthalpyleft + enthalpyright*sqrt(vectorright[0]/vectorleft[0])) /
		(1+sqrt(vectorright[0]/vectorleft[0])) ; 


		// engenvalues 
		double hx = areavectorinterface[0] /volumeinterface ;
		double hy = areavectorinterface[1] /volumeinterface ;
		double hz = areavectorinterface[2] /volumeinterface ;
		double hn = sqrt(pow(areavectorinterface[0],2) + pow(areavectorinterface[1],2) + 
			pow(areavectorinterface[2],2)) / volumeinterface ;
		double Ucont = xvelocityinterface*hx + yvelocityinterface*hy + zvelocityinterface*hz ;
		double soundvelocityinterface = sqrt((gamma -1)*(enthalpyinterface - 0.5*(pow(xvelocityinterface,2)+
			pow(yvelocityinterface,2) + pow(zvelocityinterface,2)))) ;

		eigenvalue[0] = Ucont - soundvelocityinterface*hn ;
		eigenvalue[1] = Ucont ;
		eigenvalue[2] = Ucont ;
		eigenvalue[3] = Ucont ;
		eigenvalue[4] = Ucont + soundvelocityinterface*hn ;

		// Jump vector at interface 
		vectorjumpinterface[0] = volumeinterface*(vectorright[0] - vectorleft[0]) ;
		vectorjumpinterface[1] = volumeinterface*(vectorright[1] - vectorleft[1]) ;
		vectorjumpinterface[2] = volumeinterface*(vectorright[2] - vectorleft[2]) ;
		vectorjumpinterface[3] = volumeinterface*(vectorright[3] - vectorleft[3]) ;
		vectorjumpinterface[4] = volumeinterface*(vectorright[4] - vectorleft[4]) ;

		// Defying the eigenvector matrix "R"
		double hdesx = hx/hn ; 
		double hdesy = hy/hn ; 
		double hdesz = hz/hn ; 
		double Phi = xvelocityinterface*hdesx + yvelocityinterface*hdesy + zvelocityinterface*hdesz ;
		double q = sqrt(pow(xvelocityinterface,2)+pow(yvelocityinterface,2)+pow(zvelocityinterface,2)) ;

		eigenvectormatrix[0][0] = 1 ; 
		eigenvectormatrix[0][1] = 1 ;
		eigenvectormatrix[0][2] = 0 ;
		eigenvectormatrix[0][3] = 0 ; 
		eigenvectormatrix[0][4] = 1 ;

		eigenvectormatrix[1][0] = xvelocityinterface - hdesx * soundvelocityinterface ; 
		eigenvectormatrix[1][1] = xvelocityinterface ;
		eigenvectormatrix[1][2] = hdesy ;
		eigenvectormatrix[1][3] = hdesz ; 
		eigenvectormatrix[1][4] = xvelocityinterface + hdesx*soundvelocityinterface ;

		eigenvectormatrix[2][0] = yvelocityinterface - hdesy * soundvelocityinterface ; 
		eigenvectormatrix[2][1] = yvelocityinterface ;
		eigenvectormatrix[2][2] = hdesz ;
		eigenvectormatrix[2][3] = hdesx ; 
		eigenvectormatrix[2][4] = yvelocityinterface + hdesy*soundvelocityinterface ;

		eigenvectormatrix[3][0] = zvelocityinterface - hdesz * soundvelocityinterface ; 
		eigenvectormatrix[3][1] = zvelocityinterface ;
		eigenvectormatrix[3][2] = hdesx ;
		eigenvectormatrix[3][3] = hdesy ; 
		eigenvectormatrix[3][4] = zvelocityinterface + hdesz*soundvelocityinterface ;

		eigenvectormatrix[4][0] = enthalpyinterface - (hdesx*xvelocityinterface + hdesy*yvelocityinterface + 
			hdesz*zvelocityinterface) * soundvelocityinterface ; 
		eigenvectormatrix[4][1] = 0.5 * pow(q,2);
		eigenvectormatrix[4][2] = hdesx*zvelocityinterface + hdesz*yvelocityinterface + hdesy*xvelocityinterface;
		eigenvectormatrix[4][3] = hdesy*zvelocityinterface + hdesx*yvelocityinterface + hdesz*xvelocityinterface;
		eigenvectormatrix[4][4] = enthalpyinterface + (hdesx*xvelocityinterface + hdesy*yvelocityinterface + 
			hdesz*zvelocityinterface) * soundvelocityinterface ;

		// eigenvectormatrixinvers defying
		eigenvectormatrixinvers[0][0] = 0.5*(0.5*(pow(q,2))*((gamma-1)/pow(soundvelocityinterface,2))+ 
			(Phi/soundvelocityinterface));
		eigenvectormatrixinvers[0][1] = -0.5*(xvelocityinterface*((gamma-1)/pow(soundvelocityinterface,2))+ 
			(hdesx/soundvelocityinterface));
		eigenvectormatrixinvers[0][2] = -0.5*(yvelocityinterface*((gamma-1)/pow(soundvelocityinterface,2))+ 
			(hdesy/soundvelocityinterface));
		eigenvectormatrixinvers[0][3] = -0.5*(zvelocityinterface*((gamma-1)/pow(soundvelocityinterface,2))+ 
			(hdesz/soundvelocityinterface));
		eigenvectormatrixinvers[0][4] = 0.5*(gamma -1)/ pow(soundvelocityinterface,2) ;

		eigenvectormatrixinvers[1][0] = 1 - 0.5*pow(q,2)*((gamma-1)/pow(soundvelocityinterface,2)) ; 
		eigenvectormatrixinvers[1][1] = xvelocityinterface*(gamma-1)/pow(soundvelocityinterface,2) ; 
		eigenvectormatrixinvers[1][2] = yvelocityinterface*(gamma-1)/pow(soundvelocityinterface,2) ;   
		eigenvectormatrixinvers[1][3] = zvelocityinterface*(gamma-1)/pow(soundvelocityinterface,2) ; 
		eigenvectormatrixinvers[1][4] = -(gamma-1)/pow(soundvelocityinterface,2) ;

		eigenvectormatrixinvers[2][0] = -(hdesy*xvelocityinterface + hdesz*yvelocityinterface + 
			hdesx*zvelocityinterface);
		eigenvectormatrixinvers[2][1] = hdesy ; 
		eigenvectormatrixinvers[2][2] = hdesz ;
		eigenvectormatrixinvers[2][3] = hdesx ; 
		eigenvectormatrixinvers[2][4] = 0 ;

		eigenvectormatrixinvers[3][0] = -(hdesz*xvelocityinterface + hdesx*yvelocityinterface + 
			hdesy*zvelocityinterface);
		eigenvectormatrixinvers[3][1] = hdesz ; 
		eigenvectormatrixinvers[3][2] = hdesx ;
		eigenvectormatrixinvers[3][3] = hdesy ; 
		eigenvectormatrixinvers[3][4] = 0 ;

		eigenvectormatrixinvers[4][0] = 0.5*(0.5*(pow(q,2))*((gamma-1)/pow(soundvelocityinterface,2)) - 
			(Phi/soundvelocityinterface));
		eigenvectormatrixinvers[4][1] = 0.5*( -xvelocityinterface*((gamma-1)/pow(soundvelocityinterface,2))+ 
			(hdesx/soundvelocityinterface));
		eigenvectormatrixinvers[4][2] = 0.5*( -yvelocityinterface*((gamma-1)/pow(soundvelocityinterface,2))+ 
			(hdesy/soundvelocityinterface));
		eigenvectormatrixinvers[4][3] = 0.5*( -zvelocityinterface*((gamma-1)/pow(soundvelocityinterface,2))+ 
			(hdesz/soundvelocityinterface));
		eigenvectormatrixinvers[4][4] = 0.5*(gamma -1)/ pow(soundvelocityinterface,2) ;

//alphavectorinterface defying 
 		for (int i = 0; i < 5; ++i)
		{
			alphavectorinterface[i] = 0.0 ;
			for (int l = 0; l < 5; ++l)
			{
				alphavectorinterface[i] += eigenvectormatrixinvers[i][l]*vectorjumpinterface[l] ;
			}
				// cout <<   alphavectorinterface[i] << endl ;
		}

// Defying the muvectorinterface or Zvectorinterface
	for (int i = 0; i < 5; ++i)
		{
			muvectorinterface[i] = deltat * eigenvalue[i] ;
			Zvectorinterface[i] = deltat * eigenvalue[i] ;
		}
			
// // pshivectorinterface defying (for invisid flow)
		for (int i = 0; i < 5; ++i)
		{
			pshivectorinterface[i] = pow(Zvectorinterface[i],2) + 0.25 ;
		}
// // pshivectorinterface defying (for viscus flow)
// 	double deltaf = 0.2 ; // this is given constant value (it is between 0.1 to 0.5)
// 		for (int i = 0; i < 5; ++i)
// 		{
// 			if(fabs(Zvectorinterface[i]) >= deltaf ){
// 				pshivectorinterface[i] = fabs(Zvectorinterface[i]);
// 			}
// 			else
// 			{
// 				pshivectorinterface[i] = 0.5*(pow(Zvectorinterface[i],2)+pow(deltaf,2)) / deltaf ;
// 			}
// 			// cout << "pshivectorinterface [" << i+1 <<"]     " << pshivectorinterface[i] << endl ;
// 		}

// gvectorinterface defying
		for (int i = 0; i < 5; ++i)
		{
			gvectorinterface[i] = 0.5*(pshivectorinterface[i] - 
				pow(Zvectorinterface[i],2))*alphavectorinterface[i];
		}
	 };
	// ~interface();
};