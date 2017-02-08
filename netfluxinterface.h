#include "math.h"
#include "iostream"
#include "eulerflux.h"
// #include "viscusflux.h"
#include "diffusionfluxinterface.h"
#define gamma 1.4
#define gasconstant 287.14
#define heatcapacityconstantvolume 717.5

using namespace std ;

class netfluxinterface
{
	public:
	double netflux[5] ;
	
	netfluxinterface(vector<double>& vectorleftminus, vector<double>& vectorleft, vector<double>& vectorright,
vector<double>& vectorrightplus, vector<double>& areavectorleft, vector<double>& areavectorright, vector<double>& areavectorrightplus,
double volumeleftmins, double volumeleft, double volumeright, double volumerightplus, double deltat){

eulerflux left(vectorleft);
eulerflux right(vectorright);

double volumeinterfaceleft = (volumeleftmins +volumeleft)/2 ;
double volumeinterfaceright = (volumeleft + volumeright)/2 ;
double volumeinterfacerightplus = (volumeright + volumerightplus)/2 ;

diffusionfluxinterface diffusion(vectorleftminus,vectorleft,vectorright,vectorrightplus,
	areavectorleft,areavectorright,areavectorrightplus,volumeleftmins,volumeleft,
	volumeright,volumerightplus,deltat); 

// Averaged interface Euler flux 
		double hx = areavectorright[0] / volumeinterfaceright ;
		double hy = areavectorright[1] / volumeinterfaceright ;	
		double hz = areavectorright[2] / volumeinterfaceright ;

		double xnetfluxinterface[5] ; 
		double ynetfluxinterface[5] ; 
		double znetfluxinterface[5] ; 
		for (int i = 0; i < 5; ++i)
		{
			xnetfluxinterface[i] = 0.5*(left.xeulerflux[i]+right.xeulerflux[i]) ; 
			ynetfluxinterface[i] = 0.5*(left.yeulerflux[i]+right.yeulerflux[i]) ; 
			znetfluxinterface[i] = 0.5*(left.zeulerflux[i]+right.zeulerflux[i]) ;

			netflux[i] = (xnetfluxinterface[i]*hx + ynetfluxinterface[i]*hy + znetfluxinterface[i]*hz)*volumeinterfaceright + 
			0.5*diffusion.diffusionfluxvector[i];
		}
	};
	// ~netfluxinterface();
	
};