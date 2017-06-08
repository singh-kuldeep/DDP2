#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "iostream"
#include <vector>
#include <fstream>
#include "math.h"
#include "time.h"
#include "string"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#include "BC.h"
#include "array_tester.h"
#include "colortext.h"
#include "netfluxAUSM.h"


using namespace std ;
int main()
{	
	double e = 1;
	#if 0
	// getnormal() test cases 
	vector<double> UnitNormal(3);
	getNormal(UnitNormal,{0,1,0});
	if(UnitNormal[1] == 1)
	{
		cout << blue("getNormal()") <<  pass() << endl; 
	}
	else
	{
		cout << blue("getNormal()") << fail() << endl; 
	}
	
	getNormal(UnitNormal,{1,1,1});

	if(UnitNormal[0] == 1/sqrt(3) &&
	 UnitNormal[1] == 1/sqrt(3) &&
	 UnitNormal[2] == 1/sqrt(3))
	{
		cout << blue("getNormal()") << pass() << endl; 
	}
	else
	{
		cout << blue("getNormal()") << fail() << endl; 
	}

	// WallBC() test cases 
	vector<double> Ughost(5);
	WallBC(Ughost, {1.25, 1000.0, 0.0, 0.0, 5e5}, {0.5,0.0,0.0});
	if(Ughost[0] == 1.25 && Ughost[1] == -1000.0 && Ughost[2] == 0.0 && Ughost[3] == 0 && Ughost[4] == 5e5 )
	{
		cout << blue("WallBC()")<<pass() << endl; 
	}
	else
	{
		cout << blue("WallBC()") << fail() << endl; 
	}

	WallBC(Ughost, {1.25, -1000.0, -1000.0, -1000.0, 5e5}, {0.5,0.5,0.5});
	if(
		1.25-e<Ughost[0]<1.25+e
		&& 1000-e<Ughost[1]<1000+e 
		&& 1000-e<Ughost[2]<1000+e 
		&& 1000-e<Ughost[3]<1000+e 
		&& 5e5-e<Ughost[4]<5e5+e 
	)
	{
		cout << blue("WallBC()") << pass() << endl; 
	}
	else
	{
		cout << blue("WallBC()") << fail() << endl; 
	}
	#endif
	#if 
	// SubSonicInletBC() test
	SubSonicInletBC(Ughost,{1.00,0.0, 1000.0, 1000.0, 5e5},2.25,2000.0,1.0,1.0);
	if(2.25-e<Ughost[0]<2.25+e && 2000.0*2.25-e<Ughost[1]<2000.0*2.25+e && 
		1.0*2.25-e<Ughost[2]<1.0*2.25-e && 1.0*2.25-e<Ughost[3]<1.0*2.25+e &&
		388889.333333-e<Ughost[4]<e+388889.333333 )
	{
		cout << blue("SubSonicInletBC()")<< pass() << endl; 
	}
	else
	{
		cout << blue("SubSonicInletBC()") << fail() << endl; 
	}
	#endif
	#if 0
	// SuperSonicInletBC() test
	SuperSonicInletBC(Ughost,2.25,2000.0,1.0,1.0,-200000);
	if(2.25-e<Ughost[0]<2.25+e && 2000.0*2.25-e<Ughost[1]<2000.0*2.25+e && 
		1.0*2.25-e<Ughost[2]<1.0*2.25-e && 1.0*2.25-e<Ughost[3]<1.0*2.25+e &&
		388889.333333-e<Ughost[4]<e+388889.333333 )
	{	
		cout << blue("SuperSonicInletBC()")<<pass() << endl; 
	}
	else
	{
		cout << blue("SuperSonicInletBC()") << fail() << endl; 
	}

	// SuperSonicExitBC() test
	SuperSonicExitBC(Ughost,{2.25,2000.0,1.0,1.0,5e5});
	if(
		2.25-e<Ughost[0]<2.25+e 
		&& 2000.0-e<Ughost[1]<2000.0+e 
		&& 1.0-e<Ughost[2]<1.0+e 
		&& 1.0-e<Ughost[3]<1.0+e 
		&& 5e5-e<Ughost[4]<e+5e5 
		)
	{
		cout << blue("SuperSonicExitBC()")<< pass() << endl; 
	}
	else
	{
		cout << blue("SuperSonicExitBC()") << fail() << endl; 
	}
		
	// SubSonicExitBC() test
	SubSonicExitBC(Ughost,{1.00,0.0, 1000.0, 1000.0, 5e5},1e5);
	if(2.25-e<Ughost[0]<2.25+e && 2000.0*2.25-e<Ughost[1]<2000.0*2.25+e && 
		1.0*2.25-e<Ughost[2]<1.0*2.25-e && 1.0*2.25-e<Ughost[3]<1.0*2.25+e &&
		1250000-e<Ughost[4]<e+1250000 )
	{
		cout << blue("SubSonicExitBC()")<<pass() << endl; 
	}
	else
	{
		cout << blue("SubSonicExitBC()") << fail() << endl; 
	}

	SubSonicExitBC(Ughost,{1.00,0.0, 1000.0, 1000.0, 5e5},1e5);
	if(2.25-e<Ughost[0]<2.25+e && 2000.0*2.25-e<Ughost[1]<2000.0*2.25+e && 
		1.0*2.25-e<Ughost[2]<1.0*2.25-e && 1.0*2.25-e<Ughost[3]<1.0*2.25+e &&
		1250000-e<Ughost[4]<e+1250000 )
	{
		cout << blue("SubSonicExitBC()")<<pass() << endl; 
	}
	else
	{
		cout << blue("SubSonicExitBC()") << fail() << endl; 
	}

	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> Dim4;
	
	vector<vector<vector<vector<double> > > > Coordinate;	
	Coordinate.resize(1,Dim3(1,Dim2(1,Dim1(1))));
	Coordinate[0][0][0][0] = sqrt(-1);	
 	cout << test4DArray("Coordinate",Coordinate,1,1,1,1) << endl;
	#endif 
	
 	// Flux vectro testing 
 	
 // 	netfluxAUSM face({1.23,829,0.0,0.0,529379.74},{1.23,829,0.0,0.0,529379.74},{1.0,0.,0.0});

 // 	double rho = 1.23;
 // 	double u = 829/1.23;
 // 	double p = 0.4*(529379.74-0.5*rho*u*u);
 	
 // 	std::vector<double> flux(5);
 // 	flux[0] = rho*u;
 // 	flux[1] = rho*u*u + p ;
 // 	flux[2] = 0;
 // 	flux[3] = 0;
 // 	flux[4] = u*(529379.74+p);
 // 	for (int i = 0; i < 5; ++i)
 // 	{
 // 		cout << "NetFlux[" << i << "] ->" << face.NetFlux[i] << endl; 
 // 		cout << "Flux[" << i << "] ->" << flux[i] << endl; 
 // 	}
 // 	if (-e+flux[0]<face.NetFlux[0]<e+flux[0] 
 // 		&& -e+flux[1]<face.NetFlux[1]<e+flux[1] 
 // 		&& -4+flux[2]<face.NetFlux[2]<4+flux[2]
 // 		 && -4+flux[3]<face.NetFlux[3]<4+flux[3] 
 // 		&& -12+flux[4]<face.NetFlux[4]<12+flux[4]
 // 		)
 // 	{
	// 	cout << blue("Class :: netfluxAUSM ") << blue("working fine") << endl; 
	// }
	// else
	// {
	// 	cout << blue("class :: netfluxAUSM ") << red("NOT working fine") << endl; 
	// }

	return 0;
}