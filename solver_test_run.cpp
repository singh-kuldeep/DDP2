#include "iostream"
#include <vector>
#include <fstream>
#include "math.h"
#include "time.h"
#include "BC.h"
#include "string"
#include "array_tester.h"
#include "colortext.h"
// #include "initial_condition.h"
// #include "netfluxAUSM.h" // AUSM
// #include "netfluxRoe.h" // Roe

// #include "grid.h" // Headers for grids 
// #include "ghostcell.h" // Headers for ghost cells
// string fail()
// {
// 	return "  \033[1;31m failed the test case \033[0m\n ";
// }

// string pass()
// {
// 	return "  \033[1;32m passed the test case \033[0m\n ";
// }

// string blue(string inputstring)
// {
// 	string left = "  \033[1;36m ";
// 	left.append(inputstring);
// 	left.append(" \033[0m");
// 	return left;
// }


using namespace std ;
int main()
{	
	double e = 0.00001;
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

	// SubSonicInletBC() test
	SubSonicInletBC(Ughost,{1.00,0.0, 1000.0, 1000.0, 5e5},2.25,2000.0,1.0,1.0);
	if(2.25-e<Ughost[0]<2.25+e && 2000.0*2.25-e<Ughost[1]<2000.0*2.25+e && 
		1.0*2.25-e<Ughost[2]<1.0*2.25-e && 1.0*2.25-e<Ughost[3]<1.0*2.25+e &&
		388889.333333-e<Ughost[4]<e+388889.333333 )
	{
		cout << blue("SubSonicInletBC()")<<pass() << endl; 
	}
	else
	{
		cout << blue("SubSonicInletBC()") << fail() << endl; 
	}
	
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
	
	


	return 0;
}