#include <iostream>
#include "vector"
#include "math.h"
using namespace std;
void getNormal(vector<double> & UnitNormal, vector<double> areaVector)
{
	double vectorMagnitude = sqrt(areaVector[0]*areaVector[0] + areaVector[1]*
		areaVector[1] + areaVector[2]*areaVector[2]);
	UnitNormal[0] = areaVector[0]/vectorMagnitude;
	UnitNormal[1] = areaVector[1]/vectorMagnitude;
	UnitNormal[2] = areaVector[2]/vectorMagnitude;
	areaVector[0] = 0;
	cout << areaVector[0] << " " << areaVector[1] << " " << areaVector[2] << endl;	
}

/** \brief This function implements the wall boundary condition
*\param AreaVectors Surface faces area vectors.
*\param LiveCellConservedVariables Conserved variables array for the live cell.
*\param GhostCellConservedVariables Conserved variables array for the ghost cell.
*\return void
*/
void wallBC(vector<double> AreaVectors)
{
	std::vector<double> n(3); // Unit normal vector to the face
	getNormal(n,AreaVectors);
}
int main()
{
	std::vector<double> u(3);
	std::vector<std::vector<double> > v(3,u) ;

	v[0][0] = 1;
	v[0][1] = 1;
	v[0][2] = 1;

	getNormal(u,v[0]);
	cout << u[0] << " " << u[1] << " "  << u[2] << endl;
	cout << v[0][0] << " " << v[0][1] << " "  << v[0][2] << endl;
	return 0;
}