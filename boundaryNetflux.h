#ifndef BOUNDARYNETFLUX_H
#define BOUNDARYNETFLUX_H
#include <math.h>
double Magnitude(vector<double> v)
{
	return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

void boundaryNetflux(
	vector<vector<vector<vector<double> > > >	& iFacesFlux,
	vector<vector<vector<vector<double> > > >	& jFacesFlux,
	vector<vector<vector<vector<double> > > >	& kFacesFlux,
	vector<vector<vector<vector<double> > > >	iFaceAreaVector,
	vector<vector<vector<vector<double> > > >	jFaceAreaVector,
	vector<vector<vector<vector<double> > > >	kFaceAreaVector,
	int Ni, int Nj, int Nk)
{
	double NetMassFlux=0.0;
	double NetXMomentumFlux=0.0;
	double NetYMomentumFlux = 0.0;
	double NetZMomentumFlux = 0.0;
	double NetEnergyFlux = 0.0;

	// i faces
	for (int j = 0; j < Nj; ++j)
	{
		for (int k = 0; k < Nk; ++k)
		{
			// i = 0 and i = Ni
			NetMassFlux	+= iFacesFlux[0][j][k][0]*Magnitude(iFaceAreaVector[0][j][k]);
			NetMassFlux	-= iFacesFlux[Ni][j][k][0]*Magnitude(iFaceAreaVector[Ni][j][k]);

			NetXMomentumFlux += iFacesFlux[0][j][k][1]*Magnitude(iFaceAreaVector[0][j][k]);
			NetXMomentumFlux -= iFacesFlux[Ni][j][k][1]*Magnitude(iFaceAreaVector[Ni][j][k]);

			NetYMomentumFlux += iFacesFlux[0][j][k][2]*Magnitude(iFaceAreaVector[0][j][k]);
			NetYMomentumFlux -= iFacesFlux[Ni][j][k][2]*Magnitude(iFaceAreaVector[Ni][j][k]);

			NetZMomentumFlux += iFacesFlux[0][j][k][3]*Magnitude(iFaceAreaVector[0][j][k]);
			NetZMomentumFlux -= iFacesFlux[Ni][j][k][3]*Magnitude(iFaceAreaVector[Ni][j][k]);

			NetEnergyFlux += iFacesFlux[0][j][k][4]*Magnitude(iFaceAreaVector[0][j][k]);
			NetEnergyFlux -= iFacesFlux[Ni][j][k][4]*Magnitude(iFaceAreaVector[Ni][j][k]);
		}
	}

	// j faces
	for (int i = 0; i < Ni; ++i)
	{
		for (int k = 0; k < Nk; ++k)
		{
			// j = 0 and j = Nj
			NetMassFlux	+= jFacesFlux[i][0][k][0]*Magnitude(jFaceAreaVector[i][0][k]);
			NetMassFlux	-= jFacesFlux[i][Nj][k][0]*Magnitude(jFaceAreaVector[i][Nj][k]);

			NetXMomentumFlux += jFacesFlux[i][0][k][1]*Magnitude(jFaceAreaVector[i][0][k]);
			NetXMomentumFlux -= jFacesFlux[i][Nj][k][1]*Magnitude(jFaceAreaVector[i][Nj][k]);

			NetYMomentumFlux += jFacesFlux[i][0][k][2]*Magnitude(jFaceAreaVector[i][0][k]);
			NetYMomentumFlux -= jFacesFlux[i][Nj][k][2]*Magnitude(jFaceAreaVector[i][Nj][k]);

			NetZMomentumFlux += jFacesFlux[i][0][k][3]*Magnitude(jFaceAreaVector[i][0][k]);
			NetZMomentumFlux -= jFacesFlux[i][Nj][k][3]*Magnitude(jFaceAreaVector[i][Nj][k]);

			NetEnergyFlux += jFacesFlux[i][0][k][4]*Magnitude(jFaceAreaVector[i][0][k]);
			NetEnergyFlux -= jFacesFlux[i][Nj][k][4]*Magnitude(jFaceAreaVector[i][Nj][k]);
		}
	}

	// k faces
	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			// j = 0 and j = Nj
			NetMassFlux	+= kFacesFlux[i][j][0][0]*Magnitude(kFaceAreaVector[i][j][0]);
			NetMassFlux	-= kFacesFlux[i][j][Nk][0]*Magnitude(kFaceAreaVector[i][j][Nk]);

			NetXMomentumFlux += kFacesFlux[i][j][0][1]*Magnitude(kFaceAreaVector[i][j][0]);
			NetXMomentumFlux -= kFacesFlux[i][j][Nk][1]*Magnitude(kFaceAreaVector[i][j][Nk]);

			NetYMomentumFlux += kFacesFlux[i][j][0][2]*Magnitude(kFaceAreaVector[i][j][0]);
			NetYMomentumFlux -= kFacesFlux[i][j][Nk][2]*Magnitude(kFaceAreaVector[i][j][Nk]);

			NetZMomentumFlux += kFacesFlux[i][j][0][3]*Magnitude(kFaceAreaVector[i][j][0]);
			NetZMomentumFlux -= kFacesFlux[i][j][Nk][3]*Magnitude(kFaceAreaVector[i][j][Nk]);

			NetEnergyFlux += kFacesFlux[i][j][0][4]*Magnitude(kFaceAreaVector[i][j][0]);
			NetEnergyFlux -= kFacesFlux[i][j][Nk][4]*Magnitude(kFaceAreaVector[i][j][Nk]);
		}
	}
	cout << NetMassFlux << "," << NetXMomentumFlux << "," << 
	NetYMomentumFlux << "," << NetZMomentumFlux << "," << NetEnergyFlux << endl;
}
#endif // boundaryNetflux.h ends here   