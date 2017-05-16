#ifndef BOUNDARYNETFLUX_H
#define BOUNDARYNETFLUX_H
#include <math.h>
#include "vector"
double Magnitude(vector<double> v)
{
	return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

void boundaryNetflux(
	vector<double> & iNetFlux,
	vector<double> & iNiNetFlux,
	vector<double> & jNetFlux,
	vector<double> & jNjNetFlux,
	vector<double> & kNetFlux,
	vector<double> & kNkNetFlux,
	vector<vector<vector<vector<double> > > >	& iFacesFlux,
	vector<vector<vector<vector<double> > > >	& jFacesFlux,
	vector<vector<vector<vector<double> > > >	& kFacesFlux,
	vector<vector<vector<vector<double> > > >	iFaceAreaVector,
	vector<vector<vector<vector<double> > > >	jFaceAreaVector,
	vector<vector<vector<vector<double> > > >	kFaceAreaVector,
	int Ni, int Nj, int Nk)
{
	// vector<double> iNetFlux(5);
	// vector<double> iNiNetFlux(5);
	// vector<double> jNetFlux(5);
	// vector<double> jNjNetFlux(5);
	// vector<double> kNetFlux(5);
	// vector<double> kNkNetFlux(5);

	
	// i faces
	for (int j = 0; j < Nj; ++j)
	{
		for (int k = 0; k < Nk; ++k)
		{
			for (int l = 0; l < 5; ++l)
			{
				// i = 0 
				iNetFlux[l]	+= iFacesFlux[0][j][k][l]*Magnitude(iFaceAreaVector[0][j][k]);
				//and i = Ni
				iNiNetFlux[l] -= iFacesFlux[Ni][j][k][l]*Magnitude(iFaceAreaVector[Ni][j][k]);
			}
		}
	}

	// j faces
	for (int i = 0; i < Ni; ++i)
	{
		for (int k = 0; k < Nk; ++k)
		{
			for (int l = 0; l < 5; ++l)
			{
				// j = 0 
				jNetFlux[l]	+= jFacesFlux[i][0][k][l]*Magnitude(jFaceAreaVector[i][0][k]);
				//and j = Nj
				jNjNetFlux[l] -= jFacesFlux[i][Nj][k][l]*Magnitude(jFaceAreaVector[i][Nj][k]);
			}
		}
	}

	// k faces
	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			for (int l = 0; l < 5; ++l)
			{
				// j = 0 
				kNetFlux[l]	+= kFacesFlux[i][j][0][l]*Magnitude(kFaceAreaVector[i][j][0]);
				//and k = Nk
				kNkNetFlux[l] -= kFacesFlux[i][j][Nk][l]*Magnitude(kFaceAreaVector[i][j][Nk]);
			}
		}
	}
	vector<double> NetFlux(5);
	for (int l = 0; l < 5; ++l)
	{
		NetFlux[l] = iNetFlux[l]+iNiNetFlux[l]+jNetFlux[l]+jNjNetFlux[l]+kNetFlux[l]+kNkNetFlux[l];
	}
	// cout <<iNetFlux[2]<<","<<iNiNetFlux[2]<<","<<jNetFlux[2]<<
	// ","<<jNjNetFlux[2]<<","<<kNetFlux[2]<<","<<kNkNetFlux[2] << endl;
}
#endif // boundaryNetflux.h ends here   