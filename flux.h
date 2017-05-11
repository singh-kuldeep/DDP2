#ifndef FLUX_H
#define FLUX_H

#include "vector"
#include "netfluxAUSM.h" // AUSM
#include "BC.h"

// this function calculates the flux at every face and stores it 
void flux(
	vector<vector<vector<vector<double> > > >	& iFacesFlux,
	vector<vector<vector<vector<double> > > >	& jFacesFlux,
	vector<vector<vector<vector<double> > > >	& kFacesFlux,
	vector<vector<vector<vector<double> > > >	iFaceAreaVector,
	vector<vector<vector<vector<double> > > >	jFaceAreaVector,
	vector<vector<vector<vector<double> > > >	kFaceAreaVector,
	vector<vector<vector<vector<double> > > >	ConservedVariables,
	int Ni, int Nj, int Nk)
{
	// Creating a 4D vector object
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> Dim4;

	std::vector<double> LeftConservedVariables(5);
	std::vector<double> RightConservedVariables(5);
	// i faces calculation

	/**\param i0GhostConservedVariable Ghost cell Conserved variables array at i = 0*/
	/**\param j0GhostConservedVariable Ghost cell Conserved variables array at j = 0*/
	/**\param k0GhostConservedVariable Ghost cell Conserved variables array at k = 0*/
	/**\param iNiGhostConservedVariable Ghost cell Conserved variables array at i = Ni*/
	/**\param jNjGhostConservedVariable Ghost cell Conserved variables array at j = Nj*/
	/**\param kNkGhostConservedVariable Ghost cell Conserved variables array at k = Nk*/
	Dim4 i0GhostConservedVariable(1,Dim3(Nj,Dim2(Nk,Dim1(5))));
	Dim4 iNiGhostConservedVariable(1,Dim3(Nj,Dim2(Nk,Dim1(5))));
	Dim4 j0GhostConservedVariable(Ni,Dim3(1,Dim2(Nk,Dim1(5))));
	Dim4 jNjGhostConservedVariable(Ni,Dim3(1,Dim2(Nk,Dim1(5))));
	Dim4 k0GhostConservedVariable(Ni,Dim3(Nj,Dim2(1,Dim1(5))));
	Dim4 kNkGhostConservedVariable(Ni,Dim3(Nj,Dim2(1,Dim1(5))));

	/* Before every time step we need to have proper value in the ghost 
		cells So, BC takes  care of Inlet, Exit, y-wall and Z-wall BC */
	BC(ConservedVariables,
	iFaceAreaVector,jFaceAreaVector,kFaceAreaVector,
	i0GhostConservedVariable,j0GhostConservedVariable,
	k0GhostConservedVariable,iNiGhostConservedVariable,
	jNjGhostConservedVariable,kNkGhostConservedVariable,
	Ni, Nj, Nk);
	
	#if 0 // BC() function has assigned the values correctly
	if(test4DArray("i0GhostConservedVariable",i0GhostConservedVariable,1,Nj,Nk,5)==0)
	{
		// return 0;
	}
	if(test4DArray("iNiGhostConservedVariable",iNiGhostConservedVariable,1,Nj,Nk,5)==0)
	{
		// return 0;
	}
	if(test4DArray("j0GhostConservedVariable",j0GhostConservedVariable,Ni,1,Nk,5)==0)
	{
		// return 0;
	} 
	if(test4DArray("jNjGhostConservedVariable",jNjGhostConservedVariable,Ni,1,Nk,5)==0)
	{
		// return 0;
	}
	if(test4DArray("k0GhostConservedVariable",k0GhostConservedVariable,Ni,Nj,1,5)==0)
	{
		// return 0;
	}
	if(test4DArray("kNkGhostConservedVariable",kNkGhostConservedVariable,Ni,Nj,1,5)==0)
	{
		// return 0;
	}
	#endif

	// i faces calculation
	for (int i = 0; i < Ni+1; ++i)
	{

		for (int j = 0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{	
				// i interface volume
				if(i == 0)
				{
					// iCellInterfaceVolume = 0.5*(i0GhostCellVolume[0][j][k] + 
					// 	CellVolume[i][j][k]);  
					LeftConservedVariables = i0GhostConservedVariable[0][j][k];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the -0.5 interface 
					netfluxAUSM irightface(LeftConservedVariables,
					RightConservedVariables,
					iFaceAreaVector[i][j][k]);


					for (int l = 0; l < 5; ++l)
					{
						if(isnan(irightface.NetFlux[l])==1)
						{
							// cout << "irightface.NetFlux["<<l<<"] is NaN at 
							// [" << i << "," << j << " " << k <<"]" << endl;
							// cout << "Check the line number "<< __LINE__ <<
							// " in TVDMainSolver.cpp" << endl;
							// return 0;
						}
						iFacesFlux[i][j][k][l] = irightface.NetFlux[l]; 
						// ConservedVariablesNew[i][j][k][l] +=(deltat/
						// 	iCellInterfaceVolume)*(irightface.NetFlux[l]);
					}
				}
				#if 1
				else if(i == Ni)
				{
					// iCellInterfaceVolume = 0.5*(CellVolume[i-1][j][k] +
					// 	iNiGhostCellVolume[0][j][k]);
					LeftConservedVariables = ConservedVariables[i-1][j][k];
					RightConservedVariables = iNiGhostConservedVariable[0][j][k];
					
					// Calculating the flux at the Ni-0.5 interface 
					netfluxAUSM irightface(LeftConservedVariables,
					RightConservedVariables,
					iFaceAreaVector[i][j][k]);
					
					for (int l = 0; l < 5; ++l)
					{
						if(isnan(irightface.NetFlux[l])==1)
						{
							// cout << "irightface.NetFlux["<<l<<"] is NaN at 
							// [" << i << "," << j << "," << k <<"]" << endl;
							// cout << "Check the line number "<< __LINE__ << 
							// " in TVDMainSolver.cpp" << endl;
							// return 0;
						}
						iFacesFlux[i][j][k][l] = irightface.NetFlux[l]; 
						// ConservedVariablesNew[Ni-1][j][k][l] -=(deltat/
							// iCellInterfaceVolume)*(irightface.NetFlux[l]);
					}
					
				}
				else
				{
					// iCellInterfaceVolume = 0.5*(CellVolume[i-1][j][k] + 
					// CellVolume[i][j][k]);
					LeftConservedVariables = ConservedVariables[i-1][j][k];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the i-0.5 interface 
					netfluxAUSM irightface(LeftConservedVariables,
					RightConservedVariables,
					iFaceAreaVector[i][j][k]);

					for (int l = 0; l < 5; ++l)
					{
						if(isnan(irightface.NetFlux[l])==1)
						{
							// cout << "irightface.NetFlux["<<l<<"] is NaN at 
							// [" << i << "," << j << "," << k <<"]" << endl;
							// cout << "Check the line number "<< __LINE__ << 
							// " in TVDMainSolver.cpp" << endl;
							// return 0;
						}
						iFacesFlux[i][j][k][l] = irightface.NetFlux[l]; 
						// ConservedVariablesNew[i-1][j][k][l] -=(deltat/
						// 	iCellInterfaceVolume)*(irightface.NetFlux[l]);
						// ConservedVariablesNew[i][j][k][l] +=(deltat/
						// 	iCellInterfaceVolume)*(irightface.NetFlux[l]);
					}

				}
				#endif
			}
		}
	}

	// j faces calculation
	for (int i = 0; i < Ni; ++i)
	{
		for (int j =0; j < Nj+1; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{	
				if(j == 0)
				{
					// jCellInterfaceVolume = 0.5*(j0GhostCellVolume[i][0][k] + 
					// 	CellVolume[i][j][k]);
					LeftConservedVariables = j0GhostConservedVariable[i][0][k];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the -0.5 interface 
					netfluxAUSM jrightface(LeftConservedVariables,
					RightConservedVariables,
					jFaceAreaVector[i][j][k]);

					for (int l = 0; l < 5; ++l)
					{
						if(isnan(jrightface.NetFlux[l])==1)
						{
							// cout << "jrightface.NetFlux["<<l<<"] is NaN at 
							// [" << i << "," << j << "," << k <<"]" << endl;
							// cout << "Check the line number "<< __LINE__ << 
							// " in TVDMainSolver.cpp" << endl;
							// return 0;
						}
						jFacesFlux[i][j][k][l] = jrightface.NetFlux[l]; 

						// ConservedVariablesNew[i][j][k][l] +=(deltat/
						// 	jCellInterfaceVolume)*(jrightface.NetFlux[l]);
					}
				}
				else if(j == Nj)
				{
					// jCellInterfaceVolume = 0.5*(CellVolume[i][j-1][k] +
					// jNjGhostCellVolume[i][0][k]);
					LeftConservedVariables = ConservedVariables[i][j-1][k];
					RightConservedVariables = jNjGhostConservedVariable[i][0][k];

					// Calculating the flux at the Nj-0.5 interface 
					netfluxAUSM jrightface(LeftConservedVariables,
					RightConservedVariables,
					jFaceAreaVector[i][j][k]);

					for (int l = 0; l < 5; ++l)
					{
						if(isnan(jrightface.NetFlux[l])==1)
						{
							// cout << "jrightface.NetFlux["<<l<<"] is NaN at 
							// [" << i << "," << j << "," << k <<"]" << endl;
							// cout << "Check the line number "<< __LINE__ << 
							// " in TVDMainSolver.cpp" << endl;
							// return 0;
						}
						jFacesFlux[i][j][k][l] = jrightface.NetFlux[l]; 

						// ConservedVariablesNew[i][j-1][k][l] -=(deltat/
						// 	jCellInterfaceVolume)*(jrightface.NetFlux[l]);
					}
				}
				else
				{
					// jCellInterfaceVolume = 0.5*(CellVolume[i][j-1][k] + 
					// CellVolume[i][j][k]);
					LeftConservedVariables = ConservedVariables[i][j-1][k];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the j-0.5 interface 
					netfluxAUSM jrightface(LeftConservedVariables,
					RightConservedVariables,
					jFaceAreaVector[i][j][k]);

					for (int l = 0; l < 5; ++l)
					{
						if(isnan(jrightface.NetFlux[l])==1)
						{
							// cout << "jrightface.NetFlux["<<l<<"] is NaN at 
							// [" << i << "," << j << "," << k <<"]" << endl;
							// cout << "Check the line number "<< __LINE__ << 
							// " in TVDMainSolver.cpp" << endl;
							// return 0;
						}
						jFacesFlux[i][j][k][l] = jrightface.NetFlux[l]; 

						// ConservedVariablesNew[i][j-1][k][l] -=(deltat/
						// 	jCellInterfaceVolume)*(jrightface.NetFlux[l]);
						// ConservedVariablesNew[i][j][k][l] +=(deltat/
						// 	jCellInterfaceVolume)*(jrightface.NetFlux[l]);
					}
				}

			}
		}
	}
	
	
	for (int i = 0; i < Ni; ++i)
	{
		for (int j =0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk+1; ++k)
			{
				// k interface volume
				if(k == 0)
				{						
					// kCellInterfaceVolume = 0.5*(k0GhostCellVolume[i][j][0] + 
					// 	CellVolume[i][j][k]);
					LeftConservedVariables = k0GhostConservedVariable[i][j][0];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the -0.5 interface 
					netfluxAUSM krightface(LeftConservedVariables,
					RightConservedVariables,
					kFaceAreaVector[i][j][k]);

					for (int l = 0; l < 5; ++l)
					{
						if(isnan(krightface.NetFlux[l])==1)
						{
							// cout << "krightface.NetFlux["<<l<<"] is NaN at 
							// [" << i << "," << j << "," << k <<"]" << endl;
							// cout << "Check the line number "<< __LINE__ << 
							// " in TVDMainSolver.cpp" << endl;
							// return 0;
						}
						kFacesFlux[i][j][k][l] = krightface.NetFlux[l]; 

						// ConservedVariablesNew[i][j][k][l] +=(deltat/
						// 	kCellInterfaceVolume)*(krightface.NetFlux[l]);
					}
				}
				
				else if(k == Nk)
				{
					// kCellInterfaceVolume = 0.5*(CellVolume[i][j][k-1]+
					// 	kNkGhostCellVolume[i][j][0]);
					LeftConservedVariables = ConservedVariables[i][j][k-1];
					RightConservedVariables = kNkGhostConservedVariable[i][j][0];

					// Calculating the flux at the Nk-0.5 interface 
					netfluxAUSM krightface(LeftConservedVariables,
					RightConservedVariables,
					kFaceAreaVector[i][j][k]);

					for (int l = 0; l < 5; ++l)
					{
						if(isnan(krightface.NetFlux[l])==1)
						{
							// cout << "krightface.NetFlux["<<l<<"] is NaN at 
							// [" << i << "," << j << "," << k <<"]" << endl;
							// cout << "Check the line number "<< __LINE__ << 
							// " in TVDMainSolver.cpp" << endl;
							// return 0;
						}
						kFacesFlux[i][j][k][l] = krightface.NetFlux[l]; 
						
						// ConservedVariablesNew[i][j][k-1][l] -=(deltat/
						// 	kCellInterfaceVolume)*(krightface.NetFlux[l]);
					}
				}
				else
				{
					// kCellInterfaceVolume = 0.5*(CellVolume[i][j][k-1] + 
					// CellVolume[i][j][k]);
					LeftConservedVariables = ConservedVariables[i][j][k-1];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the j-0.5 interface 
					netfluxAUSM krightface(LeftConservedVariables,
					RightConservedVariables,
					kFaceAreaVector[i][j][k]);

					for (int l = 0; l < 5; ++l)
					{
						if(isnan(krightface.NetFlux[l])==1)
						{
							// cout << "krightface.NetFlux["<<l<<"] is NaN at 
							// [" << i << "," << j << "," << k <<"]" << endl;
							// cout << "Check the line number "<< __LINE__ << 
							// " in TVDMainSolver.cpp" << endl;
							// return 0;
						}
						kFacesFlux[i][j][k][l] = krightface.NetFlux[l]; 

						// ConservedVariablesNew[i][j][k-1][l] -=(deltat/
						// 	kCellInterfaceVolume)*(krightface.NetFlux[l]);
						// ConservedVariablesNew[i][j][k][l] +=(deltat/
						// 	kCellInterfaceVolume)*(krightface.NetFlux[l]);
					}
				}
			}
		}
	}
}

#endif // FLUX_H ends here