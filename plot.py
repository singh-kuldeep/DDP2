# Post processing of the output data
# 1. Grid
# 2. Residual
# 3. Fluid dynamics parameters


import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# grid 
print('grid file reading')
grid = np.genfromtxt('./Results/outputfiles/CellCenter_ij.csv', delimiter=',')

Nx = int(grid[0,0]) # cells in 'i' direction
Ny = int(grid[0,1]) # cells in 'j' direction

print 'grid points Nx =', Nx, 'and Ny = ', Ny 

# Delete the first row after reading the grid points
grid = np.delete(grid, (0), axis=0)

x = grid[:,0]; # x coordinates of grid 
y = grid[:,1]; # y coordinates of grid 

x = np.reshape(x,(Nx,Ny))
y = np.reshape(y,(Nx,Ny))
z = np.zeros((Nx,Ny)) # zeros because grid is 2D
print('grid file reading over')


# Conserved parameters
print('Conserved parameter(U) reading from the file')
U = np.genfromtxt('./Results/outputfiles/ConservedQuantity.csv', delimiter=',')
# all the parameter are reshaped to make it 2D 
Density = np.reshape(U[:,0],(Nx,Ny))
MomentumX = np.reshape(U[:,1],(Nx,Ny)) 
MomentumY = np.reshape(U[:,2],(Nx,Ny)) 
MomentumZ = np.reshape(U[:,3],(Nx,Ny))
TotalEnergy = np.reshape(U[:,4],(Nx,Ny))

VelocityX = np.divide(MomentumX,Density)
VelocityY = np.divide(MomentumY,Density)
VelocityZ = np.divide(MomentumZ,Density)

Pressure = np.zeros((Nx,Ny))
Temperature = np.zeros((Nx,Ny))
VelocityMagnitude = np.zeros((Nx,Ny))
VelocitySound = np.zeros((Nx,Ny))
Mach = np.zeros((Nx,Ny))
TemperatureStagnation = np.zeros((Nx,Ny))
PressureStagnation = np.zeros((Nx,Ny))
for i in range(0,Nx-1):
	for j in range(0,Ny-1):
		Pressure[i,j] = 0.4*(TotalEnergy[i,j]- 0.5*Density[i,j]*
			(VelocityX[i,j]**2 + VelocityY[i,j]**2 + VelocityZ[i,j]**2));
		Temperature[i,j] = Pressure[i,j]/(287.14*Density[i,j]);
		VelocityMagnitude[i,j] =( VelocityX[i,j]**2 + 
					VelocityY[i,j]**2 + VelocityZ[i,j]**2)**0.5;
		VelocitySound[i,j] = (1.4*287.14*Temperature[i,j])**0.5;	
		Mach[i,j] = VelocityMagnitude[i,j]/VelocitySound[i,j];
		TemperatureStagnation[i,j] = Temperature[i,j]*(1+0.5*0.4*Mach[i,j]**2);
		PressureStagnation[i,j] = Pressure[i,j] * (1+0.5*0.4*Mach[i,j]**2)**3.5;

print('Conserved parameter(U) reading from the file is over')


# Residual
print('Residuals reading from the file')
Residual = np.genfromtxt('./Results/outputfiles/Residual.csv', delimiter=',')

Iterations = Residual[:,0]
DensityResidual = Residual[:,2]
MomentumXResidual = Residual[:,3]
MomentumYResidual = Residual[:,4]
MomentumZResidual = Residual[:,5]
EnergyResidual = Residual[:,6]

print('Residuals reading from the file is over')



print('Results plotting starts ...')
i = 1
################################################################################
# Geomatry Plots
################################################################################
# Mesh
if 1==1:
# if 1==1:
	fig = plt.figure(i)
	ax = fig.add_subplot(111, projection='3d')
	ax.view_init(azim=0, elev=90)
	ax.plot_wireframe(x,y,z, rstride=1, cstride=1)
	ax.view_init(azim=-90, elev=90)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Grids')
	plt.savefig('./Results/PythonPlots/Mesh2D.png')   
	# plt.show()
	plt.close()

# mesh only lines
# if 1==1:
	i = i+1;
	fig = plt.figure(i)
	plt.plot(x,y,'.')
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Grids')
	plt.savefig('./Results/PythonPlots/Mesh2DLine.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

################################################################################
# Surface Plots
################################################################################
# Mach surface
if 1==1:
# if 1==1:
	i = i+1
	fig = plt.figure(i)
	CS = plt.pcolor(x,y,Mach)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Mach')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/MachSurfacePlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

#Density Surface
# if 1==1:
	i = i+1
	fig = plt.figure(i)
	CS = plt.pcolor(x,y,Density)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Density')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/DensitySurfacePlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

# Pressure Surface
# if 1==1:
	i = i+1
	fig = plt.figure(i)
	CS = plt.pcolor(x,y,Pressure)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Pressure($N/m^2$)')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/PressureSurfacePlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

# Velocity surface
# if 1==1:
	i = i+1
	fig = plt.figure(i)
	CS = plt.pcolor(x,y,VelocityMagnitude)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Velocity Magnitude(m/s)')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/VelocityColor.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

# TemperatureStagnation Surface
# if 1==1:
	i = i+1
	fig = plt.figure(i)
	CS = plt.pcolor(x,y,TemperatureStagnation)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('TemperatureStagnation(K)')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/TemperatureStagnationSurfacePlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()	

# Temperature Surface
# if 1==1:
	i = i+1
	fig = plt.figure(i)
	CS = plt.pcolor(x,y,Temperature)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Temperature(K)')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/TemperatureSurfacePlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()	

# PressureStagnation Surface
# if 1==1:
	i = i+1
	fig = plt.figure(i)
	CS = plt.pcolor(x,y,PressureStagnation)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('PressureStagnation($N/m^2$)')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/PressureStagnationSurfacePlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

################################################################################
# Contour Plots
################################################################################
if 1==1:
# Pressure Contour
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	CS = plt.contour(x,y,Pressure,200)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Pressure Contours')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/PressureContourPlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

# Mach contour	
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	CS = plt.contour(x,y,Mach,10)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Mach Contours')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/MachContourPlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

# Density contour
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	CS = plt.contour(x,y,Density,200)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Density Contours')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/DensityContourPlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

# Velocity contour
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	CS = plt.contour(x,y,VelocityMagnitude,200)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Velocity Magnitude Contours')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/VelocityMagnitudeContourPlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

# Temperature Contour
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	CS = plt.contour(x,y,Temperature,200)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('Temperature Contours')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/TemperatureContourPlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

# TemperatureStagnation Contour
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	CS = plt.contour(x,y,TemperatureStagnation,200)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('TemperatureStagnation Contours')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/TemperatureStagnationContourPlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

# PressureStagnation Contour
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	CS = plt.contour(x,y,PressureStagnation,200)
	plt.xlabel('X(m)')
	plt.ylabel('Y(m)')
	plt.title('PressureStagnation Contours')
	plt.colorbar()
	plt.savefig('./Results/PythonPlots/PressureStagnationContourPlot.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()



################################################################################
# Residuals
################################################################################
# Density Residual
if 1==1:
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	fig = plt.semilogy(Iterations,DensityResidual,'-')
	plt.xlabel('No. of Iterations')
	plt.ylabel('Density Residual')
	plt.title('Density Residual')
	plt.savefig('./Results/PythonPlots/DensityResidual.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()

# MomentumX Residual
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	fig = plt.semilogy(Iterations,MomentumXResidual,'-')
	plt.xlabel('No. of Iterations')
	plt.ylabel('MomentumX Residual')
	plt.title('MomentumX Residual')
	plt.savefig('./Results/PythonPlots/MomentumXResidual.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()		

# MomentumY Residual
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	fig = plt.semilogy(Iterations,MomentumYResidual,'-')
	plt.xlabel('No. of Iterations')
	plt.ylabel('MomentumY Residual')
	plt.title('MomentumY Residual')
	plt.savefig('./Results/PythonPlots/MomentumYResidual.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()	

# MomentumZ Residual
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	fig = plt.semilogy(Iterations,MomentumZResidual,'-')
	plt.xlabel('No. of Iterations')
	plt.ylabel('MomentumZ Residual')
	plt.title('MomentumZ Residual')
	plt.savefig('./Results/PythonPlots/MomentumZResidual.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()		

# Energy Residual
# if 1==1:
	i = i+1;
	fig = plt.figure(i+1)
	fig = plt.semilogy(Iterations,EnergyResidual,'-')
	plt.xlabel('No. of Iterations')
	plt.ylabel('Energy Residual')
	plt.title('Energy Residual')
	plt.savefig('./Results/PythonPlots/EnergyResidual.png')   
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	# plt.show()
	plt.close()	





