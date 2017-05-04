%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS
% 1. 2D, grids data points file in .csv format. Which has (x,y) coordinates of
% all points of 2D plane
% 2. How many points are in x direction and in y direction
% 3. U = [rho, rho*u, rho*v, rho*w, e] at evrey point on 2D plane
% 4. Residual file, which has all five residual with time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
% 1. Will plot the grids point to see the geomatry
% 2. All residial plots and automatically save them   
% 3. Plot density, X-velocity, Y-velocity, Z-velocity, p, T and Mach number(M)
% at evrey point on 2D plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc ;

grids =  csvread('./Results/outputfiles/grids_2D.csv') ;

%plot the grids with ghost cell
% gridswithghost =  csvread('./Results/outputfiles/grids_2D_with_ghost.csv'); 

para =  csvread('./Results/outputfiles/2D_parameters_B.csv') ;
residual = csvread('./Results/outputfiles/Residual.csv') ;

%  Reading the first or grids points 
Nx = grids(1,1) ;
Ny = grids(1,2) ;

%  Reading the first or grids points including ghost cells 
% Nxghost = gridswithghost(1,1) ;
% Nyghost = gridswithghost(1,2) ;

% removing the first line and keeping the grids points coordinates alone 
grids(1,:) = [] ;
% gridswithghost(1,:) = [] ;

% Now re-shaping the each 1-D coloum to 2D matrix  
x = reshape(grids(:,1),[Ny,Nx]) ;
y = reshape(grids(:,2),[Ny,Nx]) ;
z = zeros(Ny,Nx);

% xghost = reshape(gridswithghost(:,1),[Nyghost,Nxghost]) ;
% yghost = reshape(gridswithghost(:,2),[Nyghost,Nxghost]) ;
% zghost = zeros(Nyghost,Nxghost);

if 1
    density = reshape(para(:,1),[Ny,Nx]);
    densityu= reshape(para(:,2),[Ny,Nx]);
    densityv= reshape(para(:,3),[Ny,Nx]);
    densityw= reshape(para(:,4),[Ny,Nx]);
    totalEnergy  = reshape(para(:,5),[Ny,Nx]);
    totalTemperature  = zeros(Ny,Nx);
    totalPressure  = zeros(Ny,Nx);
    for i = 1:Ny

        for j =	1:Nx	
            u(i,j) = densityu(i,j)/density(i,j) ; 
            v(i,j) = densityv(i,j)/density(i,j) ; 
            w(i,j) = densityw(i,j)/density(i,j) ; 
            pressure(i,j) = 0.4*(totalEnergy(i,j) - 0.5*density(i,j)* ...
            (u(i,j)*u(i,j)+v(i,j)*v(i,j)+w(i,j)*w(i,j)));
            temperature(i,j) = pressure(i,j) /(287.14*density(i,j)) ; 
            velocity(i,j) = sqrt(u(i,j)*u(i,j)+v(i,j)*v(i,j)+w(i,j)*w(i,j)) ;
            mach(i,j) = velocity(i,j) / sqrt(1.4*287.14*temperature(i,j)) ;
            totalTemperature(i,j) = temperature(i,j)*(1+0.5*(1.4-1)*mach(i,j)*mach(i,j));
            totalPressure(i,j) = pressure(i,j)*(1+0.5*(1.4-1)*mach(i,j)*mach(i,j))^(3.5);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the dimensinal parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plotting has started, Kullu...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For the geometery 

%% 1. Plotting the grid points
i=1;

% if 1
% h = figure(i) ;
% plot(x,y,'o');
% title('Nozzle geomatry')
% xlabel('x(m)')
% ylabel('y(m)')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(h,'./Results/MATLABPlots/Geomatry_Grid_points','epsc')
% end


i=i+1;
h = figure(i) ;
mesh(x,y,z,'FaceLighting','gouraud','LineWidth',0.3)
title('Nozzle geomatry mesh')
xlabel('x(m)')
ylabel('y(m)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
view(2)
saveas(h,'./Results/MATLABPlots/Geomatry_mesh','epsc')

% i=i+1;
% h = figure(i) ;
% mesh(xghost,yghost,zghost,'FaceLighting','gouraud','LineWidth',0.3)
% title('Nozzle geomatry mesh with ghost cells')
% xlabel('x(m)')
% ylabel('y(m)')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% view(2)
% saveas(h,'./Results/MATLABPlots/Geomatry_mesh_with_ghost','epsc')

% i=i+1;
% h = figure(i) ;
% plot(xghost,yghost,'o');
% title('Nozzle geomatry point swith ghost cells')
% xlabel('x(m)')
% ylabel('y(m)')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(h,'./Results/MATLABPlots/Geomatry_grid_points_with_ghost','epsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. plotting the density, velocity, mach, p, T etc. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
i=i+1;
figure(i)
h = surf(x,y,z,mach) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf {Mach}') ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Mach Number (M), flow inside nozzle')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Mach','epsc')
end 

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
h = surf(x,y,z,density) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'$rho$') ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Density(rho), flow inside nozzle')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Density','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
h = surf(x,y,z,velocity) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'$V(m/s)$') ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Velocity(V), flow inside nozzle')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/velocity','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
h = surf(x,y,z,pressure) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'$p(N/m^2)$') ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Pressure(p), flow inside nozzle')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Pressure','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h =figure(i);
surf(x,y,z,temperature) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'$T(K)$') ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Temperature (T), flow inside nozzle')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Temperature','epsc')
end

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h =figure(i);
surf(x,y,z,totalTemperature) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'$T_0(K)$') ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf TotalTemperature (T_0), flow inside nozzle')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/TotalTemperature','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h =figure(i);
surf(x,y,z,totalPressure) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'$p_0(K)$') ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf TotalPressure (P_0), flow inside nozzle')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/TotalPressure','epsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Residuals with respect to time steps
% Density residual
if 1
i=i+1;
h=figure(i);
semilogy(residual(:,1), residual(:,3),'-','LineWidth',1);
title('Density Residual')
xlabel('No. of iterations')
ylabel('Density residual')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Density_residual','epsc')
end

if 1
% x-momentum residual
i=i+1;
h=figure(i);
semilogy(residual(:,1), residual(:,4),'-','LineWidth',1);
title('x-momentum Residual')
xlabel('No. of iterations')
ylabel('x-momentum residual')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/X-momentun_residual','epsc')

% y-momentum residual
i=i+1;
h=figure(i);
semilogy(residual(:,1), residual(:,5),'-','LineWidth',1);
title('y-momentum Residual')
xlabel('No. of iterations')
ylabel('y-momentum residual')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Y-momentun_residual','epsc')

% z-momentum residual
i=i+1;
h=figure(i);
semilogy(residual(:,1), residual(:,6),'-','LineWidth',1);
title('z-momentum Residual')
xlabel('No. of iterations')
ylabel('z-momentum residual')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Z-momentun_residual','epsc')

% energy residual
i=i+1;
h=figure(i);
semilogy(residual(:,1), residual(:,7),'-','LineWidth',1);
title('Energy Residual')
xlabel('No. of iterations')
ylabel('Energy residual')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Energy_residual','epsc')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear;
% clc ;
disp('Results plotting is over, Kullu... :)')
close all;

% 
% GNU plot 
% residual (right residual by 3 deced)

% Nozzle 

% Thrust clc

% field

% Surface

% Integreted properties


% Suggestion 
%     priorties vise

%          