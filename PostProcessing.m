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

grids =  csvread('./Results/outputfiles/CellCenter_ij.csv') ;

para =  csvread('./Results/outputfiles/ConservedQuantity.csv') ;
residual = csvread('./Results/outputfiles/Residual.csv') ;

%  Reading the first or grids points 
Nx = grids(1,1) ;
Ny = grids(1,2) ;

% removing the first line and keeping the grids points coordinates alone 
grids(1,:) = [] ;

% Now re-shaping the each 1-D coloum to 2D matrix  
x = reshape(grids(:,1),[Ny,Nx]) ;
y = reshape(grids(:,2),[Ny,Nx]) ;
z = zeros(Ny,Nx);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Plotting the grid points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;

if 1
h = figure(i);
set(gca,'fontsize',18)
hold on
plot(x,y,'o');
title('\bf{Geomatry}')
xlabel('\bf x(m)')
ylabel('\bf y(m)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Geomatry_CellCenters','epsc')
end


i=i+1;
h = figure(i) ;
set(gca,'fontsize',18)
hold on
mesh(x,y,z,'FaceLighting','gouraud','LineWidth',0.3)
title('\bf Geomatry mesh')
xlabel('\bf x(m)')
ylabel('\bf y(m)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
view(2)
saveas(h,'./Results/MATLABPlots/Geomatry_mesh','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Surface plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
h = surf(x,y,z,mach) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf {Mach}','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Mach Number (M), \gamma = 1.4(Constant) ')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Mach','epsc')
end 

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
h = surf(x,y,z,density) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
set(t,'Interpreter','Latex');
t = title(hcb,'\rho','fontsize',18) ;
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf{Density(\rho)}, \gamma = 1.4(Constant)')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Density','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
h = surf(x,y,z,velocity) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'$V(m/s)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Velocity(V),  \gamma = 1.4(Constant) ')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Velocity','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
h = surf(x,y,z,pressure) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $p(N/m^2)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Pressure(p), \gamma = 1.4(Constant) ')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Pressure','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i);
set(gca,'fontsize',18)
hold on
surf(x,y,z,temperature) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $T(K)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Temperature (T), \gamma = 1.4(Constant) ')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Temperature','epsc')
end

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h =figure(i);
set(gca,'fontsize',18)
hold on
surf(x,y,z,totalTemperature) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $T_0(K)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf TotalTemperature (T_0), \gamma = 1.4(Constant) ')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/TotalTemperature','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i);
set(gca,'fontsize',18)
hold on
surf(x,y,z,totalPressure) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $p_0(K)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf TotalPressure (p_0), \gamma = 1.4(Constant) ')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/TotalPressure','epsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Contoure Plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
i=i+1;
h = figure(i)
set(gca,'fontsize',18)
hold on
% h = contour(x,y,mach,'ShowText','on');
C = contour(x,y,mach,100);
% clabel(C)
hold on;
plot(x(end,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(1,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf {Mach}','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf Mach Number (M), \gamma = 1.4(Constant) ')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/MachContourContourPlot','epsc')
end 

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i)
set(gca,'fontsize',18)
hold on
C = contour(x,y,density,100) ;
% clabel(C)
hold on;
plot(x(end,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(1,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $\rho$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf{Density(\rho)}, \gamma = 1.4(Constant) ')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/DensityContourPlot','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
C = contour(x,y,velocity,100) ;
% clabel(C)
hold on;
plot(x(end,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(1,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $V(m/s)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf Velocity(V), \gamma = 1.4(Constant) ')
% view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/VelocityContourPlot','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
C = contour(x,y,pressure,100) ;
% clabel(C)
hold on;
plot(x(end,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(1,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $p(N/m^2)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf Pressure(p), \gamma = 1.4(Constant) ')
% view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/PressureContourPlot','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i);
set(gca,'fontsize',18)
hold on
% clabel(C)
hold on;
plot(x(end,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(1,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(end,:),y(1,:),'g','LineWidth',3);

contour(x,y,temperature,100) ;
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf T(K)','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf Temperature (T), \gamma = 1.4(Constant) ')
% view(0,90)
% view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/TemperatureContourPlot','epsc')
end

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i);
set(gca,'fontsize',18)
hold on
C = contour(x,y,totalTemperature,100) ;
% clabel(C)
hold on;
plot(x(end,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(1,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $T_0(K)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf TotalTemperature (T_0), \gamma = 1.4(Constant) ')
% view(0,90)
% view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/TotalTemperatureContourPlot','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i);
set(gca,'fontsize',18)
hold on
C = contour(x,y,totalPressure,100) ;
% clabel(C)
hold on;
plot(x(end,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(1,:),y(end,:),'g','LineWidth',3);
hold on;
plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $p_0(K)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf TotalPressure (p_0), \gamma = 1.4(Constant) ')
% view(0,90)
% view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/TotalPressureContourPlot','epsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Residuals with respect to time steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Density residual
if 1
i=i+1;
h=figure(i);
% set(gca,'fontsize',18)
% hold on
semilogy(residual(:,1), residual(:,2),'-','LineWidth',1);
title('\bf Density Residual','fontsize',18)
xlabel('\bf No. of iterations','fontsize',18)
ylabel('\bf Density Residual','fontsize',18)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Density_residual','epsc')
end

if 1
% x-momentum residual
i=i+1;
h=figure(i);
% set(gca,'fontsize',18)
% hold on
semilogy(residual(:,1), residual(:,3),'-','LineWidth',1);
title('\bf x-momentum Residual','fontsize',18)
xlabel('\bf No. of iterations','fontsize',18)
ylabel('\bf x-momentum Residual','fontsize',18)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/X-momentun_residual','epsc')

% y-momentum residual
i=i+1;
h=figure(i);
% set(gca,'fontsize',18)
% hold on
semilogy(residual(:,1), residual(:,4),'-','LineWidth',1);
title('\bf y-momentum Residual','fontsize',18)
xlabel('\bf No. of iterations','fontsize',18)
ylabel('\bf y-momentum Residual','fontsize',18)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Y-momentun_residual','epsc')

% z-momentum residual
i=i+1;
h=figure(i);
% set(gca,'fontsize',18)
% hold on
semilogy(residual(:,1), residual(:,5),'-','LineWidth',1);
title('\bf z-momentum Residual','fontsize',18)
xlabel('\bf No. of iterations','fontsize',18)
ylabel('\bf z-momentum Residual','fontsize',18)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Z-momentun_residual','epsc')

% energy residual
i=i+1;
h=figure(i);
% set(gca,'fontsize',18)
% hold on
semilogy(residual(:,1), residual(:,6),'-','LineWidth',1);
title('\bf Energy Residual','fontsize',18)
xlabel('\bf No. of iterations','fontsize',18)
ylabel('\bf Energy Residual','fontsize',18)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Energy_residual','epsc')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear;
% clc ;
disp('Results plotting is over, Kullu... :)')
% close all;


% Thrust clc
% Integreted properties