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
% GIS = csvread('../FinalResults/GridIndipendence.csv');


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
        for j = 1:Nx    
            u(i,j) = densityu(i,j)/density(i,j) ; 
            v(i,j) = densityv(i,j)/density(i,j) ; 
            w(i,j) = densityw(i,j)/density(i,j) ; 
            
            Gamma(i,j) = getgamma(density(i,j),u(i,j),v(i,j),w(i,j),totalEnergy(i,j));
            % Gamma(i,j) = 1.4;
            pressure(i,j) = (Gamma(i,j)-1)*(totalEnergy(i,j) - 0.5*density(i,j)* ...
            (u(i,j)*u(i,j)+v(i,j)*v(i,j)+w(i,j)*w(i,j)));
            temperature(i,j) = pressure(i,j) /(287.14*density(i,j)) ; 
            velocity(i,j) = sqrt(u(i,j)*u(i,j)+v(i,j)*v(i,j)+w(i,j)*w(i,j)) ;
            mach(i,j) = velocity(i,j) / sqrt(Gamma(i,j)*287.14*temperature(i,j)) ;
            totalTemperature(i,j) = temperature(i,j)*(1+0.5*(Gamma(i,j)-1)*mach(i,j)*mach(i,j));
            totalPressure(i,j) = pressure(i,j)*(1+0.5*(Gamma(i,j)-1)*mach(i,j)*mach(i,j))^(Gamma(i,j)/(Gamma(i,j)-1));
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
daspect([1 1 1])
title('\bf{Geomatry}')
xlabel('\bf x(m)')
ylabel('\bf y(m)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_Geomatry_CellCenters','epsc')

i=i+1;
h = figure(i) ;
set(gca,'fontsize',18)
hold on
daspect([1 1 1])
mesh(x,y,z,'FaceLighting','gouraud','LineWidth',0.3)
title('\bf Geomatry mesh')
xlabel('\bf x(m)')
ylabel('\bf y(m)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
view(2)
saveas(h,'./Results/MATLABPlots/Gamma_Geomatry_mesh','epsc')
end

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
h = surf(x,y,z,Gamma) ;
daspect([1 1 1])
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf {Gamma}','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Specific Heat Ratio (\gamma), \gamma = \gamma(T) ')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_Mach','epsc')
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
h = surf(x,y,z,mach) ;
daspect([1 1 1])
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf {Mach}','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Mach Number (M), \gamma = \gamma(T) ')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_Mach','epsc')
end 

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
h = surf(x,y,z,density) ;
daspect([1 1 1])
colormap jet
colorbar
hcb=colorbar;
shading interp;
set(t,'Interpreter','Latex');
t = title(hcb,'\rho','fontsize',18) ;
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf{Density(\rho)}, \gamma = \gamma(T) ')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_Density','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
h = surf(x,y,z,velocity) ;
daspect([1 1 1])
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'$V(m/s)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Velocity(V),  \gamma = \gamma(T) ')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_Velocity','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
h = surf(x,y,z,pressure) ;
daspect([1 1 1])
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $p(N/m^2)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Pressure(p), \gamma = \gamma(T) ')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_Pressure','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i);
set(gca,'fontsize',18)
hold on
surf(x,y,z,temperature) ;
daspect([1 1 1])
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $T(K)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf Temperature (T), \gamma = \gamma(T) ')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_Temperature','epsc')
end

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h =figure(i);
set(gca,'fontsize',18)
hold on
surf(x,y,z,totalTemperature) ;
daspect([1 1 1])
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $T_0(K)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf TotalTemperature (T_0), \gamma = \gamma(T) ')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_TotalTemperature','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i);
set(gca,'fontsize',18)
hold on
surf(x,y,z,totalPressure) ;
daspect([1 1 1])
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $p_0(K)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)'); zlabel('\bf z(m)');
title(' \bf TotalPressure (p_0), \gamma = \gamma(T) ')
% view(0,90)
view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_TotalPressure','epsc')
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
C = contour(x,y,Gamma,100);
daspect([1 1 1])
% % clabel(C)
% hold on;
% plot(x(end,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(1,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf Gamma(\gamma)','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf Specific heat ratio(\gamma)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/GammaContourPlot','epsc')
end 

if 1
i=i+1;
h = figure(i)
set(gca,'fontsize',18)
hold on
% h = contour(x,y,mach,'ShowText','on');
C = contour(x,y,mach,100);
daspect([1 1 1])
% % clabel(C)
% hold on;
% plot(x(end,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(1,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf {Mach}','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf Mach Number (M), \gamma = \gamma(T) ')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_MachContourContourPlot','epsc')
end 

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i)
set(gca,'fontsize',18)
hold on
C = contour(x,y,density,100) ;
daspect([1 1 1])
% % clabel(C)
% hold on;
% plot(x(end,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(1,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $\rho$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf{Density(\rho)}, \gamma = \gamma(T) ')
view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_DensityContourPlot','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
C = contour(x,y,velocity,100) ;
daspect([1 1 1])
% % clabel(C)
% hold on;
% plot(x(end,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(1,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $V(m/s)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf Velocity(V), \gamma = \gamma(T) ')
% view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_VelocityContourPlot','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
figure(i)
set(gca,'fontsize',18)
hold on
C = contour(x,y,pressure,100) ;
daspect([1 1 1])
% % clabel(C)
% hold on;
% plot(x(end,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(1,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $p(N/m^2)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf Pressure(p), \gamma = \gamma(T) ')
% view(0,90)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_PressureContourPlot','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i);
set(gca,'fontsize',18)
hold on
% % clabel(C)
% hold on;
% plot(x(end,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(1,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(end,:),y(1,:),'g','LineWidth',3);

contour(x,y,temperature,100) ;
daspect([1 1 1])
colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf T(K)','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf Temperature (T), \gamma = \gamma(T) ')
% view(0,90)
% view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_TemperatureContourPlot','epsc')
end

if 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i=i+1;
% h = figure(i);
% set(gca,'fontsize',18)
% hold on
% C = contour(x,y,totalTemperature,100) ;
% % clabel(C)
% hold on;
% plot(x(end,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(1,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(end,:),y(1,:),'g','LineWidth',3);

% colormap jet
% colorbar
% hcb=colorbar;
% shading interp;
% t = title(hcb,'\bf $T_0(K)$','fontsize',18) ;
% set(t,'Interpreter','Latex');
% xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
% title(' \bf TotalTemperature (T_0), \gamma = \gamma(T) ')
% % view(0,90)
% % view(2)
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(h,'./Results/MATLABPlots/Gamma_TotalTemperatureContourPlot','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
h = figure(i);
set(gca,'fontsize',18)
hold on
C = contour(x,y,totalPressure,100) ;
daspect([1 1 1])
% % clabel(C)
% hold on;
% plot(x(end,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(1,:),y(end,:),'g','LineWidth',3);
% hold on;
% plot(x(end,:),y(1,:),'g','LineWidth',3);

colormap jet
colorbar
hcb=colorbar;
shading interp;
t = title(hcb,'\bf $p_0(K)$','fontsize',18) ;
set(t,'Interpreter','Latex');
xlabel('\bf{x(m)}'); ylabel('\bf y(m)');
title(' \bf TotalPressure (p_0), \gamma = \gamma(T) ')
% view(0,90)
% view(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./Results/MATLABPlots/Gamma_TotalPressureContourPlot','epsc')
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
saveas(h,'./Results/MATLABPlots/Gamma_Density_residual','epsc')
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
saveas(h,'./Results/MATLABPlots/Gamma_X-momentun_residual','epsc')

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
saveas(h,'./Results/MATLABPlots/Gamma_Y-momentun_residual','epsc')

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
saveas(h,'./Results/MATLABPlots/Gamma_Z-momentun_residual','epsc')

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
saveas(h,'./Results/MATLABPlots/Gamma_Energy_residual','epsc')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear;
% clc ;
disp('Results plotting is over, Kullu... :)')
% close all;

% if 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Grid Independence study
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i = i+1;
% h = figure(i)
% plot(GIS(1,:),GIS(2,:),'o-g');
% hold on;
% plot(GIS(1,:),GIS(3,:),'o-b');
% hold on;
% plot(GIS(1,:),GIS(4,:),'o-r');
% hold on;
% plot(GIS(1,:),GIS(5,:),'o-c');
% legend('Mach','Density','Temperature','Pressure')
% title('\bf Grid Independence Study','fontsize',18)
% xlabel('\bf Grid Points','fontsize',18)
% ylabel('\bf Error %','fontsize',18)
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(h,'./Results/MATLABPlots/GridIndipendence','epsc')
% end