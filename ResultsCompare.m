% 
% Aim is 
% 1) % change in the average quantities(Mach, Temperature, Pressure) at any give "x"
% 
clear all;
clc ;
close all;
para14 =  csvread('../gamma1_4/Results/outputfiles/ConservedQuantity.csv') ;
paraT14 =  csvread('../gammaT_1_4/Results/outputfiles/ConservedQuantity.csv') ;

grids =  csvread('../gamma1_4/Results/outputfiles/CellCenter_ij.csv') ;

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
    density14 = reshape(para14(:,1),[Ny,Nx]);
    densityu14= reshape(para14(:,2),[Ny,Nx]);
    densityv14= reshape(para14(:,3),[Ny,Nx]);
    densityw14= reshape(para14(:,4),[Ny,Nx]);
    totalEnergy14  = reshape(para14(:,5),[Ny,Nx]);
    totalTemperature14  = zeros(Ny,Nx);
    totalPressure14  = zeros(Ny,Nx);

    for i = 1:Ny
        for j = 1:Nx    
            u14(i,j) = densityu14(i,j)/density14(i,j) ; 
            v14(i,j) = densityv14(i,j)/density14(i,j) ; 
            w14(i,j) = densityw14(i,j)/density14(i,j) ; 
            
            % Gamma(i,j) = getgamma(density14(i,j),u14(i,j),v14(i,j),w14(i,j),totalEnergy14(i,j));
            Gamma(i,j) = 1.4;
            pressure14(i,j) = (Gamma(i,j)-1)*(totalEnergy14(i,j) - 0.5*density14(i,j)* ...
            (u14(i,j)*u14(i,j)+v14(i,j)*v14(i,j)+w14(i,j)*w14(i,j)));
            temperature14(i,j) = pressure14(i,j) /(287.14*density14(i,j)) ; 
            velocity14(i,j) = sqrt(u14(i,j)*u14(i,j)+v14(i,j)*v14(i,j)+w14(i,j)*w14(i,j)) ;
            mach14(i,j) = velocity14(i,j) / sqrt(Gamma(i,j)*287.14*temperature14(i,j)) ;
            totalTemperature14(i,j) = temperature14(i,j)*(1+0.5*(Gamma(i,j)-1)*mach14(i,j)*mach14(i,j));
            totalPressure14(i,j) = pressure14(i,j)*(1+0.5*(Gamma(i,j)-1)*mach14(i,j)*mach14(i,j))^(Gamma(i,j)/(Gamma(i,j)-1));
            KE(i,j) = 0.5*density14(i,j)*(velocity14(i,j)^2);
        end
    end
end


if 1
    densityT14 = reshape(paraT14(:,1),[Ny,Nx]);
    densityuT14= reshape(paraT14(:,2),[Ny,Nx]);
    densityvT14= reshape(paraT14(:,3),[Ny,Nx]);
    densitywT14= reshape(paraT14(:,4),[Ny,Nx]);
    totalEnergyT14  = reshape(paraT14(:,5),[Ny,Nx]);
    totalTemperatureT14  = zeros(Ny,Nx);
    totalPressureT14  = zeros(Ny,Nx);

    for i = 1:Ny
        for j = 1:Nx    
            uT14(i,j) = densityuT14(i,j)/densityT14(i,j) ; 
            vT14(i,j) = densityvT14(i,j)/densityT14(i,j) ; 
            wT14(i,j) = densitywT14(i,j)/densityT14(i,j) ; 
            
            Gamma(i,j) = getgamma(densityT14(i,j),uT14(i,j),vT14(i,j),wT14(i,j),totalEnergyT14(i,j));
            % Gamma(i,j) = 1.4;
            pressureT14(i,j) = (Gamma(i,j)-1)*(totalEnergyT14(i,j) - 0.5*densityT14(i,j)* ...
            (uT14(i,j)*uT14(i,j)+vT14(i,j)*vT14(i,j)+wT14(i,j)*wT14(i,j)));
            temperatureT14(i,j) = pressureT14(i,j) /(287.14*densityT14(i,j)) ; 
            velocityT14(i,j) = sqrt(uT14(i,j)*uT14(i,j)+vT14(i,j)*vT14(i,j)+wT14(i,j)*wT14(i,j)) ;
            machT14(i,j) = velocityT14(i,j) / sqrt(Gamma(i,j)*287.14*temperatureT14(i,j)) ;
            totalTemperatureT14(i,j) = temperatureT14(i,j)*(1+0.5*(Gamma(i,j)-1)*machT14(i,j)*machT14(i,j));
            totalPressureT14(i,j) = pressureT14(i,j)*(1+0.5*(Gamma(i,j)-1)*machT14(i,j)*machT14(i,j))^(Gamma(i,j)/(Gamma(i,j)-1));
            KET(i,j) = 0.5*densityT14(i,j)*(velocityT14(i,j)^2);
        end
    end
end

Xlocation = 1.0;  % Non dimensional "x" location  
TotalLocation = 5;
% delatChange = rand(TotalLocation,6)
% delatChange(1,1) = 'NonDimx' ;
% delatChange(1,2) = 'Mach' ;
% delatChange(1,3) = 'Density' ;
% delatChange(1,4) = 'Velocity' ;
% delatChange(1,5) = 'Temperature' ;
% delatChange(1,6) = 'Pressure' ;
for i=1:TotalLocation
	iLoc = round(Nx*i/TotalLocation);
	% xloc = Xlocation*x(end,end)
	% PercentAvgmachChange = (abs(sum(mach14(:,iLoc))) - abs(sum(machT14(:,iLoc))))*100/abs(sum(mach14(:,iLoc)))
	% PercentAvgdensityChange = (abs(sum(density14(:,iLoc))) - abs(sum(densityT14(:,iLoc))))*100/abs(sum(density14(:,iLoc)))
	% PercentAvgvelocityChange = (abs(sum(velocity14(:,iLoc))) - abs(sum(velocityT14(:,iLoc))))*100/abs(sum(velocity14(:,iLoc)))
	% PercentAvgtemperatureChange = (abs(sum(temperature14(:,iLoc))) - abs(sum(temperatureT14(:,iLoc))))*100/abs(sum(temperature14(:,iLoc)))
	% PercentAvgPressureChange = (abs(sum(pressure14(:,iLoc))) - abs(sum(pressureT14(:,iLoc))))*100/abs(sum(pressure14(:,iLoc)))
	delatChange(i,1) = i/TotalLocation;
	delatChange(i,2) = -(abs(sum(mach14(:,iLoc))) - abs(sum(machT14(:,iLoc))))*100/abs(sum(mach14(:,iLoc)));
	delatChange(i,3) = -(abs(sum(density14(:,iLoc))) - abs(sum(densityT14(:,iLoc))))*100/abs(sum(density14(:,iLoc)));
	delatChange(i,4) = -(abs(sum(velocity14(:,iLoc))) - abs(sum(velocityT14(:,iLoc))))*100/abs(sum(velocity14(:,iLoc)));
	delatChange(i,5) = -(abs(sum(temperature14(:,iLoc))) - abs(sum(temperatureT14(:,iLoc))))*100/abs(sum(temperature14(:,iLoc)));
	delatChange(i,6) = -(abs(sum(pressure14(:,iLoc))) - abs(sum(pressureT14(:,iLoc))))*100/abs(sum(pressure14(:,iLoc)));
	delatChange(i,7) = -(abs(sum(KE(:,iLoc))) - abs(sum(KET(:,iLoc))))*100/abs(sum(KE(:,iLoc)));
end


% Comparison with exact results
% Volume averaged results
	% 1+(1.4-1)/2 * Me^2 = 6
	% gamma = 1.4 Simulation results
	% AvgExitMach
	Exit(1,1) = sum(mach14(:,end))/Ny;
	% MinExitMach 
	Exit(1,2)= min(mach14(:,end));
	% MaxExitMach
	Exit(1,3) = max(mach14(:,end));
	Exit(1,4) = 5.0;
	Exit(1,5) = abs(Exit(1,1)-Exit(1,4))*100/Exit(1,4);

	% AvgExittemperature
	Exit(3,1) = sum(temperature14(:,end))/Ny;
	% MinExittemperature
	Exit(3,2) = min(temperature14(:,end));
	% MaxExittemperature
	Exit(3,3) = max(temperature14(:,end));
	Exit(3,4) = 1800/6;
	Exit(3,5) = abs(Exit(3,1)-Exit(3,4))*100/Exit(3,4);

	% AvgExitpressure
	Exit(4,1) = sum(pressure14(:,end))/Ny;
	% MinExitpressure
	Exit(4,2) = min(pressure14(:,end));
	% MaxExitpressure
	Exit(4,3) = max(pressure14(:,end));
	Exit(4,4) = 5.2909e+07/(6^3.5);
	Exit(4,5) = abs(Exit(4,1)-Exit(4,4))*100/Exit(4,4);
	
	% AvgExitdensity
	Exit(2,1) = sum(density14(:,end))/Ny;
	% MinExitdensity
	Exit(2,2) = min(density14(:,end));
	% MaxExitdensity
	Exit(2,3) = max(density14(:,end));
	Exit(2,4) = Exit(4,4)/(287.14*Exit(3,4));
	Exit(2,5) = abs(Exit(2,1)-Exit(2,4))*100/Exit(2,4);


% gamma = 1.4 Exact results


