function gamma_new = getgamma(rho,u,v,w,e)
% Gamma = getgamma(density(i,j),u(i,j),v(i,j),w(i,j),totalEnergy(i,j));
	gamma_p = 1.4;
	T_theta = 3055;
	gamma_old = 1.0; 
	gamma_new = 1.4 ; 
	gammaTolarance = 1e-5;
	R = 287.1;

	while abs(gamma_old - gamma_new)>gammaTolarance
		gamma_old = gamma_new;
		pressure = (gamma_old-1)*(e-(0.5*rho*(u^2+v^2+w^2))); 
		Temperature = pressure/(rho*R);
		TemperatureRatio = T_theta/Temperature;
		gamma_new = 1 + (gamma_p - 1)/(1 + (gamma_p - 1)*((TemperatureRatio^2)* ...
			exp(TemperatureRatio)/((exp(TemperatureRatio)-1)^2)));
	end
end