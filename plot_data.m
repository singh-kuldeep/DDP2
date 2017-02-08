clear('all');
clc ;
% grid =  csvread('grids_Bump_2D.dat') ;
grid =  csvread('grids_Nozzle_2D.csv') ;
% para =  csvread('2D_parameters_B.dat') ;

x = reshape(grid(:,1),[25,155]) ;
y = reshape(grid(:,2),[25,155]) ;
% z = reshape(grid(:,3),[26,78]) ;

% density = reshape(para(:,1),[26,78]);
% densityu= reshape(para(:,2),[26,78]);
% densityv= reshape(para(:,3),[26,78]);
% densityw= reshape(para(:,4),[26,78]);
% energy  = reshape(para(:,5),[26,78]);

% for i = 1:26
%     for j =	1:78	
%         u(i,j) = densityu(i,j)/density(i,j) ; 
%         v(i,j) = densityv(i,j)/density(i,j) ; 
%         w(i,j) = densityw(i,j)/density(i,j) ; 
%         pressure(i,j) = 0.4*(energy(i,j) - 0.5*density(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)+w(i,j)*w(i,j)));
%         temperature(i,j) = pressure(i,j) /(287.14*density(i,j)) ; 
%         velocity(i,j) = sqrt(u(i,j)*u(i,j)+v(i,j)*v(i,j)+w(i,j)*w(i,j)) ;
%         mach(i,j) = velocity(i,j) / sqrt(1.4*287.14*temperature(i,j)) ;
		
% 		% non_dimensinalisation
% 		non_dim_x(i,j) = (x(i,j) - x(1,1))/(y(26,78)-y(1,1));
%         non_dim_y(i,j) = (y(i,j) - y(1,1))/(y(26,78)-y(1,1));
% 		non_dim_z(i,j) = (z(i,j) - z(1,1))/(y(26,78)-y(1,1));
% 		non_dim_den(i,j) = density(i,j)/density(1,1);
% 		non_dim_vel(i,j) = velocity(i,j)/velocity(1,1);
% 		non_dim_press(i,j) = pressure(i,j)/pressure(1,1);
% 		non_dim_temp(i,j) = temperature(i,j)/temperature(1,1);

%     end
% end

% % plotting the non dim parameter at the surface
% figure 
%     k=24 ;
%         plot(non_dim_x(k,:),non_dim_den  (k,:),'-ro','linewidth',2); hold on ;
% 		plot(non_dim_x(k,:),non_dim_vel  (k,:),'-go','linewidth',2); hold on ;
% 		plot(non_dim_x(k,:),non_dim_press(k,:),'-bo','linewidth',2); hold on ;
% 		plot(non_dim_x(k,:),non_dim_temp (k,:),'-yo','linewidth',2); hold on ;

%         title('Nondimensional parameters v/s X/L_0')
%         xlabel('X/L_0')
%         ylabel('Nondimensional parameters') ;
%         legend('density','velocity','pressure','temperature') ;

% figure(2)
%     	plot(non_dim_x(k,:),mach         (k,:),'-o','linewidth',2); hold on ;

%         title('Nondimensional parameters v/s X/L_0')
%         xlabel('X/L_0')
%         ylabel('Mach') ;
%         legend('Mach') ;

% h = surf(non_dim_x,non_dim_y,non_dim_z,) ;
% h = surf(x,y,z,mach) ;

% h = surf(non_dim_x,non_dim_y,non_dim_z,mach) ;
% h = surf(non_dim_x,non_dim_y,non_dim_z,non_dim_den) ;
% h = surf(non_dim_x,non_dim_y,non_dim_z,non_dim_vel) ;
% h = surf(non_dim_x,non_dim_y,non_dim_z,non_dim_press) ;
% h = surf(non_dim_x,non_dim_y,non_dim_z,non_dim_temp) ;

% colormap hsv
% colorbar
% hcb=colorbar
% shading interp;

% t = title(hcb,'\bf {Mach}') ;
% set(t,'Interpreter','Latex');
% xlabel('\bf{x}'); ylabel('\bf y'); zlabel('\bf z');
% title(' \bf Mach Number (M), flow over a bump')


% t = title(hcb,'$\frac{\rho}{\rho_{\infty}}$') ;
% set(t,'Interpreter','Latex');
% xlabel('\bf{X}'); ylabel('\bf y'); zlabel('\bf z');
% title(' \bf Nondimensional density(\rho/ \rho_{\infty}), flow over a bump')

% t = title(hcb,'$\frac{V}{V_{\infty}}$') ;
% set(t,'Interpreter','Latex');
% xlabel('\bf{X}'); ylabel('\bf y'); zlabel('\bf z');
% title(' \bf Nondimensional velocity(V/V_{\infty}), flow over a bump')

% t = title(hcb,'$\frac{P}{P_{\infty}}$') ;
% set(t,'Interpreter','Latex');
% xlabel('\bf{X}'); ylabel('\bf y'); zlabel('\bf z');
% title(' \bf Nondimensional pressure(P/P_{\infty}), flow over a bump')

% t = title(hcb,'$\frac{T}{T_{\infty}}$') ;
% set(t,'Interpreter','Latex');
% xlabel('\bf{X}'); ylabel('\bf y'); zlabel('\bf z');
% title(' \bf Nondimensional temperature (T/T_{\infty}), flow over a bump')