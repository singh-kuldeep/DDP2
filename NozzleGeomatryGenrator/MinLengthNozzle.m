function MinLengthNozzle(Me)
G = 1.4 ;
% Me = 5; 
n = 150
disp('to reduce the time decrease the n ') 
%{
    Defines geometry for a minimum length nozzle based on a design exit
    mach number for a certain gas, given a finite number (n) of mach waves.
    Based on the information described in Anderson, Modern Compressible
    Flow 3rd Edition (Library of Congress CN: 2002067852).

Input parameters
    G is gamma, the ratio of specific heats (Cp/Cv)
    Me is the design exit mach number
    n is the finite number of expansion waves used in approximation
    
%}
    
%% Initialize datapoint matrices
Km = zeros(n,n);    % K- vlaues (Constant along right running characteristic lines)
Kp = zeros(n,n);    % K- vlaues (Constant along left running characteristic lines)
Theta = zeros(n,n); % Flow angles relative to the horizontal
Mu = zeros(n,n);    % Mach angles
M = zeros(n,n);     % Mach Numbers
x = zeros(n,n);     % x-coordinates
y = zeros(n,n);     % y-coordinates

%% Find NuMax (maximum angle of expansion corner)
[~, B, ~] = PMF(G,Me,0,0);
NuMax = B/2;

%% Define flow of first C+ line
y0 = 1;
x0 = 0;

dT = NuMax/n;
Theta(:,1) = (dT:dT:NuMax);

Nu = Theta;
Km = Theta + Nu;
Kp = Theta - Nu;
[M(:,1) Nu(:,1) Mu(:,1)] = PMF(G,0,Nu(:,1),0);

%% Fill in missing datapoint info along first C+ line
y(1,1) = 0;
x(1,1) = x0 - y0/tand(Theta(1,1)-Mu(1,1));
for i=2:n;
    
    s1 = tand(Theta(i,1)-Mu(i,1));
    s2 = tand((Theta(i-1,1)+Mu(i-1,1)+Theta(i,1)+Mu(i,1))/2);
    x(i,1) = ((y(i-1,1)-x(i-1,1)*s2)-(y0-x0*s1))/(s1-s2);
    y(i,1) = y(i-1) + (x(i,1)-x(i-1,1))*s2;
    
end

%% Find flow properties in characteristic web
for j=2:n;
    for i=1:1+n-j;
        
        Km(i,j) = Km(i+1,j-1);
        
        if i==1;
            
            Theta(i,j) = 0;
            Kp(i,j) = -Km(i,j);
            Nu(i,j) = Km(i,j);
            [M(i,j) Nu(i,j) Mu(i,j)] = PMF(G,0,Nu(i,j),0);
            s1 = tand((Theta(i+1,j-1)-Mu(i+1,j-1)+Theta(i,j)-Mu(i,j))/2);
            x(i,j) = x(i+1,j-1) - y(i+1,j-1)/s1;
            y(i,j) = 0;
            
        else
            
            Kp(i,j) = Kp(i-1,j);
            Theta(i,j) = (Km(i,j)+Kp(i,j))/2;
            Nu(i,j) = (Km(i,j)-Kp(i,j))/2;
            [M(i,j) Nu(i,j) Mu(i,j)] = PMF(G,0,Nu(i,j),0);
            s1 = tand((Theta(i+1,j-1)-Mu(i+1,j-1)+Theta(i,j)-Mu(i,j))/2);
            s2 = tand((Theta(i-1,j)+Mu(i-1,j)+Theta(i,j)+Mu(i,j))/2);
            x(i,j) = ((y(i-1,j)-x(i-1,j)*s2)-(y(i+1,j-1)-x(i+1,j-1)*s1))/(s1-s2);
            y(i,j) = y(i-1,j) + (x(i,j)-x(i-1,j))*s2;
            
        end
        
    end
end

%% Find wall datapoint info
xwall = zeros(1,n+1);
ywall = zeros(1,n+1);

xwall(1,1) = x0;
ywall(1,1) = y0;

walls = tand(NuMax);
webs = tand(Theta(n,1)+Mu(n,1));

xwall(1,2) = ((y(n,1)-x(n,1)*webs)-(ywall(1,1)-xwall(1,1)*walls))/(walls-webs);
ywall(1,2) = ywall(1,1)+(xwall(1,2)-xwall(1,1))*walls;

for j=3:n+1;
    
    walls = tand((Theta(n-j+3,j-2)+Theta(n-j+2,j-1))/2);
    webs = tand(Theta(n-j+2,j-1)+Mu(n-j+2,j-1));
    xwall(1,j) = ((y(n-j+2,j-1)-x(n-j+2,j-1)*webs)-(ywall(1,j-1)-xwall(1,j-1)*walls))/(walls-webs);
    ywall(1,j) = ywall(1,j-1) + (xwall(1,j)-xwall(1,j-1))*walls;
    
end

%% Provide wall geometry to user and plot
assignin('base','xwall',xwall)
assignin('base','ywall',ywall)

grid=0;
if grid == 1
    
    plot(xwall,ywall,'-')
    axis equal
    axis([0 ceil(xwall(1,length(xwall))) 0 ceil(ywall(1,length(ywall)))])
    hold on
    
    for i=1:n
        plot([0 x(i,1)],[1 y(i,1)])
        plot([x(n+1-i,i) xwall(1,i+1)],[y(n+1-i,i) ywall(1,i+1)])
    end
    
    for i=1:n-1
        plot(x(1:n+1-i,i),y(1:n+1-i,i))
    end

    for c=1:n
        for r=2:n+1-c
            plot([x(c,r) x(c+1,r-1)],[y(c,r) y(c+1,r-1)])
        end
    end

xlabel('Length [x/y0]')
ylabel('Height [y/y0]')

end

% Modification in the nozzle
ThroatArea = ywall(1,1);
dx = xwall(1,2) - xwall(1,1);

% %extending the inlet with at 10 degree angle and extra p points
for i=1:7
    xwall = [xwall(1,1)-dx xwall];
    ywall = [ywall(1,1) + tand(5)*dx ywall];
    dx = 1.02*dx;
end 

% % putting ExtraPoints extra points at the begining of the inlet
ExtraPoints = 3;
for i=1:ExtraPoints
xwall = [xwall(1,1)-1.5*dx xwall];
ywall = [ywall(1,1) ywall];
end    

h = figure(1)
set(gca,'fontsize',18)
hold on; 
% plot(xwall,ywall,'-','Linewidth',2);
plot(xwall,ywall,'b--o');
hold on; 

% % putting 10 extra points at the end of the exit
for i=1:10
xwall = [xwall xwall(1,end)+(xwall(1,end)-xwall(1,end-1))];
ywall = [ywall ywall(1,end)];
% plot(xwall(1,end),ywall(1,end),'.','MarkerEdgeColor','r');
plot(xwall(1,end),ywall(1,end),'o','MarkerEdgeColor','r');
end
title('Nozzle upper wall coordinates calculated using MOC')
xlabel('x(m)'); 
ylabel('y(m)');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h,'./NozzleUpperwall','epsc')

assignin('base','xwall',xwall)
assignin('base','ywall',ywall)

upperwall(:,1) = xwall; 
upperwall(:,2) = ywall;

assignin('base','upperwall',upperwall)

csvwrite('CoordinatesUpperWall.csv',upperwall);

% [row,col] = size(ywall) 
% csvwrite('XCoordinatesUpperWall.csv',xwall);
% csvwrite('YCoordinatesUpperWall.csv',ywall);

% % plot(xwall,zeros(row,1),'o');
% csvwrite('XCoordinatesLowerWall.csv',xwall);
% csvwrite('YCoordinatesLowerWall.csv',zeros(row,1));

InletArea = ywall(1,1);
InletAreaRatio =  InletArea/ThroatArea
getMach(InletAreaRatio);


end
% end