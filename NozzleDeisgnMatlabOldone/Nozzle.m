function Dozzle(G,Me,n)

%{
    Defines geometry for the diverging section a 2-D (rectangular) supersonic 
    nozzle based on a design exit mach number for a certain gas (defined by 
    G-value), given a finite number (n) of mach waves. Based on the information
    described in Anderson, Modern Compresisble Flow 3rd Edition (Library of 
    Congress CN: 2002067852).

    Known issue - If n is too high for a given "Me", self-intersecting
    geometry will occur at the approximate location of the end of the 
    predefined arc section (xArc(n,1), yArc(n,1)).

Input parameters
    G is gamma, the ratio of specific heats (Cp/Cv)
    Me is the design exit mach number
    n is the finite number of expansion waves used in approximation

Uses Arc and PMF functions
    [x,y] = Arc(r,t) - r = radius of curvature, t = theta [deg] (0,90)
    [ M nu mu ] = PMF(G,M,nu,mu) - G = Gamma, M = Mach number, nu =
        Prandtl-Meyer angle, mu = mach angle.

%}
n = 20;
Me = 5;
G = 1.4;

%% Initialize datapoint matrices
Km = zeros(n,n);    % K- vlaues (Constant along right running characteristic lines)
Kp = zeros(n,n);    % K- vlaues (Constant along left running characteristic lines)
Theta = zeros(n,n); % Flow angles relative to the horizontal
Mu = zeros(n,n);    % Mach angles
M = zeros(n,n);     % Mach Numbers
x = zeros(n,n);     % x-coordinates
y = zeros(n,n);     % y-coordinates

%% Find NuMax (maximum expansion angle)
[~, ThetaMax, ~] = PMF(G,Me,0,0);
NuMax = ThetaMax/2;

%% Define some flow parameters of originating characteristic lines
dT = NuMax/n;
ThetaArc(:,1) = (0:dT:NuMax);

NuArc = ThetaArc;
KmArc = ThetaArc + NuArc;
[~, ~, MuArc(:,1)] = PMF(G,0,NuArc(:,1),0);

%% Coordinates of wall along curve from throat
y0 = 1; % Define throat half-height
ThroatCurveRadius = 3*y0; % Radius of curvature just downstream of the throat
[xarc, yarc] = Arc(ThroatCurveRadius,ThetaArc); % Finds x- and y-coordinates for given theta-values
yarc(:,1) = yarc(:,1) + y0; % Defines offset due to arc being above horizontal

%% Fill in missing datapoint info along first C+ line
% First centerline datapoint done manually
Km(:,1) = KmArc(2:length(KmArc),1);
Theta(:,1) = ThetaArc(2:length(KmArc),1);
Nu(:,1) = Theta(:,1);
Kp(:,1) = Theta(:,1)-Nu(:,1);
M(1,1) = 1;
Nu(1,1) = 0;
Mu(1,1) = 90;
y(1,1) = 0;
x(1,1) = xarc(2,1) + (y(1,1) - yarc(2,1))/tand((ThetaArc(2,1) - MuArc(2,1) - MuArc(2,1))/2);

% Finds the information at interior nodes along first C+ line
for i=2:n

        [M(i,1), Nu(i,1), Mu(i,1)] = PMF(G,0,Nu(i,1),0);
        s1 = tand((ThetaArc(i+1,1) - MuArc(i+1,1) + Theta(i,1) - Mu(i,1))/2);
        s2 = tand((Theta(i-1,1) + Mu(i-1,1) + Theta(i,1) + Mu(i,1))/2);
        x(i,1) = ((y(i-1,1)-x(i-1,1)*s2)-(yarc(i+1,1)-xarc(i+1,1)*s1))/(s1-s2);
        y(i,1) = y(i-1,1) + (x(i,1)-x(i-1,1))*s2;
        
end

%% Find flow properties at remaining interior nodes
for j=2:n;
    for i=1:n+1-j;
        
        Km(i,j) = Km(i+1,j-1);
        
        if i==1;
            
            Theta(i,j) = 0;
            Kp(i,j) = -Km(i,j);
            Nu(i,j) = Km(i,j);
            [M(i,j), Nu(i,j), Mu(i,j)] = PMF(G,0,Nu(i,j),0);
            s1 = tand((Theta(i+1,j-1)-Mu(i+1,j-1)+Theta(i,j)-Mu(i,j))/2);
            x(i,j) = x(i+1,j-1) - y(i+1,j-1)/s1;
            y(i,j) = 0;
            
        else
            
            Kp(i,j) = Kp(i-1,j);
            Theta(i,j) = (Km(i,j)+Kp(i,j))/2;
            Nu(i,j) = (Km(i,j)-Kp(i,j))/2;
            [M(i,j), Nu(i,j), Mu(i,j)] = PMF(G,0,Nu(i,j),0);
            s1 = tand((Theta(i+1,j-1)-Mu(i+1,j-1)+Theta(i,j)-Mu(i,j))/2);
            s2 = tand((Theta(i-1,j)+Mu(i-1,j)+Theta(i,j)+Mu(i,j))/2);
            x(i,j) = ((y(i-1,j)-x(i-1,j)*s2)-(y(i+1,j-1)-x(i+1,j-1)*s1))/(s1-s2);
            y(i,j) = y(i-1,j) + (x(i,j)-x(i-1,j))*s2;
            
        end
        
    end
end

%% Find wall node information
xwall = zeros(2*n,1);
ywall = xwall;
ThetaWall = ywall;

xwall(1:n,1) = xarc(2:length(xarc),1);
ywall(1:n,1) = yarc(2:length(xarc),1);
ThetaWall(1:n,1) = ThetaArc(2:length(xarc),1);
for i=1:n-1
    ThetaWall(n+i,1) = ThetaWall(n-i,1);
end

for i=1:n
    
    s1 = tand((ThetaWall(n+i-1,1) + ThetaWall(n+i,1))/2);
    s2 = tand(Theta(n+1-i,i)+Mu(n+1-i,i));
    xwall(n+i,1) = ((y(n+1-i,i)-x(n+1-i,i)*s2)-(ywall(n+i-1,1)-xwall(n+i-1,1)*s1))/(s1-s2);
    ywall(n+i,1) = ywall(n+i-1,1) + (xwall(n+i,1)-xwall(n+i-1,1))*s1;
        
end

%% Provide wall geometry to user
assignin('caller','xwall',xwall)
assignin('caller','ywall',ywall)
assignin('caller','Coords',[xwall ywall])

%%
% Draw contour and characteristic web
if 0
    
    plot(xwall,ywall,'-')
    axis equal
    axis([0 ceil(xwall(length(xwall),1)) 0 ceil(ywall(length(ywall),1))])
    hold on
    plot(xarc,yarc,'k-')
    
    for i=1:n-1
        plot(x(1:n+1-i,i),y(1:n+1-i,i))
    end

    for i=1:n
        plot([xarc(i,1) x(i,1)],[yarc(i,1) y(i,1)])
        plot([x(n+1-i,i) xwall(i+n,1)],[y(n+1-i,i) ywall(i+n,1)])
    end

    for c=1:n
        for r=2:n+1-c
            plot([x(c,r) x(c+1,r-1)],[y(c,r) y(c+1,r-1)])
        end
    end

    xlabel('Length [x/y0]')
    ylabel('Height [y/y0]')
    % close all;
end

    
    % Deleteing the wrong points
    for j=1:10
        for i=1:2*n-15
            if(xwall(i+1,1)<xwall(i,1))
                % delete(ywall[i+1]);
                % delete(xwall[i+1]);
                ywall(i+1) = []; 
                xwall(i+1) = [];

                ywall(i+2) = []; 
                xwall(i+2) = []; 
            end 
        end
    end
    
    if 1
        dx = xwall(2,1) - xwall(1,1);
        
        %extending the inlet with at 10 degree angle and extra p points
        for i=1:15
            xwall = [xwall(1,1)-dx; xwall];
            ywall = [ywall(1,1) + tand(10)*dx; ywall];
            dx = 1.1*dx;
        end 
        
        % % putting 3 extra points at the begining  
        % putting InletExtraPoints extra points at the begining
        InletExtraPoints = 5;
        for i=1:InletExtraPoints
        xwall = [xwall(1,1)-3*dx; xwall];
        ywall = [ywall(1,1); ywall];
    end    
        [row,col] = size(ywall); 
        
        for i=1:row-1
            xwall(i,1) = 0.5*(xwall(i,1)+xwall(i+1,1));
            ywall(i,1) = 0.5*(ywall(i,1)+ywall(i+1,1));
        end
        figure(1)
        plot(xwall,ywall,'o');
        hold on; 
        wall = [xwall,ywall];
        csvwrite('CoordinatesUpperWall.csv',wall);
        % csvwrite('XCoordinatesUpperWall.csv',xwall);
        % csvwrite('YCoordinatesUpperWall.csv',ywall);
        % plot(xwall,zeros(row,1),'o');
        % csvwrite('XCoordinatesLowerWall.csv',xwall);
        % csvwrite('YCoordinatesLowerWall.csv',zeros(row,1));
        ThroatArea = min(ywall(:,1))
        InletArea = ywall(1,1)
        InletAreaRatio =  InletArea/ThroatArea 

end
