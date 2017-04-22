function [ x, y ] = Arc(r,t)

% Finds x and y locations of slope (t) along a semi-circular arc
% r = radius of curvature
% t = tangent angle
% x = x-coordinate
% y = y-coordinate

% y = r - sqrt(r^2 - x^2)   Equation of semi-circle 
% dy/dx = tand(t) = x/sqrt(r^2 - x^2)
% solving for x: x = r*tand(t)/sqrt(1 + tand(t)^2)

x = r.*tand(t)./sqrt(1 + tand(t).^2);

y = r - sqrt(r.^2 - x.^2);

end