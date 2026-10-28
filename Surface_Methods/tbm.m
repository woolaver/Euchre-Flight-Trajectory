function cp = tbm(M1, t)
% function [p2op1, M2] = tbm(M1, t)
% Calculates post-oblique shock properties
% INPUTS:
%  M1 : inflow Mach number
%  t  : theta (turning angle)
% OUTPUTS:
%  p2op1 : (post-shock pressure)/(inflow pressure)
%  M2    : post-shock Mach number

% solve for beta
b = fsolve(@(b) rel(t, b, M1), 45*pi/180, optimset('display', 'off'));

% g = gamma = specific heat ratio
g = 1.4;
rho=1.224;
a=343;

% calculate p2/p1
Mn1 = M1*sin(b);
Mn2 = sqrt((1+(g-1)/2*Mn1^2)/(g*Mn1^2-(g-1)/2));
M2 = Mn2/sin(b - t);


p2op1 = 1+2*g/(g+1)*(Mn1^2-1);
t2op1 = ((2.*g.*Mn1.^2)./(g+1)-(g-1)./(g+1))  .*  ((g-1)./(g+1)+   2./((g+1).*Mn1.^2));
rho2op1 = (g+1).*Mn1.^2./(2+(g-1).*Mn1.^2);

v = M1*a;
cp = 2/(g*M1^2) * (p2op1 - 1);

%---------------------------------------------
function R = rel(t,b,M1)
% residual of theta-beta-M relationship
g = 1.4;
R = tan(t) - 2*cot(b).*(M1^2*(sin(b)).^2 - 1)./(M1^2*(g+cos(2*b))+2);
