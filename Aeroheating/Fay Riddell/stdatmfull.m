function [rho,p,T] = stdatmfull(h)
% PURPOSE: Computes density (rho [kg/m^3]) at a specified altitude (h [m])
%          Also returns pressure (Pa) and temperature (K)

g  = 9.8;  % gravitational acceleration [m/s^2]
R  = 287;  % gas constant for air [J/kg.K]

hl = 1000*[0, 11, 20, 32, 47, 51, 71, 85, 90, 140, 500];  % altitudes at layer endpoints [m]
Tl = [288, 217, 217, 229, 271, 271, 215, 187, 187, 320, 700];
Ll = (Tl(2:end)-Tl(1:end-1))./(hl(2:end)-hl(1:end-1)); % lapse rate for each layer [K/m]
pl = [101325, 2.26656e+04, 5.49942e+03, 8.75227e+02, 1.12265e+02, ...
      6.78201e+01, 4.03063e+00, 3.72214e-01, 1.49376e-01, 1.51115e-04];% pressures at endpoints [Pa]

% identify the layer in which h is found
h = max(h,0); I = find(hl-h>0);
if (isempty(I)), rho = 0.; p=0; T=288; return; end;
l = I(1)-1; T = Tl(l) + Ll(l)*(h-hl(l)); % temperature at the point
p0 = pl(l);                  % pressure at the start of the layer

% pressure, density at the point
if (Ll(l) == 0), % isothermal layer
  p = p0*exp(-g*(h-hl(l))./(R*T));
else             % constant-temperature-gradient layer
  p = p0*(T/Tl(l)).^(-g/(R*Ll(l)));
end
rho = p/(R*T);

