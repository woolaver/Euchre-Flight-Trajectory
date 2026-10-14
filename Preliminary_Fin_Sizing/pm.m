function nu = pm(M)
% function nu = pm(M)
% Calculates the Prandtl-Meyer function
% INPUTS:
%   M : Mach number
% OUTPUTS:
%  nu : Prandtl-Meyer function value

g = 1.4; % specific heat ratio

% Prandtl-Meyer function
nu = sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(M.^2-1))) - atan(sqrt(M.^2-1));

