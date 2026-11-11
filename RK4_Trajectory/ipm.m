function M = ipm(nu, Mg)
% function M = ipm(nu);
% Calculates inverse of the Prandtl-Meyer function
% INPUTS:
%  nu : Prandtl-Meyer function value
%  Mg: guess for the Mach number
% OUTPUTS:
%   M : Mach number

M = fsolve(@(M) pmres(M, nu), Mg, optimset('display', 'off'));

%------------------------------------------------------
function R = pmres(M, nu)
% Prandtl-Meyer function residual
R = -nu + pm(M);
