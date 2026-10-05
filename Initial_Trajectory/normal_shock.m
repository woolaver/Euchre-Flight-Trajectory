function P2 = normal_shock(Mach, gamma, P1)
%Simple normal shock code
    P2 = (2*gamma*Mach^2 - (gamma - 1))/(gamma+1)*P1; 
end