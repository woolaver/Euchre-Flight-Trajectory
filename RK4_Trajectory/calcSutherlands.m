function mu = calcSutherlands(T)
%Calculates Sutherlands Equation to estimate dynamic viscocity
    mu_ref = 1.789*10^(-5); %kg/m*s
    T_ref = 288; %K
    S = 110; %K

    mu = mu_ref*((T/T_ref)^(3/2) * ((T_ref + S)/(T + S)));
end