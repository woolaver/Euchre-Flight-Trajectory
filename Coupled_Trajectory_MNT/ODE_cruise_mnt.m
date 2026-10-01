function state_dot = ODE_cruise_mnt(time, state)
    m = 120; %kg
    g = 9.8; %m/s^2
    
    [~, a, ~, rho] = atmoscoesa(state(3));
    gamma = 1.4;
    Mach = state(1)/a;
    Cp_max = (2/(gamma*Mach^2))*((((gamma+1)/2)*Mach^2)^((gamma)/(gamma-1))*((gamma + 1)/(2*gamma*Mach^2 - (gamma - 1)))^(1/(gamma-1))-1);
    alpha_best = optimalAlphaMNT(Cp_max);

    disp("test")

    geometry = "HARV.STL";
    S_ref = .029; 

    [CL, CD, ~, ~, ~, ~, ~] = newtonianCLCD3D(geometry, alpha_best, Cp_max, S_ref);

    D = 1/2*rho*state(1)^2*S_ref*CD;
    L = 1/2*rho*state(1)^2*S_ref*CL;

    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];

end