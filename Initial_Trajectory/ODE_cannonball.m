function state_dot = ODE_cannonball(time, state)
    m = 120; %kg
    S = (.2^2)*pi;
    g = 9.8; %m/s^2
    
    [~, ~, ~, rho] = atmoscoesa(state(3));

    CL = 0;
    CD = .05;

    D = 1/2*rho*state(1)^2*S*CD;
    L = 1/2*rho*state(1)^2*S*CL;
    
    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];
end