function state_dot = ODE_cruise(time, state)
    m = 120; %kg
    S = (.2^2)*pi;
    g = 9.8; %m/s^2
    
    [~, a, ~, rho] = atmoscoesa(state(3));
    gamma = 1.4;
    Mach = state(1)/a;
    alpha_best = optimalAlpha(Mach);
    Cp_max = (2/(gamma*Mach^2))*((((gamma+1)/2)*Mach^2)^((gamma)/(gamma-1))*((gamma + 1)/(2*gamma*Mach^2 - (gamma - 1)))^(1/(gamma-1))-1);

    %{
    % Lift required to keep gamma approximately constant
    L_req = m*g*cos(state(2));

    CL_req = 2*L_req / (rho*state(1)^2*S);

    % Search allowed alpha range
    alpha_test = deg2rad(linspace(0,40,1000));
    
    CL_test = Cp_max .* ...
          sin(alpha_test).^2 .* ...
          cos(alpha_test);

    % Pick alpha giving closest CL
    [~, idx] = min(abs(CL_test - CL_req));
    
    alpha = alpha_test(idx);
    %}

    CL = Cp_max*sin(alpha_best)^2*cos(alpha_best);
    CD0 = .05;
    CD = Cp_max*sin(alpha_best)^3 + CD0;

    D = 1/2*rho*state(1)^2*S*CD;
    L = 1/2*rho*state(1)^2*S*CL;

    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];

end