function state_dot = ODE_hold_velocity(time, state)
    m = 120; %kg
    S = (.2^2)*pi;
    g = 9.8; %m/s^2
    
    [~, a, ~, rho] = atmoscoesa(state(3));
    gamma = 1.4;
    Mach = state(1)/a;
    Cp_max = (2/(gamma*Mach^2))*((((gamma+1)/2)*Mach^2)^((gamma)/(gamma-1))*((gamma + 1)/(2*gamma*Mach^2 - (gamma - 1)))^(1/(gamma-1))-1);
    CD0 = .05;

% CD required to make V_dot = 0
CD_req = -2*m*g*sin(state(2)) / ...
         (rho*state(1)^2*S);

% Minimum achievable drag coefficient
CD_min = CD0;

if CD_req >= CD_min
    % Constant velocity is physically achievable
    alpha = asin(((CD_req - CD0)/Cp_max)^(1/3));

else
    % Cannot hold velocity without thrust
    % Go to minimum drag and allow V to decrease
    alpha = 0;
end

CL = Cp_max*sin(alpha)^2*cos(alpha);
CD = Cp_max*sin(alpha)^3 + CD0;

    D = 1/2*rho*state(1)^2*S*CD;
    L = 1/2*rho*state(1)^2*S*CL;

    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];
end