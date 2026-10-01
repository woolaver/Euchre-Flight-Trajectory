function state_dot = ODE_hold_velocity_mnt(time, state)
    m = 120; %kg
    g = 9.8; %m/s^2
    
    [~, a, ~, rho] = atmoscoesa(state(3));
    gamma = 1.4;
    Mach = state(1)/a;
    Cp_max = (2/(gamma*Mach^2))*((((gamma+1)/2)*Mach^2)^((gamma)/(gamma-1))*((gamma + 1)/(2*gamma*Mach^2 - (gamma - 1)))^(1/(gamma-1))-1);

    geometry = "HARV.STL";
    S_ref = .029; 

    % CD required to make V_dot = 0
    CD_req = -2*m*g*sin(state(2))/(rho*state(1)^2*S_ref);

    % Minimum achievable drag coefficient
    alpha = 0;
    [~, CD_min, ~, ~, ~, ~, ~] = newtonianCLCD3D(geometry, alpha, Cp_max, S_ref);


if CD_req >= CD_min
    % Constant velocity is physically achievable
    alpha = findAlphaMNT(CD_req, Cp_max);
    [CL, CD, ~, ~, ~, ~, ~] = newtonianCLCD3D(geometry, alpha, Cp_max, S_ref);
else
    % Cannot hold 
    % velocity without thrust
    % Go to minimum drag and allow V to decrease
    alpha = 0;
    [CL, CD, ~, ~, ~, ~, ~] = newtonianCLCD3D(geometry, alpha, Cp_max, S_ref);
end

    D = 1/2*rho*state(1)^2*S_ref*CD;
    L = 1/2*rho*state(1)^2*S_ref*CL;

    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];
end