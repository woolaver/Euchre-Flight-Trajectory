function state_dot = ODE_climb_mnt(time, state)
    %state = [V, alpha, h, x]'
    %V = veloctiy
    %gamma = flight path angle (not specific heat ratio ik it gets
    %confusing cause specific heat ratio is used later but deal with it,
    %you can figure it out
    %h = height
    %x = distance
    %Note: This is only for the climb section of flight, there will be two
    %other ODE's for cruise and for descent
    m = 120; %kg
    g = 9.8; %m/s^2
    
    [~, a, ~, rho] = atmoscoesa(state(3));
    gamma = 1.4;
    Mach = state(1)/a;
    Cp_max = (2/(gamma*Mach^2))*((((gamma+1)/2)*Mach^2)^((gamma)/(gamma-1))*((gamma + 1)/(2*gamma*Mach^2 - (gamma - 1)))^(1/(gamma-1))-1);

    %for now assuming constant alpha = 0, will look into changing at later
    %date
    alpha = 0;

    geometry = "../Surface_Methods/Bunny_Bomb - Revolve2.stl";
    S_ref = .47; 

    [CL, CD, ~, ~, ~, ~, ~] = newtonianCLCD3D(geometry, alpha, Cp_max, S_ref);

    D = 1/2*rho*state(1)^2*S_ref*CD;
    L = 1/2*rho*state(1)^2*S_ref*CL;
    
    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];
end