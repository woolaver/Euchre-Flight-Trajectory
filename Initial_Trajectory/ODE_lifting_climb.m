function state_dot = ODE_lifting_climb(time, state)
    %state = [V, gamma, h, x, alpha]'
    %V = veloctiy
    %gamma = flight path angle (not specific heat ratio ik it gets
    %confusing cause specific heat ratio is used later but deal with it,
    %you can figure it out
    %h = height
    %x = distance
    %alpha = angle of attack, held constant at alpha_init

    m = 120; %kg
    S = (.2^2)*pi;
    g = 9.8; %m/s^2
    
    [~, a, ~, rho] = atmoscoesa(state(3));
    gamma = 1.4;
    Mach = state(1)/a;
    Cp_max = (2/(gamma*Mach^2))*((((gamma+1)/2)*Mach^2)^((gamma)/(gamma-1))*((gamma + 1)/(2*gamma*Mach^2 - (gamma - 1)))^(1/(gamma-1))-1);
    
    CL = Cp_max*sin(state(5))^2*cos(state(5));
    CD0 = .05;
    CD = Cp_max*sin(state(5))^3 + CD0;

    D = 1/2*rho*state(1)^2*S*CD;
    L = 1/2*rho*state(1)^2*S*CL;
    
    %alpha does not change, only reason it is a state variable is because
    %ODE fucking blows
    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2));
                 0];
end