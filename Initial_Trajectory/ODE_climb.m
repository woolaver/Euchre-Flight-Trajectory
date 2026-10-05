function state_dot = ODE_climb(time, state)
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
    S = (.2^2)*pi;
    g = 9.8; %m/s^2
    
    [~, a, ~, rho] = atmoscoesa(state(3));
    gamma = 1.4;
    Mach = state(1)/a;
    Cp_max = (2/(gamma*Mach^2))*((((gamma+1)/2)*Mach^2)^((gamma)/(gamma-1))*((gamma + 1)/(2*gamma*Mach^2 - (gamma - 1)))^(1/(gamma-1))-1);

    %logic for trimming projectile at max height
    %if (state(3) > 28000)
    %    state(2) = 0;
    %end
    
    %Commented CL and CD estimations from Estimated Aerodynaics of All-Body
    %Hypersonic Aircraft Configurations, Note: CD is ignoring friction
    %drag, only accounting for estimation of induced drag plus a .05 estimation of friction drag, to get friction
    %drag need to integrate over body surface, will want to do later but
    %hoping that CFD results can be gotten and validated by the point where
    %that would be necessary

    %Actual estimations are from common CL/CD model, not as accurate for
    %winged glide vehicles but easier to test with

    %h_dot = state(1)*sin(state(2));
    %[alpha] = choose_alpha(state(1), h_dot, Mach, state(3));

    %for now assuming constant alpha = 0, will look into changing at later
    %date
    alpha = 0;

    CL = Cp_max*sin(alpha)^2*cos(alpha);
    %beta = (abs(Mach^2 - 1))^(1/2);
    %C1 = 4.17/beta - .13;
    %C2 = exp(.955 - (4.35/Mach));
    %CL = C1*sin(alpha) + C2*sin(alpha)^2;

    CD0 = .05;
    CD = Cp_max*sin(alpha)^3 + CD0;

    %Km = 1;
    %CD = Km*CL*tan(alpha) + .05;

    D = 1/2*rho*state(1)^2*S*CD;
    L = 1/2*rho*state(1)^2*S*CL;
    
    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];
end