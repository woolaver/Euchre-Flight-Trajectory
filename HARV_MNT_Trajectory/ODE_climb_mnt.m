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
    Mach = state(1)/a;

    %for now assuming constant alpha = 0, will look into changing at later
    %date
    alpha = deg2rad(2);

    S_ref = (.2^2)*pi; 

    %alpha_vec and Mach_vec must be same as the ones used to create CL_vec
    %and CD_vec
    alpha_vec = deg2rad(linspace(-30, 30, 500));
    Mach_vec = linspace(2, 8, 50);

    CL_vec = readmatrix("HARV_CL_50_500.csv");
    CD_vec = readmatrix("HARV_CD_50_500.csv");

    CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
    CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');

    %assuming flat plate assumption for skin friction calculation
    c = 2; %m
    mu = 1.422*10^(-5); %N*s/m^2
    
    Cf_turb = @(x) .0592./(((rho.*state(1).*x)./mu).^(1/5));
    F_turb = (2/c) * integral(Cf_turb, 0, c);

    CD_turb = F_turb*cos(alpha);
    CL_turb = F_turb*sin(alpha);

    CL = CL + CL_turb;
    CD = CD + CD_turb + .05;

    disp(state(3))

    D = 1/2*rho*state(1)^2*S_ref*CD;
    L = 1/2*rho*state(1)^2*S_ref*CL;

    %V_dot = -D/m - g*sin(state(2));
    %gamma_dot = L/(m*state(1)) - (g/state(1))*cos(state(2));
    
    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];
end