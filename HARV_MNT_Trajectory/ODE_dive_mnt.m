function state_dot = ODE_dive_mnt(time, state)
    m = 120; %kg
    g = 9.8; %m/s^2
    
    [~, a, ~, rho] = atmoscoesa(state(3));
    Mach = state(1)/a;

    %for now assuming constant alpha = 0, will look into changing at later
    %date
    alpha = 0;

    disp(state(3))

    S_ref = (.2^2)*pi; 

    %alpha_vec and Mach_vec must be same as the ones used to create CL_vec
    %and CD_vec
    alpha_vec = deg2rad(linspace(-30, 30, 500));
    Mach_vec = linspace(2, 8, 50);

    CL_vec = readmatrix("HARV_CL_50_500.csv");
    CD_vec = readmatrix("HARV_CD_50_500.csv");

    CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
    CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');

    D = 1/2*rho*state(1)^2*S_ref*CD;
    L = 1/2*rho*state(1)^2*S_ref*CL;

    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];
end