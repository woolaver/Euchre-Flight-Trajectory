function state_dot = ODE_cruise_mnt(time, state)
    m = 120; %kg
    g = 9.8; %m/s^2
    S_ref = (.2^2)*pi; %m^2, taken as maximum diameter of projectile
    
    [~, a, ~, rho] = atmoscoesa(state(3));
    Mach = state(1)/a;
    
    if(state(3) < 27500)
        alpha_best = optimalAlphaMNT(Mach);
    else
        %{
        L_req = g*m*cos(state(2));
        CL_req = L_req*2/(rho*state(1)^2*S_ref)
        alpha_best = findGammaZero(CL_req, Mach);
        %}
        alpha_best = deg2rad(-5);
    end

    %alpha_vec and Mach_vec must be same as the ones used to create CL_vec
    %and CD_vec
    alpha_vec = deg2rad(linspace(-30, 30, 100));
    Mach_vec = linspace(2, 8, 100);

    CL_vec = readmatrix("Ice_Cream_Cone_CL.csv");
    CD_vec = readmatrix("Ice_Cream_Cone_CD.csv");

    CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha_best, Mach, 'linear');
    CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha_best, Mach, 'linear');

    D = 1/2*rho*state(1)^2*S_ref*CD;
    L = 1/2*rho*state(1)^2*S_ref*CL;

    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];

end