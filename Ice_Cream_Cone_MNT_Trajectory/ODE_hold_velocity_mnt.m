function state_dot = ODE_hold_velocity_mnt(time, state)
    m = 120; %kg
    g = 9.8; %m/s^2
    
    [~, a, ~, rho] = atmoscoesa(state(3));
    Mach = state(1)/a;

    disp("Height: " + state(3))
    disp("Velocity: " + state(1))
   
    S_ref = (.2^2)*pi; 

    %alpha_vec and Mach_vec must be same as the ones used to create CL_vec
    %and CD_vec
    alpha_vec = deg2rad(linspace(-30, 30, 100));
    Mach_vec = linspace(2, 8, 100);

    CL_vec = readmatrix("Ice_Cream_Cone_CL.csv");
    CD_vec = readmatrix("Ice_Cream_Cone_CD.csv");

    % CD required to make V_dot = 0
    CD_req = -2*m*g*sin(state(2))/(rho*state(1)^2*S_ref);

    % Minimum achievable drag coefficient
    alpha = deg2rad(-1);
    CD_min = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');

    if(state(3) < 27500)
        alpha = deg2rad(-1);
    else
        %{
        L_req = g*m*cos(state(2));
        CL_req = L_req*2/(rho*state(1)^2*S_ref)
        alpha_best = findGammaZero(CL_req, Mach);
        %}
        alpha = deg2rad(-5);
    end

    %{
    if CD_req >= CD_min
        % Constant velocity is physically achievable
        alpha = findAlphaMNT(CD_req, Mach);
        CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
        CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');
    else
        % Cannot hold 
        % velocity without thrust
        % Go to minimum drag and allow V to decrease
        alpha = 0;
        CL = interp2(Mach_vec, alpha_vec, CL_vec, Mach, alpha, 'linear');
        CD = interp2(Mach_vec, alpha_vec, CD_vec, Mach, alpha, 'linear');
    end
    %}

    CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
    CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');

    D = 1/2*rho*state(1)^2*S_ref*CD;
    L = 1/2*rho*state(1)^2*S_ref*CL;

    state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];
end