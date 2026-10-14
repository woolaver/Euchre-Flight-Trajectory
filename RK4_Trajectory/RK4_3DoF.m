function [t, state, usefulData] = RK4_3DoF(state0, dt, CL_vec_file, CD_vec_file, Sref_lift, Sref_drag)
%Same physics as ODE model but using RK4 because ODE45 sucks, also all ODE
%models are combined into this one so the whole trajectory can be ran at
%once
%Inputs:
%   state0: initial state of phase to run (V, gamma, h, x)
%   dt: time step to run
%   CL_vec_file: file to use for CL values
%   CD_vec_file: file to use for CD values
%   Sref_lift: reference area to use for lift (bottom surface)
%   Sref_drag: reference area to use for drag (maximum cross-sectional area
%   in axial direction)
%Outputs:
%   t: time vector (1xN)
%   state: state array (4xN)
%   usefulData: Vector containing angle of attack (rad), CL, CD, rho, and Mach vector (5xN)
    
    climb_end = 0;
    cruise_end = 1;
    state(:, 1) = state0;
    t(1) = 0;

    m = 120; %kg
    g = 9.8; %m/s^2

    %populate first column of usefulData
    [~, a, ~, rho] = atmoscoesa(state0(3));
        
    Mach = state(1, 1)/a;

    %assume CL and CD are 0 for the first input
    CL = 0;
    CD = 0;
    
    usefulData(:, 1) = [0; CL; CD; rho; Mach];

    i = 1;
    
    while(climb_end == 0) 

        [~, a, ~, rho] = atmoscoesa(state(3, i));
        Mach = state(1, i)/a;
    
        %for now hard set alpha constant during climb
        alpha = deg2rad(4.75);

        %alpha_vec and Mach_vec must be the same size and steps used to
        %create CL_vec_file and CD_vec_file for interpolation to work
        alpha_vec = deg2rad(linspace(-30, 30, 200));
        Mach_vec = linspace(2, 8, 200);

        CL_vec = readmatrix(CL_vec_file);
        CD_vec = readmatrix(CD_vec_file);

        CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
        CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');

        D = 1/2*rho*state(1, i)^2*Sref_drag*CD;
        L = 1/2*rho*state(1, i)^2*Sref_lift*CL;

        %{
        state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];
        %}

        %RK-4 iteration
        %k1 = dt*f(t, yn)
        %k2 = dt*f(t + dt/2, yn + k1/2)
        %k3 = dt*f(t + dt/2, yn + k2/2)
        %k4 = dt*f(t + dt, yn + k3)
        %yn+1 = yn + dt/6(k1 + 2k2 + 2k3 + k4)

        k1 = [-D/m - g*sin(state(2, i));
                 L/(m*state(1, i)) - (g/state(1, i))*cos(state(2, i));
                 state(1, i)*sin(state(2, i));
                 state(1, i)*cos(state(2, i))];

        k2 = [-D/m - g*sin(state(2, i) + k1(2)/2);
                 L/(m*(state(1, i) + k1(1)/2)) - (g/(state(1, i) + k1(1)/2))*cos((state(2, i) + k1(2)/2));
                 (state(1, i) + k1(1)/2)*sin((state(2, i) + k1(2)/2));
                 (state(1, i) + k1(1)/2)*cos((state(2, i) + k1(2)/2))];

        k3 = [-D/m - g*sin(state(2, i) + k2(2)/2);
                 L/(m*(state(1, i) + k2(1)/2)) - (g/(state(1, i) + k2(1)/2))*cos((state(2, i) + k2(2)/2));
                 (state(1, i) + k2(1)/2)*sin((state(2, i) + k2(2)/2));
                 (state(1, i) + k2(1)/2)*cos((state(2, i) + k2(2)/2))];

        k4 = [-D/m - g*sin(state(2, i) + k3(2));
                 L/(m*(state(1, i) + k3(1))) - (g/(state(1, i) + k3(1)))*cos((state(2, i) + k3(2)));
                 (state(1, i) + k3(1))*sin((state(2, i) + k3(2)));
                 (state(1, i) + k3(1))*cos((state(2, i) + k3(2)))];

        state(:, i+1) = state(:, i) + (dt./6).*(k1 + 2.*k2 + 2.*k3 + k4);

        usefulData(:, i+1) = [alpha, CL, CD, rho, Mach];

        t(i+1) = t(i)+dt;

        if(state(2, i+1) < 0)
            climb_end = 1;
            cruise_end = 0;
            disp("Climb End")
        end

        i = i+1;
    end

    i = i-1;

    LD_vec = CL_vec./CD_vec;

    %logic variable to force projectile to enter glide at optimal L/D once
    %it first loses enough speed to hold altitude
    glide = 0;
    dive = 0;

    while (cruise_end == 0)

        disp(state(3, i))

        [~, a, ~, rho] = atmoscoesa(state(3, i));
        Mach = state(1, i)/a;
        
        %attempt to hold cruise altitude as long as possible
        CL_req = 2*g*m*cos(state(2, i))/(rho*state(1, i)^2*Sref_lift);

        CL_max = interp2(alpha_vec, Mach_vec, CL_vec, deg2rad(30), Mach, 'linear');
        
        if (Mach < 3 || dive == 1)
            alpha = 0;
            dive = 1;
            glide = 0;
        %elseif ((CL_req + .05) > CL_max || glide == 1)
        elseif(Mach < 5 || glide == 1)
            %alpha = optimal_LD(Mach, LD_vec);
            alpha = deg2rad(4);
            glide = 1;
        else
            alpha = findAlphaReq(CL_req, Mach, CL_vec);
        end

        CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
        CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');

        D = 1/2*rho*state(1, i)^2*Sref_drag*CD;
        L = 1/2*rho*state(1, i)^2*Sref_lift*CL;

        %{
        state_dot = [-D/m - g*sin(state(2));
                 L/(m*state(1)) - (g/state(1))*cos(state(2));
                 state(1)*sin(state(2));
                 state(1)*cos(state(2))];
        %}

        %RK-4 iteration
        %k1 = dt*f(t, yn)
        %k2 = dt*f(t + dt/2, yn + k1/2)
        %k3 = dt*f(t + dt/2, yn + k2/2)
        %k4 = dt*f(t + dt, yn + k3)
        %yn+1 = yn + dt/6(k1 + 2k2 + 2k3 + k4)

        k1 = [-D/m - g*sin(state(2, i));
                 L/(m*state(1, i)) - (g/state(1, i))*cos(state(2, i));
                 state(1, i)*sin(state(2, i));
                 state(1, i)*cos(state(2, i))];

        k2 = [-D/m - g*sin(state(2, i) + k1(2)/2);
                 L/(m*(state(1, i) + k1(1)/2)) - (g/(state(1, i) + k1(1)/2))*cos((state(2, i) + k1(2)/2));
                 (state(1, i) + k1(1)/2)*sin((state(2, i) + k1(2)/2));
                 (state(1, i) + k1(1)/2)*cos((state(2, i) + k1(2)/2))];

        k3 = [-D/m - g*sin(state(2, i) + k2(2)/2);
                 L/(m*(state(1, i) + k2(1)/2)) - (g/(state(1, i) + k2(1)/2))*cos((state(2, i) + k2(2)/2));
                 (state(1, i) + k2(1)/2)*sin((state(2, i) + k2(2)/2));
                 (state(1, i) + k2(1)/2)*cos((state(2, i) + k2(2)/2))];

        k4 = [-D/m - g*sin(state(2, i) + k3(2));
                 L/(m*(state(1, i) + k3(1))) - (g/(state(1, i) + k3(1)))*cos((state(2, i) + k3(2)));
                 (state(1, i) + k3(1))*sin((state(2, i) + k3(2)));
                 (state(1, i) + k3(1))*cos((state(2, i) + k3(2)))];

        state(:, i+1) = state(:, i) + (dt./6).*(k1 + 2.*k2 + 2.*k3 + k4);

        usefulData(:, i+1) = [alpha, CL, CD, rho, Mach];

        t(i+1) = t(i)+dt;

        if(state(3, i+1) < 0)
            cruise_end = 1;
            disp("Cruise End")
        end

        if(isnan(state(1, i+1)))
            break
        end

        i = i+1;
    end

end

