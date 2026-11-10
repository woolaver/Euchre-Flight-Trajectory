function [t, state, usefulData] = RK4_3DoF(state0, vehicleProperties)
%Same physics as ODE model but using RK4 because ODE45 sucks, also all ODE
%models are combined into this one so the whole trajectory can be ran at
%once
%Inputs:
%   state0: initial state of phase to run (V, gamma, h, x)
%   vehicleProperties: struct containing all vehicle properties needed
%   (dt, CL_vec, CD_vec, Sref, length, height, m_skin, m_structure, thermal
%   coefficients, ...)
%Outputs:
%   t: time vector (1xN)
%   state: state array (4xN)
%   usefulData: Array containing angle of attack (rad), CL, CD, rho, Mach, Wall Temperature, Internal Temperature, and Nose Temperature (8xN)
    
    dt = vehicleProperties.dt;
    
    climb_end = 0;
    cruise_end = 1;
    state(:, 1) = state0;
    t(1) = 0;

    m = vehicleProperties.m; 
    g = 9.8; %m/s^2

    CL_vec = readmatrix(vehicleProperties.CL_vec_file);
    CD_vec = readmatrix(vehicleProperties.CD_vec_file);

    %populate first column of usefulData
    [~, a, ~, rho] = atmoscoesa(state0(3));
        
    Mach = state(1, 1)/a;

    %assume CL and CD are 0 for the first input
    CL = 0;
    CD = 0;

    Tw = vehicleProperties.Tw; %K, room temperature starting temperature
    T_int = vehicleProperties.T_int; %K, room temperature starting temperature
    T_nose = vehicleProperties.T_nose; %K, room temperature starting temperature

    Sref = vehicleProperties.Sref;
    length = vehicleProperties.length;
    height = vehicleProperties.height;
    
    usefulData(:, 1) = [0; CL; CD; rho; Mach; Tw; T_int; T_nose];

    i = 1;

    %alpha_vec and Mach_vec must be the same size and steps used to
    %create CL_vec_file and CD_vec_file for interpolation to work
    alpha_vec = deg2rad(linspace(-30, 30, 200));
    Mach_vec = linspace(1.5, 9, 200); 
    
    while(climb_end == 0)

        disp("Height: " + state(3, i) + " m")
        disp("Velocity: " + state(1, i) + " m/s")

        [~, a, ~, rho] = atmoscoesa(state(3, i));
        Mach = state(1, i)/a;
    
        %for now hard set alpha constant during climb
        if (state(3, i) < 15000)
            alpha = vehicleProperties.alpha_climb;
        else
            alpha = vehicleProperties.alpha_level;
        end

        CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
        CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');

        if(isnan(CL))
            disp("CL Search Failed")
            break;
        end

        D = 1/2*rho*state(1, i)^2*Sref*CD;
        L = 1/2*rho*state(1, i)^2*Sref*CL;

        [Fx, Fy, Q_conv_rate] = refTempMethod(length, height, state(3, i), Mach, alpha, Tw);
        Q_conv_nose_rate = fayRiddellMethod(vehicleProperties, Mach, state(3, i), T_nose);

        if(isnan(Fx) || isnan(Fy) || isnan(Q_conv_rate))
            disp("Reference Temperature Method Failed")
            break;
        end

        D = D + Fx;
        L = L + Fy;

        [Tw, T_int, T_nose] = updateTemperatures(vehicleProperties, state(3, i), Tw, T_int, T_nose, Q_conv_rate, Q_conv_nose_rate);
        
        if(isnan(Tw))
            disp("Temperature Update Failed")
            break;
        end

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

        if(isnan(state(:, :)))
            disp("State Update Failed")
            break;
        end

        usefulData(:, i+1) = [alpha, CL, CD, rho, Mach, Tw, T_int, T_nose];

        t(i+1) = t(i)+dt;

        if(state(2, i+1) < 0)
            climb_end = 1;
            cruise_end = 0;
            disp("Climb End")
        end

        i = i+1;
    end

    i = i-1;

    %LD_vec = CL_vec./CD_vec;

    %logic variable to force projectile to enter glide at optimal L/D once
    %it first loses enough speed to hold altitude
    glide = 0;
    dive = 0;

    pause(3)

    while (cruise_end == 0)

        disp("Height: " + state(3, i) + " m")
        disp("Velocity: " + state(1, i) + " m/s")

        [~, a, ~, rho] = atmoscoesa(state(3, i));
        Mach = state(1, i)/a;
        
        %attempt to hold cruise altitude as long as possible
        CL_req = 2*g*m*cos(state(2, i))/(rho*state(1, i)^2*Sref);

        CL_max = interp2(alpha_vec, Mach_vec, CL_vec, deg2rad(30), Mach, 'linear');
        
        if (state(1, i) < vehicleProperties.dive_speed || dive == 1)
            if(dive == 0)
                disp("Glide End")
                pause(3)
            end

            if(state(2, i) <= deg2rad(-85))
                %hold glide path angle
                CL_req = (2*m*g*cos(state(2, i)))/(rho*state(1, i)^2*Sref);
                if (CL_req > CL_max)
                    alpha = 0;
                else
                    alpha = findAlphaReq(CL_req, Mach, CL_vec);
                end
            else
                alpha = vehicleProperties.alpha_dive;
                dive = 1;
                glide = 0;
            end
        elseif((CL_req > CL_max) || glide == 1)
            alpha = vehicleProperties.alpha_glide;
            glide = 1;
        else
            alpha = findAlphaReq(CL_req, Mach, CL_vec);
        end

        %{
        if(Mach < 4 || dive == 1)
            alpha = vehicleProperties.alpha_dive;
            dive = 1;
            glide = 0;
        else
            glide = 1;
            alpha = vehicleProperties.alpha_glide;
        end
        %}

        CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
        CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');

        D = 1/2*rho*state(1, i)^2*Sref*CD;
        L = 1/2*rho*state(1, i)^2*Sref*CL;

        [Fx, Fy, Q_conv_rate] = refTempMethod(length, height, state(3, i), Mach, alpha, Tw);
        Q_conv_nose_rate = fayRiddellMethod(vehicleProperties, Mach, state(3, i), T_nose);

        D = D + Fx;
        L = L + Fy;

        [Tw, T_int, T_nose] = updateTemperatures(vehicleProperties, state(3, i), Tw, T_int, T_nose, Q_conv_rate, Q_conv_nose_rate);

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

        usefulData(:, i+1) = [alpha, CL, CD, rho, Mach, Tw, T_int, T_nose];

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

