function [t, state, usefulData] = RK4_6DoF(vehicleProperties, state0, dt)
%RK4 based 6DoF model, for now roll and yaw initial states will be set to 0
%and no side forces will be applied but once fin sizing and stability is
%made we will want to have this to run Monte Carlo simulations with
%Inputs:
%   vehicleProperties: struct containing any necessary vehicle properties
%       (CL_vec, CD_vec, Sref_lift, Sref_drag, mass, Ixx, Iyy, Izz, Ixy,
%       Iyz, Ixz)
%   state0: initial state of phase to run (x, y, z, phi, theta, psi, U, V, W, P, Q, R)
%   dt: time step to run
%   
%Outputs:
%   t: time vector (1xN)
%   state: state array (4xN)
%   usefulData: Vector containing: angle of attack (rad), CL, CD, rho, and Mach vector (5xN)

%Dynamics from AE341
% x is out the noze of the vehicle, y is out the right side, z is down
% state_dot = [(cos(theta)*cos(psi)*U + (sin(phi)*sin(theta)*cos(psi))*V + (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*W;
%              (cos(theta)*sin(psi)*U + (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi))*V + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*W;
%              -(sin(theta))*U + (sin(phi)*cos(theta))*V + (cos(phi)*cos(theta))*W;
%              P + (sin(phi)*tan(theta)*Q + (cos(phi)*tan(theta))*R;
%              (cos(phi))*Q - (sin(pih))*R;
%              (sin(phi)*sec(theta))*Q + (cos(phi)*sec(theta))*R;
%              R*V - Q*W + Fx/m;
%              -R*U + P*W + Fy/m;
%              Q*U - P*V + Fz/m;
%              1/(Ixx*Izz - Ixz^2) * ((Iyy*Izz - Izz^2 - Ixz^2)*Q*R + Ixz*(Ixx - Iyy + Izz)*P*Q + Izz*L + Ixz*N);
%              1/Iyy * ((Izz - Ixx)*P*R + Ixz*(R^2 - P^2) + M);
%              1/(Ixx*Izz - Ixz^2) * ((-Ixx*Iyy + Ixx^2 + Ixz^2)*P*Q + Ixz*(-Ixx + Iyy - Izz)*Q*R + Ixz*L + Ixx*N);

    %pull out vehicle properties and initalize variables
    g = 9.8; %m/s^2
    mass = vehicleProperties.mass;

    alpha_climb = vehicleProperties.alpha_climb;
    Sref_lift = vehicleProperties.Sref_lift;
    Sref_drag = vehicleProperties.Sref_drag;

    width = vehicleProperties.width;
    length = vehicleProperties.length;

    %initalize vectors for storage
    state(:, 1) = state0;
    t(1) = 0;
    
    %call in coefficient matricies 
    % (will need to get Ce, Cl, Cm, Cn for full trajectory later)
    CL_vec = readmatrix(vehicleProperties.CL_vec_file);
    CD_vec = readmatrix(vehicleProperties.CD_vec_file);

    %alpha_vec and Mach_vec must be the same size and steps used to
    %create CL_vec_file and CD_vec_file for interpolation to work
    alpha_vec = deg2rad(linspace(-30, 30, 200));
    Mach_vec = linspace(2, 8, 200);

    %populate first column of usefulData
    [~, a, ~, rho] = atmoscoesa(state0(3));

    speed = sqrt(state0(7)^2 + state0(8)^2 + state0(9)^2);
        
    Mach = speed/a;

    %assume CL and CD are 0 for the first input
    CL = 0;
    CD = 0;
    
    usefulData(:, 1) = [0; CL; CD; rho; Mach];

    i = 1;

    %start in climb phase
    phase = "climb";

    while (phase == "climb")
        disp(state(3, i))

        [~, a, ~, rho] = atmoscoesa(state(3, i));

        speed = sqrt(state(7, i)^2 + state(8, i)^2 + state(9, i)^2);
        Mach = speed/a;
    
        %for now hard set alpha constant during climb
        %4.75 for virginia, .5 for Ice Cream Cone, 0 for HARV, 3 for Bunny Bomb 

        alpha = alpha_climb;

        %set beta to 0 for 3DoF check
        beta = 0;

        CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
        CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');
        CE = 0;

        Cl_mom = 0;
        CM_mom = 0;
        CN_mom = 0;

        if(isnan(CL))
            CL = 0;
            disp("CL returned Nan, default value of 0 is used")
        end

        if(isnan(CD))
            CD = .1;
            disp("CD returned NaN, default value of .1 is used")
        end

        D = 1/2*rho*speed^2*Sref_drag*CD;
        E = 1/2*rho*speed^2*Sref_side*CE; %will be 0 for 3DoF
        L = 1/2*rho*speed^2*Sref_lift*CL;

        %Force Calculations
        Fx = -D*cos(alpha)*cos(beta) + E*cos(alpha)*sin(beta) + L*sin(alpha) - m*g*sin(state(5, i));
        Fy = -D*sin(beta) - E*cos(beta) + m*g*sin(state(4, i))*cos(state(5, i));
        Fz = -D*sin(alpha)*cos(beta) + E*sin(alpha)*sin(beta) - L*cos(alpha) + m*g*cos(state(4, i))*cos(state(5, i));
        
        %Moment Calculations
        l_mom = 1/2*rho*speed^2*Sref_drag*width*Cl_mom;
        M_mom = 1/2*rho*speed^2*Sref_lift*length*CM_mom;
        N_mom = 1/2*rho*speed^2*Sref_side*width*CN_mom;

        state(i+1, :) = RK4_iteration(state(i, 0), dt, vehicleProperties, Fx, Fy, Fz, l_mom, M_mom, N_mom);

        usefulData(:, i+1) = [alpha, CL, CD, rho, Mach];

        t(i+1) = t(i)+dt;

        if(state(7, i+1) < 0)
            phase = "hold_velocity";
            disp("Climb End")
        end

        i = i+1;
    end
%              Q*U - P*V + Fz/m;

    while (phase == "hold_velocity")
        disp(state(3, i))

        [~, a, ~, rho] = atmoscoesa(state(3, i));

        speed = sqrt(state(7, i)^2 + state(8, i)^2 + state(9, i)^2);
        Mach = speed/a;

        Fz_req = (-state(11, i)*state(7, i) + state(10, i)*state(8, i))/m;

        %assuming beta is 0 for alphaReq
        

        CL_max = interp2(alpha_vec, Mach_vec, CL_vec, deg2rad(30), Mach, 'linear');

        alpha = findAlphaReq(CL_req, Mach, CL_vec);

        %set beta to 0 for 3DoF check
        beta = 0;

        CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
        CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');
        CE = 0;

        Cl_mom = 0;
        CM_mom = 0;
        CN_mom = 0;

        if(isnan(CL))
            CL = 0;
            disp("CL returned Nan, default value of 0 is used")
        end

        if(isnan(CD))
            CD = .1;
            disp("CD returned NaN, default value of .1 is used")
        end

        D = 1/2*rho*speed^2*Sref_drag*CD;
        E = 1/2*rho*speed^2*Sref_side*CE; %will be 0 for 3DoF
        L = 1/2*rho*speed^2*Sref_lift*CL;

        %Force Calculations
        Fx = -D*cos(alpha)*cos(beta) + E*cos(alpha)*sin(beta) + L*sin(alpha);
        Fy = -D*sin(beta) - E*cos(beta);
        Fz = -D*sin(alpha)*cos(beta) + E*sin(alpha)*sin(beta) - L*cos(alpha) + m*g;
        
        %Moment Calculations
        l_mom = 1/2*rho*speed^2*Sref_drag*width*Cl_mom;
        M_mom = 1/2*rho*speed^2*Sref_lift*length*CM_mom;
        N_mom = 1/2*rho*speed^2*Sref_side*width*CN_mom;

        state(i+1, :) = RK4_iteration(state(i, 0), dt, vehicleProperties, Fx, Fy, Fz, l_mom, M_mom, N_mom);

        usefulData(:, i+1) = [alpha, CL, CD, rho, Mach];

        t(i+1) = t(i)+dt;

        if(state(7, i+1) < 0)
            phase = "hold_velocity";
            disp("Climb End")
        end

        i = i+1;
    end
    
    
    while (phase == "cruise")
    
    end
    
    while (phase == "dive")
    
    end
    
    end

function next_state = RK4_iteration(state, dt, vehicleProperties, Fx, Fy, Fz, L, M, N)
    phi = state(4);
    theta = state(5);
    psi = state(6);
    U = state(7);
    V = state(8);
    W = state(9);
    P = state(10);
    Q = state(11);
    R = state(12);

    m = vehicleProperties.mass;
    Ixx = vehicleProperties.Ixx;
    Iyy = vehicleProperties.Iyy;
    Izz = vehicleProperties.Izz;
    Ixz = vehicleProperties.Ixz;

    k1 = [cos(theta)*cos(psi)*U + (sin(phi)*sin(theta)*cos(psi))*V + (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*W;
          cos(theta)*sin(psi)*U + (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi))*V + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*W;
          -(sin(theta))*U + (sin(phi)*cos(theta))*V + (cos(phi)*cos(theta))*W;
          P + (sin(phi)*tan(theta))*Q + (cos(phi)*tan(theta))*R;
          (cos(phi))*Q - (sin(phi))*R;
          (sin(phi)*sec(theta))*Q + (cos(phi)*sec(theta))*R;
          R*V - Q*W + Fx/m;
          -R*U + P*W + Fy/m;
          Q*U - P*V + Fz/m;
          1/(Ixx*Izz - Ixz^2) * ((Iyy*Izz - Izz^2 - Ixz^2)*Q*R + Ixz*(Ixx - Iyy + Izz)*P*Q + Izz*L + Ixz*N);
          1/Iyy * ((Izz - Ixx)*P*R + Ixz*(R^2 - P^2) + M);
          1/(Ixx*Izz - Ixz^2) * ((-Ixx*Iyy + Ixx^2 + Ixz^2)*P*Q + Ixz*(-Ixx + Iyy - Izz)*Q*R + Ixz*L + Ixx*N)];

    k2 = [cos(theta + k1(5)/2)*cos(psi + k1(6)/2)*(U + k1(7)/2) + (sin(phi + k1(4)/2)*sin(theta + k1(5)/2)*cos(psi + k1(6)/2))*(V + k1(8)/2) + (cos(phi + k1(4)/2)*sin(theta + k1(5)/2)*cos(psi + k1(6)/2) + sin(phi + k1(4)/2)*sin(psi + k1(6)/2))*(W + k1(9)/2);
          cos(theta + k1(5)/2)*sin(psi + k1(6)/2)*(U + k1(7)/2) + (sin(phi + k1(4)/2)*sin(theta + k1(5)/2)*sin(psi + k1(6)/2) + cos(phi + k1(4)/2)*cos(psi + k1(6)/2))*(V + k1(8)/2) + (cos(phi + k1(4)/2)*sin(theta + k1(5)/2)*sin(psi + k1(6)/2) - sin(phi + k1(4)/2)*cos(psi + k1(6)/2))*(W + k1(9)/2);
          -(sin(theta + k1(5)/2))*(U + k1(7)/2) + (sin(phi + k1(4)/2)*cos(theta + k1(5)/2))*(V + k1(8)/2) + (cos(phi + k1(4)/2)*cos(theta + k1(5)/2))*(W + k1(9)/2);
          (P + k1(10)/2) + (sin(phi + k1(4)/2)*tan(theta + k1(5)/2))*(Q + k1(11)/2) + (cos(phi + k1(4)/2)*tan(theta + k1(5)/2))*(R + k1(12)/2);
          (cos(phi + k1(4)/2))*(Q + k1(11)/2) - (sin(phi + k1(4)/2))*(R + k1(12)/2);
          (sin(phi + k1(4)/2)*sec(theta + k1(5)/2))*(Q + k1(11)/2) + (cos(phi + k1(4)/2)*sec(theta + k1(5)/2))*(R + k1(12)/2);
          (R + k1(12)/2)*(V + k1(8)/2) - (Q + k1(11)/2)*(W + k1(9)/2) + Fx/m;
          -(R + k1(12)/2)*(U + k1(7)/2) + (P + k1(10)/2)*(W + k1(9)/2) + Fy/m;
          (Q + k1(11)/2)*(U + k1(7)/2) - (P + k1(10)/2)*(V + k1(8)/2) + Fz/m;
          1/(Ixx*Izz - Ixz^2) * ((Iyy*Izz - Izz^2 - Ixz^2)*(Q + k1(11)/2)*(R + k1(12)/2) + Ixz*(Ixx - Iyy + Izz)*(P + k1(10)/2)*(Q + k1(11)/2) + Izz*L + Ixz*N);
          1/Iyy * ((Izz - Ixx)*(P + k1(10)/2)*(R + k1(12)/2) + Ixz*((R + k1(12)/2)^2 - (P + k1(10)/2)^2) + M);
          1/(Ixx*Izz - Ixz^2) * ((-Ixx*Iyy + Ixx^2 + Ixz^2)*(P + k1(10)/2)*(Q + k1(11)/2) + Ixz*(-Ixx + Iyy - Izz)*(Q + k1(11)/2)*(R + k1(12)/2) + Ixz*L + Ixx*N)];

     k3 = [cos(theta + k2(5)/2)*cos(psi + k2(6)/2)*(U + k2(7)/2) + (sin(phi + k2(4)/2)*sin(theta + k2(5)/2)*cos(psi + k2(6)/2))*(V + k2(8)/2) + (cos(phi + k2(4)/2)*sin(theta + k2(5)/2)*cos(psi + k2(6)/2) + sin(phi + k2(4)/2)*sin(psi + k2(6)/2))*(W + k2(9)/2);
          cos(theta + k2(5)/2)*sin(psi + k2(6)/2)*(U + k2(7)/2) + (sin(phi + k2(4)/2)*sin(theta + k2(5)/2)*sin(psi + k2(6)/2) + cos(phi + k2(4)/2)*cos(psi + k2(6)/2))*(V + k2(8)/2) + (cos(phi + k2(4)/2)*sin(theta + k2(5)/2)*sin(psi + k2(6)/2) - sin(phi + k2(4)/2)*cos(psi + k2(6)/2))*(W + k2(9)/2);
          -(sin(theta + k2(5)/2))*(U + k2(7)/2) + (sin(phi + k2(4)/2)*cos(theta + k2(5)/2))*(V + k2(8)/2) + (cos(phi + k2(4)/2)*cos(theta + k2(5)/2))*(W + k2(9)/2);
          (P + k2(10)/2) + (sin(phi + k2(4)/2)*tan(theta + k2(5)/2))*(Q + k2(11)/2) + (cos(phi + k2(4)/2)*tan(theta + k2(5)/2))*(R + k2(12)/2);
          (cos(phi + k2(4)/2))*(Q + k2(11)/2) - (sin(phi + k2(4)/2))*(R + k2(12)/2);
          (sin(phi + k2(4)/2)*sec(theta + k2(5)/2))*(Q + k2(11)/2) + (cos(phi + k2(4)/2)*sec(theta + k2(5)/2))*(R + k2(12)/2);
          (R + k2(12)/2)*(V + k2(8)/2) - (Q + k2(11)/2)*(W + k2(9)/2) + Fx/m;
          -(R + k2(12)/2)*(U + k2(7)/2) + (P + k2(10)/2)*(W + k2(9)/2) + Fy/m;
          (Q + k2(11)/2)*(U + k2(7)/2) - (P + k2(10)/2)*(V + k2(8)/2) + Fz/m;
          1/(Ixx*Izz - Ixz^2) * ((Iyy*Izz - Izz^2 - Ixz^2)*(Q + k2(11)/2)*(R + k2(12)/2) + Ixz*(Ixx - Iyy + Izz)*(P + k2(10)/2)*(Q + k2(11)/2) + Izz*L + Ixz*N);
          1/Iyy * ((Izz - Ixx)*(P + k2(10)/2)*(R + k2(12)/2) + Ixz*((R + k2(12)/2)^2 - (P + k2(10)/2)^2) + M);
          1/(Ixx*Izz - Ixz^2) * ((-Ixx*Iyy + Ixx^2 + Ixz^2)*(P + k2(10)/2)*(Q + k2(11)/2) + Ixz*(-Ixx + Iyy - Izz)*(Q + k2(11)/2)*(R + k2(12)/2) + Ixz*L + Ixx*N)];

     k4 = [cos(theta + k3(5))*cos(psi + k3(6))*(U + k3(7)) + (sin(phi + k3(4))*sin(theta + k3(5))*cos(psi + k3(6)))*(V + k3(8)) + (cos(phi + k3(4))*sin(theta + k3(5))*cos(psi + k3(6)) + sin(phi + k3(4))*sin(psi + k3(6)))*(W + k3(9));
          cos(theta + k3(5))*sin(psi + k3(6))*(U + k3(7)) + (sin(phi + k3(4))*sin(theta + k3(5))*sin(psi + k3(6)) + cos(phi + k3(4))*cos(psi + k3(6)))*(V + k3(8)) + (cos(phi + k3(4))*sin(theta + k3(5))*sin(psi + k3(6)) - sin(phi + k3(4))*cos(psi + k3(6)))*(W + k3(9));
          -(sin(theta + k3(5)))*(U + k3(7)) + (sin(phi + k3(4))*cos(theta + k3(5)))*(V + k3(8)) + (cos(phi + k3(4))*cos(theta + k3(5)))*(W + k3(9));
          (P + k3(10)) + (sin(phi + k3(4))*tan(theta + k3(5)))*(Q + k3(11)) + (cos(phi + k3(4))*tan(theta + k3(5)))*(R + k3(12));
          (cos(phi + k3(4)))*(Q + k3(11)) - (sin(phi + k3(4)))*(R + k3(12));
          (sin(phi + k3(4))*sec(theta + k3(5)))*(Q + k3(11)) + (cos(phi + k3(4))*sec(theta + k3(5)))*(R + k3(12));
          (R + k3(12))*(V + k3(8)) - (Q + k3(11))*(W + k3(9)) + Fx/m;
          -(R + k3(12))*(U + k3(7)) + (P + k3(10))*(W + k3(9)) + Fy/m;
          (Q + k3(11))*(U + k3(7)) - (P + k3(10))*(V + k3(8)) + Fz/m;
          1/(Ixx*Izz - Ixz^2) * ((Iyy*Izz - Izz^2 - Ixz^2)*(Q + k3(11))*(R + k3(12)) + Ixz*(Ixx - Iyy + Izz)*(P + k3(10))*(Q + k3(11)) + Izz*L + Ixz*N);
          1/Iyy * ((Izz - Ixx)*(P + k3(10))*(R + k3(12)) + Ixz*((R + k3(12))^2 - (P + k3(10))^2) + M);
          1/(Ixx*Izz - Ixz^2) * ((-Ixx*Iyy + Ixx^2 + Ixz^2)*(P + k3(10))*(Q + k3(11)) + Ixz*(-Ixx + Iyy - Izz)*(Q + k3(11))*(R + k3(12)) + Ixz*L + Ixx*N)];

     next_state = state + (dt/6).*(k1 + k2.*2 + k3.*2 + k4);
end

function alpha = findAlphaReq6DoF(vehicleProperties, FzReq)