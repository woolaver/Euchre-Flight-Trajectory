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
%   usefulData: Vector containing angle of attack (rad), CL, CD, rho, and Mach vector (5xN)

%Dynamics from AE341
% state_dot = [(cos(theta)*cos(psi)*U + (sin(phi)*sin(theta)*cos(psi))*V + (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*W;
%              (cos(theta)*sin(psi)*U + (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi))*V + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*W;
%              -(sin(theta))*U + (sin(phi)*cos(theta))*V + (cos(phi)*cos(theta))*W;
%              P + (sin(phi)*tan(theta)*Q + (cos(phi)*tan(theta))*R;
%              (cos(phi))*Q - (sin(pih))*R;
%              (sin(phi)*sec(theta))*Q + (cos(phi)*sec(theta))*R;
%              R*V - Q*W + Fx/m;
%              -R*U + P*W + Fy/m;
%              Q*U - P*V + Fx/m;
%              1/(Ixx*Izz - Ixz^2) * ((Iyy*Izz - Izz^2 - Ixz^2)*Q*R + Ixz*(Ixx - Iyy + Izz)*P*Q + Izz*L + Ixz*N);
%              1/Iyy * ((Izz - Ixx)*P*R + Ixz*(R^2 - P^2) + M);
%              1/(Ixx*Izz - Ixz^2) * ((-Ixx*Iyy + Ixx^2 + Ixz^2)*P*Q + Ixz*(-Ixx + Iyy - Izz)*Q*R + Ixz*L + Ixx*N);



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
          Q*U - P*V + Fx/m;
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
          (Q + k1(11)/2)*(U + k1(7)/2) - (P + k1(10)/2)*(V + k1(8)/2) + Fx/m;
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
          (Q + k2(11)/2)*(U + k2(7)/2) - (P + k2(10)/2)*(V + k2(8)/2) + Fx/m;
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
          (Q + k3(11))*(U + k3(7)) - (P + k3(10))*(V + k3(8)) + Fx/m;
          1/(Ixx*Izz - Ixz^2) * ((Iyy*Izz - Izz^2 - Ixz^2)*(Q + k3(11))*(R + k3(12)) + Ixz*(Ixx - Iyy + Izz)*(P + k3(10))*(Q + k3(11)) + Izz*L + Ixz*N);
          1/Iyy * ((Izz - Ixx)*(P + k3(10))*(R + k3(12)) + Ixz*((R + k3(12))^2 - (P + k3(10))^2) + M);
          1/(Ixx*Izz - Ixz^2) * ((-Ixx*Iyy + Ixx^2 + Ixz^2)*(P + k3(10))*(Q + k3(11)) + Ixz*(-Ixx + Iyy - Izz)*(Q + k3(11))*(R + k3(12)) + Ixz*L + Ixx*N)];

     next_state = state + (dt/6).*(k1 + k2.*2 + k3.*2 + k4);
end