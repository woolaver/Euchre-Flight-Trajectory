%Simple unpowered projectile model to get simple range of projectile down

clear 
close all
clc

mass = 120; %kg
g = 9.8; %m/s
D = .3; %m
S_ref = pi/4*D^2; %m^2
ref_length = 2; %m

grav_Force = mass*g;

Mach = 8;
alt = 1; %m
[T, a, P, rho] = atmoscoesa(alt); %[K, m/s, Pa, kg/m^3]

alpha = 20; %degrees
total_vel = Mach*a;

dt = .001; %sec

x_accel = 0;
y_accel = 0;

x_pos_vec = 0;
y_pos_vec = alt;
time_vec = 0;

x_vel = total_vel*cosd(45);
y_vel = total_vel*sind(45);

index = 1;

while alt > 0
    %NOTE: gamma is calculated using emperical formula from NASA Glenn
    %based on temperature, the temperature used is the previous post shock
    %temperature, for the first time point a value of T = 3500K was
    %estimated by claude, blame it not me if this is wrong
    %NOTE 2: modified_newtonian treats the body as a flat plate with one
    %constant angle (as opposed to the reality of a conic structure), this
    %should hold up fine for now but will definetly need to be changed
    %later
    %[CL, CDi] = modified_newtonian(alpha, total_vel, alt);
    %NOTE: Assuming that entire body is in turbulent flow, this may not
    %hold as the projectile reaches high altitudes ( > 30km)
    %CD0 = friction_drag(total_vel, alt, ref_length);

    CL = .01;
    CD = .01;
    Lift = .5*rho*x_vel^2*S_ref*CL;
    Drag = .5*rho*x_vel^2*S_ref*CD;
    %Drag = .5*rho*total_vel^2*S_ref*CDi + .5*rho*total_vel^2*S_ref*CD0;

    x_accel = x_accel - Drag/mass*dt;
    y_accel = y_accel + (Lift/mass)*dt - g*dt;

    x_vel = x_vel + x_accel*dt;
    y_vel = y_vel + y_accel*dt;

    x_pos_vec(index + 1) = x_pos_vec(index) + x_vel*dt;
    y_pos_vec(index + 1) = y_pos_vec(index) + y_vel*dt;
    time_vec(index + 1) = time_vec(index) + dt;

    alt = y_pos_vec(index);
    total_vel = sqrt(x_vel^2 + y_vel^2);
    %alpha = tan(y_vel/x_vel);
    index = index+1;
    [T, a, P, rho] = atmoscoesa(alt); %[K, m/s, Pa, kg/m^3]
end

figure()
plot(x_pos_vec, y_pos_vec, 'LineWidth', 1, 'Color', 'Blue')
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Projectile Trajectory')

figure()
plot(time_vec, x_pos_vec, 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('X Position (m)')
title('Time')