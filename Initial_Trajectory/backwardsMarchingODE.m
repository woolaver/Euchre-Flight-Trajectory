%Testing script for backwards marching from a Mach 3 impact to get an idea
%of height required for straight line dive or for constant Mach 3 flight

clear
close all
clc

Mach_end = 3;
%straight dive
launch_angle = deg2rad(0);
ground = 0;
[~, a, ~, ~] = atmoscoesa(ground); %m/s
V0 = Mach_end*a;

state_0 = [V0; launch_angle; ground; 0];
%integrate for 30 seconds
time_range = [0 -10];

%ode_hold_velocity to see when to begin velocity hold, alpha is determined
%by finding drag value that corresponds to V_dot = 0
[t, state] = ode45(@ODE_hold_velocity, time_range, state_0);

figure()
plot(t, state(:,1), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('2D Trajectory Velocity')

figure()
plot(t, rad2deg(state(:,2)), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Flight Path Angle (deg)')
title('2D Trajectory Flight Path Angle')

figure()
plot(t, state(:,3), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Height (m)')
title('2D Trajectory Height')

figure()
plot(t, state(:,4), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Position (m)')
title('2D Trajectory Position')

figure()
plot(state(:,4), state(:,3), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Position (m)')
ylabel('Height (m)')
title('Trajectory Simulation')