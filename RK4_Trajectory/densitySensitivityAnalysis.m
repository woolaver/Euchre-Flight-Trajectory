%Complete a density sensitivity analysis on the trajectory to determine how
%much changes in density can effect the trajectory of the vehicle

clc
clear
close all

vehicleProperties = half_cone_properties();

height_init = 1;

[~, a, ~, ~] = atmoscoesa(height_init);

Mach_init = 8;
V_init = a*Mach_init;

launch_angle = vehicleProperties.launch_angle; 
x_init = 0;

state0 = [V_init; launch_angle; height_init; x_init];

%nominal run
tempOffset = 0;
[t_nom, state_nom, usefulData_nom] = RK4_3DoF(state0, vehicleProperties, tempOffset);

%low density
tempOffset = 10; %K
[t_low, state_low, usefulData_low] = RK4_3DoF(state0, vehicleProperties, tempOffset);

%high density
tempOffset = -10; %K
[t_high, state_high, usefulData_high] = RK4_3DoF(state0, vehicleProperties, tempOffset);

figure()
hold on
plot(t_nom, state_nom(1, :), 'LineWidth', 1)
plot(t_low, state_low(1, :), 'LineWidth', 1)
plot(t_high, state_high(1, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity vs. Time')
hold off

figure()
hold on
plot(t_nom, rad2deg(state_nom(2, :)), 'LineWidth', 1)
plot(t_low, rad2deg(state_low(2, :)), 'LineWidth', 1)
plot(t_high, rad2deg(state_high(2, :)), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Flight Path Angle (deg)')
title('Flight Path Angle vs. Time')
hold off

figure()
hold on
plot(t_nom, state_nom(3, :), 'LineWidth', 1)
plot(t_low, state_low(3, :), 'LineWidth', 1)
plot(t_high, state_high(3, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Height (m)')
title('Height vs. Time')
hold off

figure()
hold on
plot(t_nom, state_nom(4, :), 'LineWidth', 1)
plot(t_low, state_low(4, :), 'LineWidth', 1)
plot(t_high, state_high(4, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Distance (m)')
title('Distance vs. Time')
hold off

figure()
hold on
plot(state_nom(4, :), state_nom(3, :), 'LineWidth', 1)
plot(state_low(4, :), state_low(3, :), 'LineWidth', 1)
plot(state_high(4, :), state_high(3, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Distance (m)')
ylabel('Height (m)')
title('Projectile Trajectory')
hold off

figure()
hold on
plot(t_nom, rad2deg(usefulData_nom(1, :)), 'LineWidth', 1)
plot(t_low, rad2deg(usefulData_low(1, :)), 'LineWidth', 1)
plot(t_high, rad2deg(usefulData_high(1, :)), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Body Angle of Attack (deg)')
title('Body Angle of Attack vs. Time')
hold off

figure()
hold on
plot(t_nom, rad2deg(usefulData_nom(2, :)), 'LineWidth', 1)
plot(t_low, rad2deg(usefulData_low(2, :)), 'LineWidth', 1)
plot(t_high, rad2deg(usefulData_high(2, :)), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Fin Angle of Attack (deg)')
title('Fin Angle of Attack vs. Time')
hold off

figure()
hold on
plot(t_nom, usefulData_nom(3, :), 'LineWidth', 1)
plot(t_low, usefulData_low(3, :), 'LineWidth', 1)
plot(t_high, usefulData_high(3, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('CL')
title('CL vs. Time')
hold off

figure()
hold on
plot(t_nom, usefulData_nom(4, :), 'LineWidth', 1)
plot(t_low, usefulData_low(4, :), 'LineWidth', 1)
plot(t_high, usefulData_high(4, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('CD')
title('CD vs. Time')
hold off

figure()
hold on
plot(t_nom, usefulData_nom(5, :), 'LineWidth', 1)
plot(t_low, usefulData_low(5, :), 'LineWidth', 1)
plot(t_high, usefulData_high(5, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Density (kg/m^3)')
title('Density vs. Time')
hold off

figure()
hold on
plot(t_nom, usefulData_nom(6, :), 'LineWidth', 1)
plot(t_low, usefulData_low(6, :), 'LineWidth', 1)
plot(t_high, usefulData_high(6, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Mach Number')
title('Mach Number vs. Time')
hold off

figure()
hold on
plot(t_nom, usefulData_nom(7, :), 'LineWidth', 1)
plot(t_low, usefulData_low(7, :), 'LineWidth', 1)
plot(t_high, usefulData_high(7, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Wall Temperature (K)')
title('Wall Temperature vs. Time')
hold off

figure()
hold on
plot(t_nom, usefulData_nom(8, :), 'LineWidth', 1)
plot(t_low, usefulData_low(8, :), 'LineWidth', 1)
plot(t_high, usefulData_high(8, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Internal Temperature (K)')
title('Internal Temperature vs. Time')
hold off

figure()
hold on
plot(t_nom, usefulData_nom(9, :), 'LineWidth', 1)
plot(t_low, usefulData_low(9, :), 'LineWidth', 1)
plot(t_high, usefulData_high(9, :), 'LineWidth', 1)
legend('Nominal Density', 'Low Density', 'High Density')
xlabel('Time (s)')
ylabel('Nose Temperature (K)')
title('Nose Temperature vs. Time')
hold off

disp("Final Results: ")
disp("Final Results: ")
fprintf('%-30s %16s %16s %16s\n', '', 'Nominal Density', 'Low Density', 'High Density');
fprintf('%-30s %16.2f %16.2f %16.2f\n', 'Total Time (s)', t_nom(end), t_low(end), t_high(end));
fprintf('%-30s %16.2f %16.2f %16.2f\n', 'Total Distance (m)', state_nom(4, end), state_low(4, end), state_high(4, end));
fprintf('%-30s %16.2f %16.2f %16.2f\n', 'Maximum Height (m)', max(state_nom(3, :)), max(state_low(3, :)), max(state_high(3, :)));
fprintf('%-30s %16.2f %16.2f %16.2f\n', 'Landing Mach Number', usefulData_nom(6, end), usefulData_low(6, end), usefulData_high(6, end));
fprintf('%-30s %16.2f %16.2f %16.2f\n', 'Maximum Wall Temperature (K)', max(usefulData_nom(7, :)), max(usefulData_low(7, :)), max(usefulData_high(7, :)));
fprintf('%-30s %16.2f %16.2f %16.2f\n', 'Maximum Internal Temp (K)', max(usefulData_nom(8, :)), max(usefulData_low(8, :)), max(usefulData_high(8, :)));
fprintf('%-30s %16.2f %16.2f %16.2f\n', 'Maximum Nose Temperature (K)', max(usefulData_nom(9, :)), max(usefulData_low(9, :)), max(usefulData_high(9, :)));