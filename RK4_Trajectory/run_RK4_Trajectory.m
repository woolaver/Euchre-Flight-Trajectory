%RK4 Trajectory Driver function

clear
clc
close all

height_init = 1;

[~, a, ~, ~] = atmoscoesa(height_init);

Mach_init = 8;
V_init = a*Mach_init;

launch_angle = deg2rad(20);
x_init = 0;

state0 = [V_init; launch_angle; height_init; x_init];

dt = .1; %s

CL_vec_file = "Virginia_blended_CL.csv";
CD_vec_file = "Virginia_blended_CD.csv";

Sref_lift = .439; %m^2
Sref_drag = .035; %m^2

[t, state, usefulData] = RK4_3DoF(state0, dt, CL_vec_file, CD_vec_file, Sref_lift, Sref_drag);

% Plot results

figure()
plot(t, state(1, :), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity vs. Time')

figure()
plot(t, rad2deg(state(2, :)), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Flight Path Angle (deg)')
title('Flight Path Angle vs. Time')

figure()
plot(t, state(3, :), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Height (m)')
title('Height vs. Time')

figure()
plot(t, state(4, :), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Distance (m)')
title('Distance vs. Time')

figure()
plot(state(4, :), state(3, :), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Distance (m)')
ylabel('Height (m)')
title('Projectile Trajectory')

figure()
plot(t, rad2deg(usefulData(1, :)), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Angle of Attack (deg)')
title('Angle of Attack vs. Time')

figure()
plot(t, usefulData(2, :), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('CL')
title('CL vs. Time')

figure()
plot(t, usefulData(3, :), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('CD')
title('CD vs. Time')

figure()
plot(t, usefulData(4, :), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Density (kg/m^3)')
title('Density vs. Time')

figure()
plot(t, usefulData(5, :), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Mach Number')
title('Mach Number vs. Time')


disp("Final Results: ")
disp("Total Time: " + t(end) + " s")
disp("Total Distance: " + state(4, end) + " m")
disp("Maximum Height: " + max(state(3, :)) + " m")
disp("Landing Mach Number: " + usefulData(5, end))

saveArray = [t state usefulData];

saveTable = array2table(saveArray, 'VariableNames', {'Time (s)', ...
    'Velocity (m/s)', 'Flight Path Angle (rad)', 'Height (m)', ...
    'Distance (m)', 'Angle of Attack (rad)', 'CL', 'CD', ...
    'Density (kg/m^3)', 'Mach Number'});

writeTable(saveTable, 'Virginia_Trajectory.csv');