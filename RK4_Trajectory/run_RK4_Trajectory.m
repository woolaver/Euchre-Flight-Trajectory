%% RK4 Trajectory Driver function, ENSURE CORRECT CSV FILES ARE BEING USED/SAVED BEFORE RUNNING
% Also you have to change the angle of attacks within the function before
% running

clear
clc
close all

vehicleProperties = blend_round_properties();

height_init = 1;

[~, a, ~, ~] = atmoscoesa(height_init);

Mach_init = 8;
V_init = a*Mach_init;

launch_angle = vehicleProperties.launch_angle; 
x_init = 0;

state0 = [V_init; launch_angle; height_init; x_init];

[t, state, usefulData] = RK4_3DoF(state0, vehicleProperties);

%sometimes final row will have NaN's due to certain values going
%negative or something, if this happens just remove the last row
if(isnan(state(4, end)))
    t = t(:, 1:end-1);
    state = state(:, 1:end-1);
    usefulData = usefulData(:, 1:end-1);
end

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

figure()
plot(t, usefulData(6, :), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Wall Temperature (K)')
title('Wall Temperature vs. Time')

figure()
plot(t, usefulData(7, :), 'LineWidth', 1, 'Color', 'Blue')
xlabel('Time (s)')
ylabel('Internal Temperature (K)')
title('Internal Temperature vs. Time')

disp("Final Results: ")
disp("Total Time: " + t(end) + " s")
disp("Total Distance: " + state(4, end) + " m")
disp("Maximum Height: " + max(state(3, :)) + " m")
disp("Landing Mach Number: " + usefulData(5, end))
disp("Maximum Wall Temperature: " + max(usefulData(6, :)) + " K")
disp("Maximum Internal Temperature: " + max(usefulData(7, :)) + " K")

saveArray = [t' state' usefulData'];

saveTable = array2table(saveArray, 'VariableNames', {'Time', ...
    'Velocity', 'FlightPathAngle', 'Height', ...
    'Distance', 'AngleOfAttack', 'CL', 'CD', ...
    'Density', 'Mach Number', 'WallTemperature', 'InternalTemperature'});

writetable(saveTable, Trajectory_Save_File);
