%% RK4 Trajectory Driver function, ENSURE CORRECT PROPERTIES FILES ARE BEING USED BEFORE RUNNING

clear
clc
close all

vehicleProperties = blend_sharp_properties();

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

%find indicies where climb, level, glide, and dive end
climbEndIndex = find(state(3, :) >= 15000, 1, 'first');
levelEndIndex = find(state(2, :) <= 0.0015, 1, 'first');
glideEndIndex = find(state(1, levelEndIndex:end) <= vehicleProperties.dive_speed, 1, 'first') + levelEndIndex;

figure()
hold on
plot(t, state(1, :), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), state(1, climbEndIndex), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), state(1, levelEndIndex), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), state(1, glideEndIndex), 'o', 'LineWidth', 3)
legend('Velocity', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity vs. Time')
hold off

figure()
hold on
plot(t, rad2deg(state(2, :)), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), rad2deg(state(2, climbEndIndex)), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), rad2deg(state(2, levelEndIndex)), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), rad2deg(state(2, glideEndIndex)), 'o', 'LineWidth', 3)
legend('Glide Path Angle', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Flight Path Angle (deg)')
title('Flight Path Angle vs. Time')
hold off

figure()
hold on
plot(t, state(3, :), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), state(3, climbEndIndex), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), state(3, levelEndIndex), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), state(3, glideEndIndex), 'o', 'LineWidth', 3)
legend('Height', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Height (m)')
title('Height vs. Time')
hold off

figure()
hold on
plot(t, state(4, :), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), state(4, climbEndIndex), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), state(4, levelEndIndex), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), state(4, glideEndIndex), 'o', 'LineWidth', 3)
legend('Distance', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Distance (m)')
title('Distance vs. Time')
hold off

figure()
hold on
plot(state(4, :), state(3, :), 'LineWidth', 1, 'Color', 'Blue')
plot(state(4, climbEndIndex), state(3, climbEndIndex), 'o', 'LineWidth', 3)
plot(state(4, levelEndIndex), state(3, levelEndIndex), 'o', 'LineWidth', 3)
plot(state(4, glideEndIndex), state(3, glideEndIndex), 'o', 'LineWidth', 3)
legend('Trajectory', 'Climb End', 'Level End', 'Glide End')
xlabel('Distance (m)')
ylabel('Height (m)')
title('Projectile Trajectory')
hold off

figure()
hold on
plot(t, rad2deg(usefulData(1, :)), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), rad2deg(usefulData(1, climbEndIndex)), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), rad2deg(usefulData(1, levelEndIndex)), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), rad2deg(usefulData(1, glideEndIndex)), 'o', 'LineWidth', 3)
legend('Angle of Attack', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Angle of Attack (deg)')
title('Angle of Attack vs. Time')
hold off

figure()
hold on
plot(t, usefulData(2, :), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), usefulData(2, climbEndIndex), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), usefulData(2, levelEndIndex), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), usefulData(2, glideEndIndex), 'o', 'LineWidth', 3)
legend('CL', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('CL')
title('CL vs. Time')
hold off

figure()
hold on
plot(t, usefulData(3, :), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), usefulData(3, climbEndIndex), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), usefulData(3, levelEndIndex), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), usefulData(3, glideEndIndex), 'o', 'LineWidth', 3)
legend('CD', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('CD')
title('CD vs. Time')
hold off

figure()
hold on
plot(t, usefulData(4, :), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), usefulData(4, climbEndIndex), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), usefulData(4, levelEndIndex), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), usefulData(4, glideEndIndex), 'o', 'LineWidth', 3)
legend('Density', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Density (kg/m^3)')
title('Density vs. Time')
hold off

figure()
hold on
plot(t, usefulData(5, :), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), usefulData(5, climbEndIndex), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), usefulData(5, levelEndIndex), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), usefulData(5, glideEndIndex), 'o', 'LineWidth', 3)
legend('Mach Number', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Mach Number')
title('Mach Number vs. Time')
hold off

figure()
hold on
plot(t, usefulData(6, :), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), usefulData(6, climbEndIndex), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), usefulData(6, levelEndIndex), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), usefulData(6, glideEndIndex), 'o', 'LineWidth', 3)
legend('Wall Temperature', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Wall Temperature (K)')
title('Wall Temperature vs. Time')
hold off

figure()
hold on
plot(t, usefulData(7, :), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), usefulData(7, climbEndIndex), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), usefulData(7, levelEndIndex), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), usefulData(7, glideEndIndex), 'o', 'LineWidth', 3)
legend('Internal Temperature', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Internal Temperature (K)')
title('Internal Temperature vs. Time')
hold off

figure()
hold on
plot(t, usefulData(8, :), 'LineWidth', 1, 'Color', 'Blue')
plot(t(climbEndIndex), usefulData(8, climbEndIndex), 'o', 'LineWidth', 3)
plot(t(levelEndIndex), usefulData(8, levelEndIndex), 'o', 'LineWidth', 3)
plot(t(glideEndIndex), usefulData(8, glideEndIndex), 'o', 'LineWidth', 3)
legend('Nose Temperature', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Nose Temperature (K)')
title('Nose Temperature vs. Time')
hold off

disp("Final Results: ")
disp("Total Time: " + t(end) + " s")
disp("Total Distance: " + state(4, end) + " m")
disp("Maximum Height: " + max(state(3, :)) + " m")
disp("Landing Mach Number: " + usefulData(5, end))
disp("Maximum Wall Temperature: " + max(usefulData(6, :)) + " K")
disp("Maximum Internal Temperature: " + max(usefulData(7, :)) + " K")
disp("Maximum Nose Temperature: " + max(usefulData(8, :)) + " K")

Trajectory_Save_File = vehicleProperties.Trajectory_Save_File;

saveArray = [t' state' usefulData'];

saveTable = array2table(saveArray, 'VariableNames', {'Time', ...
    'Velocity', 'FlightPathAngle', 'Height', ...
    'Distance', 'AngleOfAttack', 'CL', 'CD', ...
    'Density', 'Mach Number', 'WallTemperature', 'InternalTemperature', 'NoseTemperature'});

writetable(saveTable, Trajectory_Save_File);
