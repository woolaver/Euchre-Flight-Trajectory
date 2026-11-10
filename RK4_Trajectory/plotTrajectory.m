%Read in saved trajectory file and create plots

clear
clc
close all

vehicleProperties = half_cone_properties();

trajectoryData = readtable('csv_files\Half_Cone_Trajectory.csv');

time = trajectoryData.Time;
velocity = trajectoryData.Velocity;
flightPathAngle = rad2deg(trajectoryData.FlightPathAngle);
height = trajectoryData.Height;
distance = trajectoryData.Distance;
angleOfAttack = rad2deg(trajectoryData.AngleOfAttack);
CL = trajectoryData.CL;
CD = trajectoryData.CD;
density = trajectoryData.Density;
MachNumber = trajectoryData.MachNumber;
wallTemperature = trajectoryData.WallTemperature;
internalTemperature = trajectoryData.InternalTemperature;
noseTemperature = trajectoryData.NoseTemperature;

%find indicies where climb, level, glide, and dive end
climbEndIndex = find(height >= 15000, 1, 'first');
levelEndIndex = find(flightPathAngle <= rad2deg(0.0015), 1, 'first');
glideEndIndex = find(velocity(levelEndIndex:end) <= vehicleProperties.dive_speed, 1, 'first') + levelEndIndex;

figure()
hold on
plot(time, velocity, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), velocity(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), velocity(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), velocity(glideEndIndex), 'o', 'LineWidth', 3)
legend('Velocity', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity vs. Time')
hold off

figure()
hold on
plot(time, flightPathAngle, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), flightPathAngle(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), flightPathAngle(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), flightPathAngle(glideEndIndex), 'o', 'LineWidth', 3)
legend('Glide Path Angle', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Flight Path Angle (deg)')
title('Flight Path Angle vs. Time')
hold off

figure()
hold on
plot(time, height, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), height(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), height(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), height(glideEndIndex), 'o', 'LineWidth', 3)
legend('Height', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Height (m)')
title('Height vs. Time')
hold off

figure()
hold on
plot(time, distance, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), distance(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), distance(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), distance(glideEndIndex), 'o', 'LineWidth', 3)
legend('Distance', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Distance (m)')
title('Distance vs. Time')
hold off

figure()
hold on
plot(distance, height, 'LineWidth', 1, 'Color', 'Blue')
plot(distance(climbEndIndex), height(climbEndIndex), 'o', 'LineWidth', 3)
plot(distance(levelEndIndex), height(levelEndIndex), 'o', 'LineWidth', 3)
plot(distance(glideEndIndex), height(glideEndIndex), 'o', 'LineWidth', 3)
legend('Trajectory', 'Climb End', 'Level End', 'Glide End')
xlabel('Distance (m)')
ylabel('Height (m)')
title('Projectile Trajectory')
hold off

figure()
hold on
plot(time, angleOfAttack, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), angleOfAttack(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), angleOfAttack(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), angleOfAttack(glideEndIndex), 'o', 'LineWidth', 3)
legend('Angle of Attack', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Angle of Attack (deg)')
title('Angle of Attack vs. Time')
hold off

figure()
hold on
plot(time, CL, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), CL(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), CL(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), CL(glideEndIndex), 'o', 'LineWidth', 3)
legend('CL', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('CL')
title('CL vs. Time')
hold off

figure()
hold on
plot(time, CD, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), CD(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), CD(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), CD(glideEndIndex), 'o', 'LineWidth', 3)
legend('CD', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('CD')
title('CD vs. Time')
hold off

figure()
hold on
plot(time, density, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), density(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), density(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), density(glideEndIndex), 'o', 'LineWidth', 3)
legend('Density', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Density (kg/m^3)')
title('Density vs. Time')
hold off

figure()
hold on
plot(time, MachNumber, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), MachNumber(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), MachNumber(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), MachNumber(glideEndIndex), 'o', 'LineWidth', 3)
legend('Mach Number', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Mach Number')
title('Mach Number vs. Time')
hold off

figure()
hold on
plot(time, wallTemperature, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), wallTemperature(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), wallTemperature(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), wallTemperature(glideEndIndex), 'o', 'LineWidth', 3)
legend('Wall Temperature', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Wall Temperature (K)')
title('Wall Temperature vs. Time')
hold off

figure()
hold on
plot(time, internalTemperature, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), internalTemperature(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), internalTemperature(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), internalTemperature(glideEndIndex), 'o', 'LineWidth', 3)
legend('Internal Temperature', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Internal Temperature (K)')
title('Internal Temperature vs. Time')
hold off

figure()
hold on
plot(time, noseTemperature, 'LineWidth', 1, 'Color', 'Blue')
plot(time(climbEndIndex), noseTemperature(climbEndIndex), 'o', 'LineWidth', 3)
plot(time(levelEndIndex), noseTemperature(levelEndIndex), 'o', 'LineWidth', 3)
plot(time(glideEndIndex), noseTemperature(glideEndIndex), 'o', 'LineWidth', 3)
legend('Nose Temperature', 'Climb End', 'Level End', 'Glide End')
xlabel('Time (s)')
ylabel('Nose Temperature (K)')
title('Nose Temperature vs. Time')
hold off

disp("Final Results: ")
disp("Total Time: " + time(end) + " s")
disp("Total Distance: " + distance(end) + " m")
disp("Maximum Height: " + max(height) + " m")
disp("Landing Mach Number: " + MachNumber(end))
disp("Maximum Wall Temperature: " + max(wallTemperature(:)) + " K")
disp("Maximum Internal Temperature: " + max(internalTemperature(:)) + " K")
disp("Maximum Nose Temperature: " + max(noseTemperature(:)) + " K")
