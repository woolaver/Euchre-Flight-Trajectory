%% RK4 Trajectory Driver function, ENSURE CORRECT CSV FILES ARE BEING USED/SAVED BEFORE RUNNING
% Also you have to change the angle of attacks within the function before
% running

clear
clc
close all

height_init = 1;

[~, a, ~, ~] = atmoscoesa(height_init);

Mach_init = 8;
V_init = a*Mach_init;

launch_angle = deg2rad(18.5); %20 for virginia, 17 for Ice Cream Cone, 18.5 for HARV
x_init = 0;

state0 = [V_init; launch_angle; height_init; x_init];

dt = .01; %s

%{
CL_vec_file = "Virginia_blended_CL.csv";
CD_vec_file = "Virginia_blended_CD.csv";
Trajectory_Save_File = "Virginia_Trajectory.csv";

Sref_lift = .439; %m^2
Sref_drag = .035; %m^2
%}

%{
CL_vec_file = "Ice_Cream_Cone_blended_CL.csv";
CD_vec_file = "Ice_Cream_Cone_blended_CD.csv";
Trajectory_Save_File = "Ice_Cream_Cone_Trajectory.csv";

Sref_lift = .253; %m^2
Sref_drag = .0299; %m^2 
%}

%{
CL_vec_file = "HARV_blended_CL.csv";
CD_vec_file = "HARV_blended_CD.csv";
Trajectory_Save_File = "HARV_Trajectory.csv";

Sref_lift = .359; %m^2
Sref_drag = .029; %m^2
%}

CL_vec_file = "Bunny_Bomb_blended_CL.csv";
CD_vec_file = "Bunny_Bomb_blended_CD.csv";
Trajectory_Save_File = "Bunny_Bomb_Trajectory.csv";

Sref_lift = .467; %m^2
Sref_drag = .049; %m^2

[t, state, usefulData] = RK4_3DoF(state0, dt, CL_vec_file, CD_vec_file, Sref_lift, Sref_drag);

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


disp("Final Results: ")
disp("Total Time: " + t(end) + " s")
disp("Total Distance: " + state(4, end) + " m")
disp("Maximum Height: " + max(state(3, :)) + " m")
disp("Landing Mach Number: " + usefulData(5, end))

saveArray = [t' state' usefulData'];

saveTable = array2table(saveArray, 'VariableNames', {'Time', ...
    'Velocity', 'Flight Path Angle', 'Height', ...
    'Distance', 'Angle of Attack', 'CL', 'CD', ...
    'Density', 'Mach Number'});

writetable(saveTable, Trajectory_Save_File);
