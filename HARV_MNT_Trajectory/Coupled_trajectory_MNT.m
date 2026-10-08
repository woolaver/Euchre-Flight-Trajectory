%Coupled trajectory model with MNT for more accurate trajectory model

clear
clc
close all

function [value, isterminal, direction] = cruise_start(t, state)
    value = state(1)*sin(state(2));
    isterminal = 1;
    direction = -1;
end

function [value, isterminal, direction] = ground(t, state)
    value = state(3);
    isterminal = 1;
    direction = -1;
end

function [value, isterminal, direction] = dive(t, state)
    value = state(3) - 10000;
    isterminal = 1;
    direction = -1;
end

function [value, isterminal, direction] = hold_speed(t, state)

    [~, a, ~, ~] = atmoscoesa(0);
    impact_speed = 3*a;

    value = [ ...
        state(1) - impact_speed;  % Event 1: Mach 3
        state(3)                  % Event 2: ground
    ];

    isterminal = [1; 1];

    direction = [ ...
         0;   % reaching the threshold from either direction
        -1    % descending through h = 0
    ];

end

Mach_init = 8;
launch_angle = deg2rad(15);
init_height = 1;
[~, a, ~, ~] = atmoscoesa(init_height); %m/s
V0 = Mach_init*a;

state_0_climb= [V0; launch_angle; init_height; 0];
time_range_climb = [0 500];
solver_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 1);
options_climb = odeset(solver_options, 'Events', @cruise_start);

disp("before climb")

[t_climb, state_climb, te_climb] = ode45(@ODE_climb_mnt, time_range_climb, state_0_climb, options_climb);

if isempty(te_climb)
    error('Climb never reached the apex before time limit.');
end

disp("climb end")

state_0_cruise = state_climb(end, :);
time_range_cruise = [t_climb(end), t_climb(end) + 500];
options_cruise = odeset(solver_options, 'Events', @hold_speed);

[t_cruise, state_cruise, te_cruise] = ode45(@ODE_cruise_mnt, time_range_cruise, state_0_cruise, options_cruise);

if isempty(te_cruise)
    error('Cruise never reached the 7000 m hold speed altitude.');
end

disp("cruise end")

state_0_hold = state_cruise(end, :);
time_range_hold = [t_cruise(end), t_cruise(end) + 500];
options_hold = odeset(solver_options, 'Events', @ground);
[t_hold, state_hold, te_hold] = ode45( ...
    @ODE_hold_velocity_mnt, time_range_hold, state_0_hold, options_hold);

disp("hold velocity end")


%{
state_0_dive = state_hold(end, :);
state_0_dive(2) = deg2rad(-90);
time_range_dive = [t_hold(end), t_hold(end) + 500];
options_dive = odeset(solver_options, 'Events', @ground);
[t_dive, state_dive, te_dive] = ode45(@ODE_dive_mnt, time_range_dive, state_0_dive, options_dive);

disp("dive end")
%}

t_total = [t_climb; t_cruise(2:end); t_hold(2:end, :)];
state_total = [state_climb; state_cruise(2:end,:); state_hold(2:end, :)];

phase_names = {'Climb', 'Cruise', 'Hold'};
phase_times = {t_climb, t_cruise, t_hold};
phase_states = {state_climb, state_cruise, state_hold};

fprintf('Phase end:        t (s)      V (m/s)   gamma (deg)      h (m)        x (m)\n');
for k = 1:numel(phase_names)
    last = phase_states{k}(end,:);
    fprintf('%-10s %11.3f %12.3f %12.3f %12.3f %12.3f\n', phase_names{k}, phase_times{k}(end), last(1), rad2deg(last(2)), last(3), last(4));
end
%fprintf('%s\n', termination_reason);

% Use the same state column for each curve and its phase-end markers.
y_labels = {'Speed (m/s)', 'Flight Path Angle (deg)', 'Height (m)', 'Horizontal Position (m)'};
for column = 1:4
    scale = 1;
    if column == 2
        scale = 180/pi;
    end
    figure()
    hold on
    plot(t_total, scale*state_total(:,column), 'b','DisplayName', y_labels{column});
    for k = 1:numel(phase_names)
        plot(phase_times{k}(end), scale*phase_states{k}(end,column), 'o', 'LineWidth', 2, 'DisplayName', [phase_names{k} ' End']);
    end
    xlabel('Time (s)');
    ylabel(y_labels{column});
    title(['2D Trajectory: ' y_labels{column}]);
    legend('show');
    hold off
end

figure()
hold on
plot(state_total(:,4), state_total(:,3), 'b', 'DisplayName', 'Trajectory');
for k = 1:numel(phase_names)
    plot(phase_states{k}(end,4), phase_states{k}(end,3), ...
         'o', 'LineWidth', 2, 'DisplayName', [phase_names{k} ' End']);
end
xlabel('Horizontal Position (m)');
ylabel('Height (m)');
title('Trajectory Simulation');
legend('show', Location='best');
hold off

disp("Max Height: " + max(state_total(:, 3)) + " m")
disp("Final Range: " + state_total(end, 4) + " m")
disp("Total Time: " + t_total(end) + " s")
disp("Mach Number at Impact: " + state_total(end, 1)/a)

csv_file = [t_total'; state_total'];

writematrix(csv_file, "HARV_MNT_Trajectory.csv")

[alpha, Mach, CL, CD] = getTrajectoryData("HARV_MNT_Trajectory.csv", t_climb(end), t_cruise(end));

index_climb_end = length(t_climb);
index_cruise_end = length(t_cruise) + index_climb_end - 1;
index_hold_end = length(t_hold) + index_cruise_end - 1;

figure()
hold on
plot(t_total, rad2deg(alpha), 'LineWidth', 1, 'Color', 'Blue')
plot(t_climb(end), rad2deg(alpha(index_climb_end)), 'o', 'LineWidth', 2)
plot(t_cruise(end), rad2deg(alpha(index_cruise_end)), 'o', 'LineWidth', 2)
plot(t_hold(end), rad2deg(alpha(index_hold_end)), 'o', 'LineWidth', 2)
xlabel('Time (s)')
ylabel("\alpha (deg)")
legend('Angle of Attack', 'Climb End', 'Cruise End', 'Hold End')
title('Angle of Attack vs. Time')
hold off

figure()
hold on
plot(t_total, Mach, 'LineWidth', 1, 'Color', 'Blue')
plot(t_climb(end), Mach(index_climb_end), 'o', 'LineWidth', 2)
plot(t_cruise(end), Mach(index_cruise_end), 'o', 'LineWidth', 2)
plot(t_hold(end), Mach(index_hold_end), 'o', 'LineWidth', 2)
xlabel('Time (s)')
ylabel("Mach Number")
legend('Mach Number', 'Climb End', 'Cruise End', 'Hold End')
title('Mach Number vs. Time')
hold off

figure()
hold on
plot(t_total, CL, 'LineWidth', 1, 'Color', 'Blue')
plot(t_climb(end), CL(index_climb_end), 'o', 'LineWidth', 2)
plot(t_cruise(end), CL(index_cruise_end), 'o', 'LineWidth', 2)
plot(t_hold(end), CL(index_hold_end), 'o', 'LineWidth', 2)
xlabel('Time (s)')
ylabel("CL")
legend('CL', 'Climb End', 'Cruise End', 'Hold End')
title('CL vs. Time')
hold off

figure()
hold on
plot(t_total, CD, 'LineWidth', 1, 'Color', 'Blue')
plot(t_climb(end), CD(index_climb_end), 'o', 'LineWidth', 2)
plot(t_cruise(end), CD(index_cruise_end), 'o', 'LineWidth', 2)
plot(t_hold(end), CD(index_hold_end), 'o', 'LineWidth', 2)
xlabel('Time (s)')
ylabel("CD")
legend('CD', 'Climb End', 'Cruise End', 'Hold End')
title('CD vs. Time')
hold off