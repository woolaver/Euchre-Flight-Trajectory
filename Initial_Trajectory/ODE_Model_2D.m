%2D ODE Model for simple trajectory analysis, model is a 3DOF point mass glider
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
launch_angle = deg2rad(28);
init_height = 1;
[~, a, ~, ~] = atmoscoesa(init_height); %m/s
V0 = Mach_init*a;

state_0_climb= [V0; launch_angle; init_height; 0];
time_range_climb = [0 200];
solver_options = odeset('RelTol', 1e-7, 'AbsTol', [1e-6 1e-9 1e-5 1e-5], ...
                        'MaxStep', 1);
options_climb = odeset(solver_options, 'Events', @cruise_start);

[t_climb, state_climb, te_climb] = ode45(@ODE_climb, time_range_climb, state_0_climb, options_climb);

if isempty(te_climb)
    error('Climb never reached the apex before time limit.');
end

state_0_cruise = state_climb(end, :);
time_range_cruise = [t_climb(end), t_climb(end) + 240];
options_cruise = odeset(solver_options, 'Events', @hold_speed);

[t_cruise, state_cruise, te_cruise] = ode45(@ODE_cruise, time_range_cruise, state_0_cruise, options_cruise);

if isempty(te_cruise)
    error('Cruise never reached the 7000 m hold speed altitude.');
end

state_0_hold = state_cruise(end, :);
time_range_hold = [t_cruise(end), t_cruise(end) + 240];
options_hold = odeset(solver_options, 'Events', @dive);
[t_hold, state_hold, te_hold] = ode45(@ODE_hold_velocity, time_range_hold, state_0_hold, options_hold);

state_0_dive = state_hold(end, :);
state_0_dive(2) = deg2rad(-90);
time_range_dive = [t_hold(end), t_hold(end) + 240];
options_dive = odeset(solver_options, 'Events', @ground);
[t_dive, state_dive, te_dive] = ode45(@ODE_dive, time_range_dive, state_0_dive, options_dive);

t_total = [t_climb; t_cruise(2:end); t_hold(2:end); t_dive(2:end)];
state_total = [state_climb; state_cruise(2:end,:); ...
               state_hold(2:end,:); state_dive(2:end,:)];

phase_names = {'Climb', 'Cruise', 'Hold', 'Dive'};
phase_times = {t_climb, t_cruise, t_hold, t_dive};
phase_states = {state_climb, state_cruise, state_hold, state_dive};
if ~isempty(t_dive)
    phase_names{end+1} = 'Dive';
    phase_times{end+1} = t_dive;
    phase_states{end+1} = state_dive;
end
fprintf('Phase end:        t (s)      V (m/s)   gamma (deg)      h (m)        x (m)\n');
for k = 1:numel(phase_names)
    last = phase_states{k}(end,:);
    fprintf('%-10s %11.3f %12.3f %12.3f %12.3f %12.3f\n', ...
        phase_names{k}, phase_times{k}(end), last(1), rad2deg(last(2)), last(3), last(4));
end
%fprintf('%s\n', termination_reason);

% Use the same state column for each curve and its phase-end markers.
y_labels = {'Speed (m/s)', 'Flight Path Angle (deg)', ...
            'Height (m)', 'Horizontal Position (m)'};
for column = 1:4
    scale = 1;
    if column == 2
        scale = 180/pi;
    end
    figure()
    hold on
    plot(t_total, scale*state_total(:,column), 'b', ...
         'DisplayName', y_labels{column});
    for k = 1:numel(phase_names)
        plot(phase_times{k}(end), scale*phase_states{k}(end,column), ...
             'o', 'LineWidth', 2, 'DisplayName', [phase_names{k} ' End']);
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

%write final state vector to a csv file

csv_file = [t_total'; state_total'];

writematrix(csv_file, "justin_sucks_cock.csv")

%{
% Events detect zero crossings, so check an already-low entry speed before
% integrating. Keep one boundary sample for the existing phase bookkeeping.
[entry_values, ~, ~] = dive(time_range_hold(1), state_0_hold);
if entry_values(2) <= 0 || entry_values(1) <= 0
    t_hold = time_range_hold(1);
    state_hold = state_0_hold;
    te_hold = t_hold;
    ye_hold = state_hold;
    if entry_values(2) <= 0
        ie_hold = 2; % Ground contact takes priority over starting a dive.
    else
        ie_hold = 1; % At or below target speed: start the dive immediately.
    end
else
    [t_hold, state_hold, te_hold, ye_hold, ie_hold] = ...
        ode45(@ODE_hold_velocity, ...
        time_range_hold, ...
        state_0_hold, ...
        options_hold);
end

if isempty(te_hold)
    error('Trajectory:HoldTimeout', ...
          'Hold reached its time limit before the target speed or ground.');
end

t_dive = zeros(0,1);
state_dive = zeros(0,4);
te_dive = [];
ye_dive = zeros(0,4);
ie_dive = [];
if ie_hold(end) == 2 || state_hold(end,3) <= 0
    termination_reason = 'Ground reached during hold; dive skipped';
else
    state_0_dive = state_hold(end, :);
    state_0_dive(2) = deg2rad(-40);
    time_range_dive = [t_hold(end), t_hold(end) + 240];
    options_dive = odeset(solver_options, 'Events', @ground);
    [t_dive, state_dive, te_dive, ye_dive, ie_dive] = ...
        ode45(@ODE_dive, time_range_dive, state_0_dive, options_dive);
    if isempty(te_dive)
        error('Trajectory:DiveTimeout', ...
              'Dive reached its time limit before ground contact.');
    end
    termination_reason = 'Ground reached during dive';
end


t_total = [t_climb; t_cruise(2:end); t_hold(2:end); t_dive(2:end)];
state_total = [state_climb; state_cruise(2:end,:); ...
               state_hold(2:end,:); state_dive(2:end,:)];

phase_names = {'Climb', 'Cruise', 'Hold'};
phase_times = {t_climb, t_cruise, t_hold};
phase_states = {state_climb, state_cruise, state_hold};
if ~isempty(t_dive)
    phase_names{end+1} = 'Dive';
    phase_times{end+1} = t_dive;
    phase_states{end+1} = state_dive;
end
fprintf('Phase end:        t (s)      V (m/s)   gamma (deg)      h (m)        x (m)\n');
for k = 1:numel(phase_names)
    last = phase_states{k}(end,:);
    fprintf('%-10s %11.3f %12.3f %12.3f %12.3f %12.3f\n', ...
        phase_names{k}, phase_times{k}(end), last(1), rad2deg(last(2)), last(3), last(4));
end
fprintf('%s\n', termination_reason);

% Use the same state column for each curve and its phase-end markers.
y_labels = {'Speed (m/s)', 'Flight Path Angle (deg)', ...
            'Height (m)', 'Horizontal Position (m)'};
for column = 1:4
    scale = 1;
    if column == 2
        scale = 180/pi;
    end
    figure()
    hold on
    plot(t_total, scale*state_total(:,column), 'b', ...
         'DisplayName', y_labels{column});
    for k = 1:numel(phase_names)
        plot(phase_times{k}(end), scale*phase_states{k}(end,column), ...
             'o', 'LineWidth', 2, 'DisplayName', [phase_names{k} ' End']);
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
%}