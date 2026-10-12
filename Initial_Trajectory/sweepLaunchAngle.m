%Sweeps initial launch angle using ODE_Model_2D functions

clear
clc
close all

%function to stop integration when h_dot is negative for the first time
%(switch from climb to cruise state)
function [value, isterminal, direction] = cruise_start(t, state)
    value = (state(1)*sin(state(2)) <= 0);
    isterminal = 1;
    direction = 0;
end

%event to stop integration when the ground is hit
function [value, isterminal, direction] = ground(t, state)
    value = (state(3) <= 0);
    isterminal = 1;
    direction = 0;
end

Mach_init = 8;
init_height = 1;
[~, a, ~, ~] = atmoscoesa(init_height); %m/s
V0 = Mach_init*a;

results.time = zeros(10);
results.V = zeros(10);
results.flight_path = zeros(10);
results.height = zeros(10);
results.position = zeros(10);

for i = 1:10
    launch_angle = deg2rad(i*3);

    state_0_climb= [V0; launch_angle; init_height; 0];
    time_range_climb = [0 100];
    options_climb = odeset('Events', @cruise_start);

    [t_climb, state_climb] = ode45(@ODE_climb, time_range_climb, state_0_climb, options_climb);

    state_0_cruise = state_climb(end, :);
    time_range_cruise = [t_climb(end), 240];
    options_cruise = odeset('Events', @ground);

    [t_cruise, state_cruise] = ode45(@ODE_cruise, time_range_cruise, state_0_cruise, options_cruise);

    t_total = [t_climb; t_cruise(2:end)];        
    state_total = [state_climb; state_cruise(2:end,:)];

    % Store results for analysis
    results(i).time = t_total;
    results(i).V = state_total(:,1);
    results(i).flight_path = state_total(:,2);
    results(i).height = state_total(:,3);
    results(i).position = state_total(:,4);
end

% Plot the results for each launch angle
figure()
hold on;
for i = 1:10
    plot(results(i).time, results(i).V, 'DisplayName', sprintf('Launch Angle: %d°', i*3));
end
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity vs. Time for Different Launch Angles');
legend show;
hold off

figure()
hold on;
for i = 1:10
    plot(results(i).time, rad2deg(results(i).flight_path), 'DisplayName', sprintf('Launch Angle: %d°', i*3));
end
xlabel('Time (s)');
ylabel('Flight Path Angle (deg)');
title('Flight Path Angle vs. Time for Different Launch Angles');
legend show;
hold off

figure()
hold on;
for i = 1:10
    plot(results(i).time, results(i).height, 'DisplayName', sprintf('Launch Angle: %d°', i*3));
end
xlabel('Time (s)');
ylabel('Height (m)');
title('Height vs. Time for Different Launch Angles');
legend show;
hold off

figure()
hold on;
for i = 1:10
    plot(results(i).time, results(i).position, 'DisplayName', sprintf('Launch Angle: %d°', i*3));
end
xlabel('Time (s)');
ylabel('Position (m)');
title('Position vs. Time for Different Launch Angles');
legend show;
hold off

figure()
hold on;
for i = 1:10
    plot(results(i).position, results(i).height, 'DisplayName', sprintf('Launch Angle: %d°', i*3));
end
xlabel('Position (m)');
ylabel('Height (m)');
title('Trajectory for Different Launch Angles');
legend show;
hold off
