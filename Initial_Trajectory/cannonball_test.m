%Sweeps initial launch angle using ODE_Model_2D functions

clear
clc
close all

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
    launch_angle = deg2rad(i*4);

    state_0 = [V0; launch_angle; init_height; 0];
    time_range = [0 1000];
    options = odeset('Events', @ground);

    [t, state] = ode45(@ODE_cannonball, time_range, state_0, options);

    % Store results for analysis
    results(i).time = t;
    results(i).V = state(:,1);
    results(i).flight_path = state(:,2);
    results(i).height = state(:,3);
    results(i).position = state(:,4);
end

% Plot the results for each launch angle
figure()
hold on;
for i = 1:10
    plot(results(i).time, results(i).V, 'DisplayName', sprintf('Launch Angle: %d°', i*4));
end
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity vs. Time for Different Launch Angles');
legend show;
hold off

figure()
hold on;
for i = 1:10
    plot(results(i).time, rad2deg(results(i).flight_path), 'DisplayName', sprintf('Launch Angle: %d°', i*4));
end
xlabel('Time (s)');
ylabel('Flight Path Angle (deg)');
title('Flight Path Angle vs. Time for Different Launch Angles');
legend show;
hold off

figure()
hold on;
for i = 1:10
    plot(results(i).time, results(i).height, 'DisplayName', sprintf('Launch Angle: %d°', i*4));
end
xlabel('Time (s)');
ylabel('Height (m)');
title('Height vs. Time for Different Launch Angles');
legend show;
hold off

figure()
hold on;
for i = 1:10
    plot(results(i).time, results(i).position, 'DisplayName', sprintf('Launch Angle: %d°', i*4));
end
xlabel('Time (s)');
ylabel('Position (m)');
title('Position vs. Time for Different Launch Angles');
legend show;
hold off

figure()
hold on;
for i = 1:10
    plot(results(i).position, results(i).height, 'DisplayName', sprintf('Launch Angle: %d°', i*4));
end
xlabel('Position (m)');
ylabel('Height (m)');
title('Trajectory for Different Launch Angles');
legend show;
hold off
