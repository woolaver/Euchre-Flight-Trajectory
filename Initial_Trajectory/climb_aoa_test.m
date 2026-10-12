%loop through multiple angles of attack for climb period to determine
%climb properties
clear
clc
close all

N_alpha = 100;
alpha = deg2rad(linspace(0, 15, N_alpha));
mass = 120; %kg
g = 9.8; %m/s^2

% Initialize arrays to store results
results = zeros(N_alpha, 6);

function [value, isterminal, direction] = climb_end(t, state)
    value = (state(1)*sin(state(2)) <= 0);
    isterminal = 1;
    direction = 0;
end

for i = 1:100

    h0 = 1;
    Mach = 8;
    [~, a, ~, ~] = atmoscoesa(h0);

    V0 = Mach*a;
    launch_angle = deg2rad(25);
    alpha0 = alpha(i);
    state0 = [V0; launch_angle; h0; 0; alpha0];

    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'Events', @climb_end);

    tspan = [0 240];

    [t, state] = ode45(@ODE_lifting_climb, tspan, state0, options);

    energy = .5*state(end, 1)^2*mass + mass*g*state(end, 3);

    results(i, :) = [t(end), state(end, 1), state(end, 2), state(end, 3), state(end, 4), energy];
end

% Analyze and plot the results for each angle of attack
figure()
plot(rad2deg(alpha), results(:, 1), 'LineWidth', 1, 'Color', 'Blue');
xlabel('Angle of Attack (deg)');
ylabel('Climb Time (s)');
title('Climb Time vs Angle of Attack');

figure()
plot(rad2deg(alpha), results(:, 2), 'LineWidth', 1, 'Color', 'Blue');
xlabel('Angle of Attack (deg)');
ylabel('Final Velocity (m/s)');
title('Final Velocity vs Angle of Attack');

figure()
plot(rad2deg(alpha), rad2deg(results(:, 3)), 'LineWidth', 1, 'Color', 'Blue');
xlabel('Angle of Attack (deg)');
ylabel('Final Flight Path Angle (deg)');
title('Final Flight Path Angle vs Angle of Attack');

figure()
plot(rad2deg(alpha), results(:, 4), 'LineWidth', 1, 'Color', 'Blue');
xlabel('Angle of Attack (deg)');
ylabel('Final Height (m)');
title('Final Height vs Angle of Attack');

figure()
plot(rad2deg(alpha), results(:, 5), 'LineWidth', 1, 'Color', 'Blue');
xlabel('Angle of Attack (deg)');
ylabel('Final Position (m)');
title('Final Position vs Angle of Attack');

figure()
plot(rad2deg(alpha), results(:, 6), 'LineWidth', 1, 'Color', 'Blue');
xlabel('Angle of Attack (deg)');
ylabel('Final Energy (J)');
title('Final Energy vs Angle of Attack');

disp('Maximum final velocity and alpha: ')
[max_V, alpha_V] = max(results(:,2));
disp(max_V)
disp(rad2deg(alpha(alpha_V)))

disp('Maximum final position and alpha: ')
[max_x, alpha_x] = max(results(:,5));
disp(max_x)
disp(rad2deg(alpha(alpha_x)))

disp('Maximum final height and alpha: ')
[max_h, alpha_h] = max(results(:,4));
disp(max_h)
disp(rad2deg(alpha(alpha_h)))

disp('Maximum final energy and alpha: ')
[max_E, alpha_E] = max(results(:,6));
disp(max_E)
disp(rad2deg(alpha(alpha_E)))