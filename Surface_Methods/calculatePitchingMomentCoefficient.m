%Use Cp distribution from Newtownian Theory to calculate pitching moment
%coefficient from Mach 2 to Mach 8

clc
close all
clear

Mach_vec = linspace(2, 8, 7);
alpha_vec = deg2rad(linspace(-30, 30, 200));
Sref_lift = .439; %m^2
Sref_drag = .035; %m^2
c = 2; %m, reference length

x_cg = .6*c; %assumption for now, 60% of body

gamma = 1.4;

stlFile = "csv_files\Virginia.stl";

Cm_cg = zeros([length(Mach_vec), length(alpha_vec)]);

for i = 1:length(Mach_vec)
    for j = 1:length(alpha_vec)
        [CL, ~ , ~,x_cp] = blended_driver(stlFile, alpha_vec(j), Mach_vec(i), Sref_lift, Sref_drag, gamma);
        Cm_cg(i, j) = -1*((x_cg - x_cp)/c) * CL;
    end
end

figure()
hold on
plot(rad2deg(alpha_vec), Cm_cg(1, :), 'LineWidth', 1)
plot(rad2deg(alpha_vec), Cm_cg(2, :), 'LineWidth', 1)
plot(rad2deg(alpha_vec), Cm_cg(3, :), 'LineWidth', 1)
plot(rad2deg(alpha_vec), Cm_cg(4, :), 'LineWidth', 1)
plot(rad2deg(alpha_vec), Cm_cg(5, :), 'LineWidth', 1)
plot(rad2deg(alpha_vec), Cm_cg(6, :), 'LineWidth', 1)
plot(rad2deg(alpha_vec), Cm_cg(7, :), 'LineWidth', 1)
title('Moment Coefficient vs. Angle of Attack')
xlabel('Angle of Attack (deg)')
ylabel('Moment Coefficient')
legend('Mach 2', 'Mach 3', 'Mach 4', 'Mach 5', 'Mach 6', 'Mach 7')
hold off


