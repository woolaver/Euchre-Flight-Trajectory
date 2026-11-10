%test script to test calcFinForces function

clear
clc
close all

vehicleProperties = half_cone_properties();

M_inf = 3;
alt = 3000;
[~, ~, p_inf, ~] = atmoscoesa(alt);

theta_OS = tan(vehicleProperties.fin_y_location/vehicleProperties.fin_x_location);

[ps_pinf, Ms, ~] = tbm(M_inf, theta_OS);
ps = ps_pinf*p_inf;

alpha_fin = deg2rad(20);

[CL_fin, CD_fin] = calcFinForces(vehicleProperties, ps, Ms, alpha_fin, M_inf, p_inf);

disp("CL_fin: " + CL_fin)
disp("CD_fin: " + CD_fin)