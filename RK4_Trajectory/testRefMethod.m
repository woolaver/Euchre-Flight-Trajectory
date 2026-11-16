%test reference method function

clear
clc
close all

vehicleProperties = half_cone_properties();

length = vehicleProperties.length();
height = vehicleProperties.height(); %m
alt = 0; %m
[~, a_inf, ~, rho_inf] = atmoscoesa(alt);
alpha = deg2rad(0);
M_inf = 3;
V_inf = a_inf*M_inf;
Sref = vehicleProperties.Sref;
T_wall = vehicleProperties.Tw; %K, room temperature assumption at launch
T_internal = 293; %K room temperature assumption at launch minus 5 to see conduction temperature change
T_nose = 302; %K room temperature assumption at launch plus 5 to see conduction temperature change

[Ffx, Ffy, Q_conv_rate] = refTempMethod(length, height, alt, M_inf, alpha, T_wall);

CL_vec = readmatrix(vehicleProperties.CL_vec_file);
CD_vec = readmatrix(vehicleProperties.CD_vec_file);

alpha_vec = deg2rad(linspace(-30, 30, 200));
Mach_vec = linspace(1.5, 9, 200);

CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, M_inf, 'linear');
CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, M_inf, 'linear');

Lp = .5*rho_inf*V_inf^2*Sref*CL;
Dp = .5*rho_inf*V_inf^2*Sref*CD;

%test fayRiddell Method

vehicleProperties = half_cone_properties();

Q_stag_rate = fayRiddellMethod(vehicleProperties, M_inf, alt, T_nose);

[T_wall_new, T_internal_new, T_nose_new] = updateTemperatures(vehicleProperties, alt, T_wall, T_internal, T_nose, Q_conv_rate, Q_stag_rate);

disp("Friction Force X Direction: " + Ffx + " N")
disp("Friction Force Y Direction: " + Ffy + " N")
disp("Pressure Force X Direction: " + Dp + " N")
disp("Pressure Force Y Direction: " + Lp + " N")
disp("New Tw: " + T_wall_new + " K")
disp("New T_internal: " + T_internal_new + " K")
disp("New T_nose: " + T_nose_new + " K")