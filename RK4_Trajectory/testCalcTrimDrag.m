%test script to run calcTrimDrag

tic

clear
clc
close all

vehicleProperties = half_cone_properties();
M_inf = 8;
alpha_body = vehicleProperties.alpha_climb;
alt = 0;

[~, a_inf, ~, rho_inf] = atmoscoesa(alt);
V_inf = M_inf*a_inf;

Sref = vehicleProperties.Sref;
length = vehicleProperties.length;
height = vehicleProperties.height;

Tw = vehicleProperties.Tw;
T_nose = vehicleProperties.T_nose;

CL_vec = readmatrix(vehicleProperties.CL_vec_file);
CD_vec = readmatrix(vehicleProperties.CD_vec_file);

alpha_vec = deg2rad(linspace(-30, 30, 200));

Mach_vec = linspace(1.5, 9, 200); 
CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha_body, M_inf, 'linear');
CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha_body, M_inf, 'linear');

D = 1/2*rho_inf*V_inf^2*Sref*CD;
L = 1/2*rho_inf*V_inf^2*Sref*CL;

[Fx, Fy, Q_conv_rate] = refTempMethod(length, height, alt, M_inf, alpha_body, Tw);
Q_conv_nose_rate = fayRiddellMethod(vehicleProperties, M_inf, alt, T_nose);

D = D + Fx;
L = L + Fy;

[trimDrag, trimLift, alpha_fin] = calcTrimDrag(vehicleProperties, alpha_body, L, D, alt, M_inf);

total_time = toc;

disp("Trim Drag: " + trimDrag + " N")
disp("Trim Lift: " + trimLift + " N")
disp("Fin Angle of Attack: " + rad2deg(alpha_fin) + " deg")
disp("Time to compute: " + total_time + " s")

total_lift = L + trimLift;

disp("Vehicle G Force: " + total_lift/(120*9.8) + " g's")