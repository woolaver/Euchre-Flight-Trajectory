%test reference method function

clear
clc
close all

%blend_sharp geometry
length = 1.491; %m
height = .126; %m
alt = 1; %m
alpha = deg2rad(0);
M_inf = 8;
T_wall = 298; %K, room temperature assumption at launch
A = 1.679; %m^2
cp_Ti = 523; %J/kg*K for TI6ALV4
cp_Al = 870; %J/kg*K for 7075-T6 Al
m_skin = 60; %kg, estimate
m_internal = 120 - m_skin;
dt = .01;
sigma = 5.67*10^(-8); %W/m^2*K
emissivity = .4; %estimate for TI6ALV4
wall_thickness = .0254; %m, assuming 1 inch thickness for testing purpose, need structure for real value
k_Ti = 7; %W/m*K, estimate for TI6ALV4
k_Al = 130; %W/m*K
T_internal = 293; %K room temperature assumption at launch minus 5 to see conduction temperature change

[Ffx, Ffy, Q_conv_rate] = refTempMethod(length, height, alt, M_inf, alpha, T_wall);

[T_env, ~, ~, ~] = atmoscoesa(alt);

%radiative heat transfer
Q_rad_rate = emissivity*sigma*A*(T_env);

%conductive heat transfer using steady state assumption for each time point
Q_cond_rate = (k_Ti*A*(T_wall - T_internal))/(wall_thickness);

Q_wall_total_rate = Q_rad_rate + Q_conv_rate - Q_cond_rate;

Q_wall_total = Q_wall_total_rate * dt;
Q_internal_total = Q_cond_rate*dt;

delta_T_wall = Q_wall_total/(m_skin * cp_Ti);
delta_T_internal = Q_internal_total/(m_internal * cp_Al);

T_wall_new = T_wall + delta_T_wall;
T_internal_new = T_internal + delta_T_internal;

disp("Friction Force X Direction: " + Ffx + " N")
disp("Friction Force Y Direction: " + Ffy + " N")
disp("New Tw: " + T_wall_new + " K")
disp("New T_internal: " + T_internal_new + " K")