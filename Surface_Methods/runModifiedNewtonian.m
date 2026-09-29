%driver script for running modified Newtonian function with given .stl file

clear
clc
close all

%Geometry to import, must be .stl file
geometry = 'Bunny_Bomb - Revolve2.stl';
%geometry = 'HARV.STL';

%helper function to ensure .stl file is correct
%plotSTLRaw(geometry, true)

M_inf = 8;
alpha = deg2rad(40);

gamma = 1.4;

p_pinf = (((gamma+1)^2*M_inf^2)/((4*gamma*M_inf^2) - 2*(gamma - 1)))^(gamma/(gamma-1))*((1-gamma+(2*gamma*M_inf^2))/(gamma+1)); 
Cp_max = 2/(gamma*M_inf^2)*(p_pinf - 1);

%for cruise take bottom surface from CAD file
S_ref = .47;

[CL, CD, CA, CN, CY, Cp_vec, sinTheta] = newtonianCLCD3D(geometry, alpha, Cp_max, S_ref);

angles = rad2deg(sinTheta);

plotVehicleCp3D(geometry, Cp_vec, false)