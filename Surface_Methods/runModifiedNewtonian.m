%driver script for running modified Newtonian function with given .stl file

clear
clc
close all

%Geometry to import, must be .stl file
geometry = 'Bunny_Bomb - Revolve2.stl';
%geometry = 'HARV.STL';

%helper function to ensure .stl file is correct
%plotSTLRaw(geometry, true)

M_inf = 5;
alpha = deg2rad(20);
gamma = 1.4;


%for cruise take bottom surface from CAD file
sref_lift = .42;
sref_drag = .03;

[CL, CD, CA, CN, CY, Cp_vec, sinTheta,x_cp] = newtonianCLCD3D(geometry, alpha, M_inf, sref_lift,sref_drag,gamma);

angles = rad2deg(sinTheta);
fprintf('Newtoian Cl = %4f, Cd = %4f x_cp percent of length = %4f \n',CL, CD,x_cp);
plotVehicleCp3D(geometry, Cp_vec,sinTheta , false,[])