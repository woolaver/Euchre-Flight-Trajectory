%driver script for running modified Newtonian function with given .stl file

clear
clc
close all

%Geometry to import, must be .stl file
%geometry = 'csv_files/Bunny_Bomb - Revolve2.stl';
%geometry = 'csv_files/HARV.STL';
%geometry = "csv_files/blend_sharp.stl";
%geometry = "csv_files/Virginia.STL";
%geometry = "csv_files/blend_round.STL";
geometry = "csv_files/blend_sharp.STL";

%helper function to ensure .stl file is correct
%plotSTLRaw(geometry, true)

M_inf = 8;
alpha = deg2rad(25);
gamma = 1.4;

%for cruise take bottom surface from CAD file

%sref = .455; %m^2

% Virginia
sref= .439; %m^2, take flat surface

[CL, CD, CA, CN, CY, Cp_vec, sinTheta, x_cp] = newtonianCLCD3D(geometry, alpha, M_inf, sref,gamma);

angles = rad2deg(sinTheta);
fprintf('Newtonian Cl = %4f, Cd = %4f x_cp percent of length = %4f \n', CL, CD, x_cp);
plotVehicleCp3D(geometry, Cp_vec, sinTheta, false,[])