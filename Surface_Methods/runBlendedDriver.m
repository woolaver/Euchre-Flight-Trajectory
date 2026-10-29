%driver script for running blended driver with given .stl file

clear
%clc
close all

%Geometry to import, must be .stl file
%geometry = 'csv_files/Bunny_Bomb - Revolve2.stl';
%geometry = 'csv_files/HARV.STL';
%geometry = "csv_files/blend_sharp.stl";
geometry = 'csv_files/Virginia.STL';
%geometry = "csv_files/blend_round.STL";
%geometry = "csv_files/blend_flat.STL";

%helper function to ensure .stl file is correct
%plotSTLRaw(geometry, true)

M_inf = 5;
alpha = deg2rad(25);
gamma = 1.4;

%for cruise take bottom surface from CAD file
%sref = .445; %m^2

sref = .637; %m^2


% Virginia
%sref = .439; %m^2, take flat surface

[CL, CD, Cp_vec, x_cp] = blended_driver(geometry, alpha, M_inf, sref, gamma);

fprintf('Blended Cl = %4f, Cd = %4f x_cp percent of length = %4f \n',CL, CD,x_cp);
%plotVehicleCp3D(geometry, Cp_vec, false,[])