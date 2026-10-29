%test driver for small dist theory
%clc
%clear
%close all

minf = 2;
alpha = deg2rad(-20);
sref = .445; %m^2
y = 1.4;
stlfile = 'csv_files/Virginia.STL';
showVectors = false;

mu = rad2deg(asin(1/minf));
[CL, CD, CA, CN, CY, Cp, theta, centroid, x_cp] = smalldist_supersonic(stlfile, alpha, minf, sref,y);
fprintf('Smal DistCl = %4f, Cd = %4f x_cp percent of length = %4f \n',CL, CD,x_cp);
plotVehicleCp3D(stlfile, Cp, rad2deg(theta) , showVectors, [])