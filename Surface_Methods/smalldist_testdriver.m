%test driver for small dist theory
clc
clear
close all

minf = 5;
alpha = deg2rad(20);
sref_lift = .42;
sref_drag = .03;
y = 1.4;
stlfile = 'Bunny_Bomb - Revolve2.stl';
showVectors = false;

mu = rad2deg(asin(1/minf));
[CL, CD, CA, CN, CY, Cp, theta, centroid, x_cp] = smalldist_supersonic(stlfile, alpha, minf, sref_lift, sref_drag,y);
fprintf('Smal DistCl = %4f, Cd = %4f x_cp percent of length = %4f \n',CL, CD,x_cp);
%plotVehicleCp3D(stlfile, Cp, rad2deg(theta) , showVectors, [])