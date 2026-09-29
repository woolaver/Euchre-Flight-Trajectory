%script to estimate CD and CL of vehicle at supersonic speeds > M = 5 using
%shock-expansion theory

clear
clc
close all

%INPUTS:

Mach = 3;
alpha = deg2rad(20);
geometry = 'Bunny_Bomb - Revolve2.stl';

% 1. Read the STL file
STL = stlread(geometry);

% 2. Extract vertices and faces if needed
vertices = STL.Points;
faces = STL.ConnectivityList;

% 3. Visualize the 3D object
figure;
trimesh(STL); % Or use trisurf(TR) for a shaded surface
axis equal;
title('Imported STL Mesh');

%save all vertices on the y-axis in seperate vectors
