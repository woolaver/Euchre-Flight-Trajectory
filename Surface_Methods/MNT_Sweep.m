%sweep over range of Mach numbers and angles of attack and write data to a
%data file

clear
clc
close all

%Geometry to import, must be .stl file
%geometry = 'Bunny_Bomb - Revolve2.stl';
geometry = 'HARV.STL';  %NOTE if using HARV you have to divide points
%geometry = 'Ice_Cream_Cone_STL.STL';
%geometry = 'Virginia.stl';

%helper function to ensure .stl file is correct
%plotSTLRaw(geometry, true)

M_inf = linspace(2, 8, 50);
alpha = deg2rad(linspace(-30, 30, 500));

CL_vec = zeros([50, 500]);
CD_vec = zeros([50, 500]);

gamma = 1.4;

for i = 1:50
    disp(i)
    for j = 1:500
        p_pinf = (((gamma+1)^2*M_inf(i)^2)/((4*gamma*M_inf(i)^2) - 2*(gamma - 1)))^(gamma/(gamma-1))*((1-gamma+(2*gamma*M_inf(i)^2))/(gamma+1)); 
        Cp_max = 2/(gamma*M_inf(i)^2)*(p_pinf - 1);
        
        %for cruise take bottom surface from CAD file
        %S_ref = .029; %For HARV
        %S_ref = .049; %For Bunny_bomb
        S_ref = (.2^2)*pi;
        
        [CL, CD, CA, CN, CY, Cp_vec, sinTheta] = newtonianCLCD3D(geometry, alpha(j), Cp_max, S_ref);
        CL_vec(i, j) = CL;
        CD_vec(i, j) = CD;
    end
end

datafile_CL = "HARV_CL_50_500.csv";
writematrix(CL_vec, datafile_CL)

datafile_CD = "HARV_CD_50_500.csv";
writematrix(CD_vec, datafile_CD)
