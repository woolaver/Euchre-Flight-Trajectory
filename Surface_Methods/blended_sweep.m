%sweep over range of Mach numbers and angles of attack and write data to a
%data file

clear
clc
close all

%stlfile = 'csv_files/Virginia.STL';
%stlfile = 'csv_files/HARV.STL';
%stlfile = 'csv_files/Ice_Cream_Cone_STL.STL';
stlfile = 'csv_files/Bunny_Bomb - Revolve2.STL';
%stlfile = 'csv_files/blend_sharp.STL';
%stlfile = 'csv_files/blend_round.STL';
%stlfile = 'csv_files/blend_flat.STL';

%helper function to ensure .stl file is correct
%plotSTLRaw(geometry, true)

M_inf = linspace(2, 8, 200);
alpha = deg2rad(linspace(-30, 30, 200));

CL_vec = zeros([200, 200]);
CD_vec = zeros([200, 200]);

gamma = 1.4;

% Virginia
%sref = .439; %m^2, take flat surface

%HARV
%sref = .359; %m^2

%Ice Cream Cone
%sref = .253; %m^2

%Bunny Bomb
sref = .467; %m^2

%Blend Sharp/Round
%sref = .445; %m^2

%Blend Flat
%sref = .637; %m^2


for i = 1:200
    disp(i)
    for j = 1:200
        [Cl_s, Cd_s, ~, ~, ~, ~, ~, ~, ~] = smalldist_supersonic(stlfile, alpha(j), M_inf(i), sref, gamma);
        [Cl_n, Cd_n, ~, ~, ~, ~, ~, ~] = newtonianCLCD3D(stlfile, alpha(j), M_inf(i), sref, gamma);

        %interpolation
        M1 = 3;   % below small dist
        M2 = 5;   % above newtonian

        t = (M_inf(i) - M1) / (M2 - M1);
        t = min(max(t, 0), 1);          % clamp

        w = 3*t.^2 - 2*t.^3;             % smoothstep:

        Cl = (1 - w).*Cl_s + w.*Cl_n;
        Cd = (1 - w).*Cd_s + w.*Cd_n;

        CL_vec(i, j) = Cl;
        CD_vec(i, j) = Cd;
    end
end

datafile_CL = "Bunny_Bomb_blended_CL.csv";
writematrix(CL_vec, datafile_CL)

datafile_CD = "Bunny_Bomb_blended_CD.csv";
writematrix(CD_vec, datafile_CD)
