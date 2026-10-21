%sweep over range of Mach numbers and angles of attack and write data to a
%data file

clear
clc
close all


%stlfile = 'csv_files\Virginia.stl';
%stlfile = 'csv_files\HARV.stl';
stlfile = 'csv_files\Ice_Cream_Cone_STL.STL';
%stlfile = 'csv_files\Bunny_Bomb - Revolve2.stl';

%helper function to ensure .stl file is correct
%plotSTLRaw(geometry, true)

M_inf = linspace(2, 8, 200);
alpha = deg2rad(linspace(-30, 30, 200));

CL_vec = zeros([200, 200]);
CD_vec = zeros([200, 200]);

gamma = 1.4;

% Virginia
%{
sref_lift = .439; %m^2, take flat surface
sref_drag = .035; %m^2 take biggest cross sectional area from front face (i.e. back surface)
%}

%HARV
%{
sref_lift = .359; %m^2
sref_drag = .029; %m^2
%}

%Ice Cream Cone

sref_lift = .253; %m^2
sref_drag = .0299; %m^2 
%Note: Took the diameter of the thickest part and
%assumed it is a circle then found area, not exactly the thickets cross
%sectional area but pretty close for sref_drag


%Bunny Bomb
%{
sref_lift = .467; %m^2
sref_drag = .049; %m^2
%}

for i = 1:200
    disp(i)
    for j = 1:200
        [Cl_s, Cd_s, ~, ~, ~, ~, ~, ~, ~] = smalldist_supersonic(stlfile, alpha(j), M_inf(i), sref_lift, sref_drag, gamma);
        [Cl_n, Cd_n, ~, ~, ~, ~, ~, ~] = newtonianCLCD3D(stlfile, alpha(j), M_inf(i), sref_lift, sref_drag, gamma);

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

datafile_CL = "Ice_Cream_Cone_blended_CL.csv";
writematrix(CL_vec, datafile_CL)

datafile_CD = "Ice_Cream_Cone_blended_CD.csv";
writematrix(CD_vec, datafile_CD)
