%% Example driver for hypersonicHeatFlux.m
% Fill in your actual trajectory data and geometry/material properties.

clear; clc; close all;

% --- Trajectory data (replace with real data) ---

traj.time     = linspace(0, 80, 800);              % s
traj.altitude = linspace(0000, 0000, 800);         % m, ascending
traj.velocity = linspace(343*8, 343*3, 800);           % m/s, decelerating
traj.alpha    = linspace(4, 4, 800);           % deg, example AoA history

trajectory_data = readmatrix("Heating Trajectory.csv");
trajectory_data = trajectory_data(:,1:(length(trajectory_data)-1));
traj.time = trajectory_data(1, :);
traj.velocity = trajectory_data(2,:);
traj.alpha = rad2deg(trajectory_data(6,:));
traj.altitude = trajectory_data(4,:);

%[~, max_alt_i] = max(traj.altitude);
%i_dive = 700;
%traj.alpha = [zeros(size(traj.alpha))];
%traj.alpha(1, max_alt_i:i_dive) = 15;



% --- Geometry (placeholder values -- update with actual vehicle data) ---
geom.x          = 0.1;                  % m, station for the time-history run
geom.xStations  = linspace(0.1, 1, 100); % m, stations for the position-profile run
geom.halfAngle  = 8;                    % deg, body half-angle (set to [] to skip shock calc)
geom.emissivity = 0.85;                 % wall emissivity

% --- Wall/material properties (placeholder values) ---
ss.density   = 8000;     % kg/m^3
ss.cp        = 480;      % J/kg-K (Aluminum = 900)

aluminum.density   = 2700;     % kg/m^3
aluminum.cp        = 900;      % J/kg-K (Aluminum = 900)

mat = aluminum;
mat.thickness = 0.04;    % m
mat.Tw0       = 288;      % K, initial wall temperature

% --- Model options ---
opts.Re_trans         = 5e5;
opts.includeRadiation = true;
opts.makePlots        = true;

%% Run 1: quantities vs. TIME at the single station geom.x
opts.mode = 'time';
resultsTime = hypersonicHeatFlux(traj, geom, mat, opts);

%% Run 2: quantities vs. POSITION (0.1 to 1 m) at a single snapshot time
opts.mode = 'position';
opts.snapshotTime = 10;   % s, must lie within traj.time
resultsPos = hypersonicHeatFlux(traj, geom, mat, opts);
