%% Example driver for fayRiddellStagnation.m
% Stagnation-point heat flux on a blunted nose behind the bow shock.
% All values below are placeholders -- replace with your own.

clear; clc;

% --- Freestream (altitude is passed through stdatmfull) ---
fs.altitude = 0000;     % m
fs.velocity = 2720;      % m/s
% To use a trajectory instead, give vectors of equal length, e.g.:
%   fs.time = traj.time;  fs.altitude = traj.altitude;  fs.velocity = traj.velocity;
% To bypass stdatmfull, omit fs.altitude and set fs.T [K] and fs.p [Pa].


trajectory_data = readmatrix("Trajectory.csv");
trajectory_data = trajectory_data(:,1:(length(trajectory_data)-1));
fs.time = trajectory_data(1, :);
fs.velocity = trajectory_data(2,:);
traj.altitude = trajectory_data(4,:);
%}


% --- Geometry ---
geom.Rn = 0.05;          % m, nose radius

% --- Wall ---
wall.Tw = 288;           % K, isothermal wall temperature

% --- Model options ---
opts.gasModel      = 'equilibrium';   % 'equilibrium', 'vibEq' or 'perfect'
opts.Pr            = 0.71;
opts.Le            = 1.4;
opts.boundaryLayer = 'equilibrium';   % 'equilibrium' (a = 0.52) or 'frozen' (a = 0.63)
opts.hD            = [];              % J/kg, [] = auto from equilibrium composition (0 for vibEq/perfect)
opts.makePlots     = true;

% --- Run ---
results = fayRiddellStagnation(fs, geom, wall, opts);

%{
Rns = linspace(0.01,0.2,100);
qs = size(Rns);
for i = 1:(length(Rns))
    geom.Rn = Rns(i);
    results = fayRiddellStagnation(fs, geom, wall, opts);
    qs(i) = results.q_w;
end
figure();
plot(Rns,qs/1e4, LineWidth=2)
xlabel('Nose Radius [m]');
ylabel('Fay Riddell Heat Flux [W/cm^2]');
title(['Heat flux on Blunted Nose (M_inf = 8, Sea level)'])
%}