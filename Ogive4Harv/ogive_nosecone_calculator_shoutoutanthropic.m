%% ====================================================================
%  Von Karman Ogive Nosecone Generator (with blunted spherical tip)
%  ====================================================================
%  Generates an X,Y,Z point cloud (Z = 0, planar profile) for a von
%  Karman (LD-Haack, C = 0) ogive nosecone whose sharp mathematical tip
%  is replaced with a tangent spherical cap of specified radius.
%
%  The profile is exported as a CSV that can be imported into SolidWorks
%  via Insert > Curve > Curve Through XYZ Points, sketched, and revolved
%  about the X axis to produce the solid nosecone.
%
%  IMPORTANT: The overall length you specify (L_target) is treated as
%  the FINAL, BLUNTED overall length (tip of sphere to base). Because
%  blunting truncates a small amount of material off the sharp tip, the
%  script solves iteratively for the underlying (longer) sharp-ogive
%  length so that, after blunting, the finished part is exactly
%  L_target inches long.
%
%  Units: inches (change R, L_target, r_n as needed; units are
%  otherwise unitless/consistent).
% =====================================================================

clear; clc; close all;

%% -------- USER INPUTS --------
R          = 0.5*7;              % Base (max) radius of nosecone      [in]
L_target   = 5.1913*7;           % Desired FINAL overall length        [in]
r_n        = 0.05*7;             % Spherical tip blunt radius          [in]
N_ogive    = 5000;              % # points along ogive curve
N_sphere   = 5000;               % # points along spherical cap
outputFile     = 'vonKarman_nosecone.csv';
outputFileTxt  = 'vonKarman_nosecone.txt';

%% -------- SHAPE FUNCTIONS --------
y_of_x  = @(x,L) vk_y(x,L,R);
dy_of_x = @(x,L) vk_dy(x,L,R);

%% -------- SOLVE FOR UNDERLYING SHARP-OGIVE LENGTH --------
% blunted_OAL(L) = L - x_trunc(L). We need blunted_OAL(L) == L_target.
oal_error = @(L) blunted_OAL(L,r_n,y_of_x,dy_of_x) - L_target;

L_lo = L_target;
L_hi = L_target*1.05;
iter = 0;
while oal_error(L_lo)*oal_error(L_hi) > 0
    L_hi = L_hi*1.05;
    iter = iter+1;
    if iter > 50
        error('Could not bracket a solution for ogive length. Check inputs.');
    end
end
L_ogive = fzero(oal_error, [L_lo, L_hi]);

% Recover full geometry at the solution
[x_trunc, x0, y0, xc] = blunted_geom(L_ogive, r_n, y_of_x, dy_of_x);

fprintf('---------------------------------------------------\n');
fprintf('Solved sharp-ogive length : %.6f in\n', L_ogive);
fprintf('Tangency point (x0)       : %.6f in\n', x0);
fprintf('Tangency radius (y0)      : %.6f in\n', y0);
fprintf('Sphere center (xc)        : %.6f in\n', xc);
fprintf('Truncated tip length      : %.6f in\n', x_trunc);
fprintf('Final blunted OAL         : %.6f in\n', L_ogive - x_trunc);
fprintf('Base radius (check)       : %.6f in\n', y_of_x(L_ogive, L_ogive));
fprintf('---------------------------------------------------\n');

%% -------- BUILD POINT CLOUD --------
% Spherical cap: parametrize by angle phi about sphere center (xc,0)
%   x = xc + r_n*cos(phi),  y = r_n*sin(phi)
% Forward-most point (nose tip) is at phi = pi -> x = xc - r_n, y = 0
% Tangency point corresponds to phi = phi0
phi0 = atan2(y0, x0 - xc);
phi  = linspace(pi, phi0, N_sphere)';
x_sph = xc + r_n*cos(phi);
y_sph = r_n*sin(phi);

% Ogive curve: from tangency point x0 out to the base at L_ogive
x_og = linspace(x0, L_ogive, N_ogive)';
y_og = y_of_x(x_og, L_ogive);

% Combine (drop duplicate point at the sphere/ogive junction)
x_all = [x_sph; x_og(2:end)];
y_all = [y_sph; y_og(2:end)];

% Shift so the nose tip sits at X = 0
x_all = x_all - (xc - r_n);
z_all = zeros(size(x_all));

%% -------- WRITE CSV & TXT (no header -- plain numeric XYZ for SolidWorks) --------
M = [x_all, y_all, z_all];
writematrix(M, outputFile);
writematrix(M, outputFileTxt);   % identical data, .txt extension (also SW-compatible)

fprintf('Wrote %d points to %s and %s\n', size(M,1), outputFile, outputFileTxt);
fprintf('Nose tip at X = 0, base at X = %.4f, base radius = %.4f\n', ...
        max(x_all), R);

%% -------- PLOT (sanity check) --------
figure; hold on; grid on; axis equal;
plot(x_all, y_all, 'b-', 'LineWidth', 1.5);
plot(x_all, -y_all, 'b-', 'LineWidth', 1.5); % mirror for visualization only
xlabel('X (in)'); ylabel('Y (in)');
title('Von Karman Ogive Nosecone Profile (Blunted Tip)');

%% ===================== LOCAL FUNCTIONS =====================

function y = vk_y(x, L, R)
% Von Karman (Haack series, C = 0) ogive radius as a function of
% axial position x, for a curve of length L and base radius R.
    s = 1 - 2.*x./L;
    s = min(max(s,-1),1);            % guard against roundoff at endpoints
    theta = acos(s);
    val = theta - sin(2*theta)/2;
    val(val<0) = 0;                  % guard against tiny negative roundoff
    y = (R/sqrt(pi)) .* sqrt(val);
end

function dy = vk_dy(x, L, R)
% Analytical derivative dy/dx of the von Karman ogive.
    s = 1 - 2.*x./L;
    s = min(max(s,-1),1);
    theta = acos(s);
    sinth = sin(theta);
    dy = zeros(size(x));
    for i = 1:numel(x)
        if sinth(i) < 1e-9
            dy(i) = 1e6;              % ~vertical tangent at the sharp tip
        else
            dtheta_dx = 2/(L*sinth(i));
            denom = sqrt(max(theta(i) - sin(2*theta(i))/2, 0));
            if denom < 1e-9
                dydtheta = 0;
            else
                dydtheta = (R/sqrt(pi)) * sinth(i)^2/denom;
            end
            dy(i) = dydtheta*dtheta_dx;
        end
    end
end

function [x_trunc, x0, y0, xc] = blunted_geom(L, r_n, y_of_x, dy_of_x)
% For a sharp ogive of length L, find the point (x0,y0) where a sphere
% of radius r_n is tangent to the curve, the resulting sphere center
% xc (on axis), and the length x_trunc truncated off the original tip.
    g = @(x) y_of_x(x,L).*sqrt(1+dy_of_x(x,L).^2) - r_n;

    x_lo = 1e-8*L;
    x_hi = 0.3*L;
    grow = 0;
    while g(x_lo)*g(x_hi) > 0 && grow < 20
        x_hi = min(x_hi*1.5, 0.9*L);
        grow = grow+1;
    end
    x0 = fzero(g, [x_lo, x_hi]);

    y0 = y_of_x(x0,L);
    m  = dy_of_x(x0,L);
    xc = x0 + y0*m;
    x_trunc = xc - r_n;
end

function oal = blunted_OAL(L, r_n, y_of_x, dy_of_x)
% Final overall length after blunting, for a given sharp-ogive length L.
    [x_trunc,~,~,~] = blunted_geom(L, r_n, y_of_x, dy_of_x);
    oal = L - x_trunc;
end