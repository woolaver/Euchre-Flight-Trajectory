function plotVehicleCp(segs, normals, Cp, scaleFactor)
% PLOTVEHICLECP Plot the sliced 2D vehicle outline with the local
% pressure coefficient shown directly on the surface as a vector: each
% arrow's direction is the panel's outward normal, and its length is
% proportional to Cp at that panel (Cp*normal, elementwise). Positive Cp
% (windward, pushing outward) and non-positive Cp (shadowed/leeward, per
% the basic Newtonian shadowing rule) are colored separately so the
% windward/leeward split is visible at a glance.
%
% plotVehicleCp(segs, normals, Cp, scaleFactor)
%
% Inputs:
%   segs        - Ns x 4 matrix from sliceSTLPlane: [c1a c2a c1b c2b]
%   normals     - Ns x 2 matrix from sliceSTLPlane: outward unit normal
%                 per panel (NaN rows are skipped)
%   Cp          - Ns x 1 vector from mntCLCD (or any other per-panel Cp)
%   scaleFactor - arrow-length multiplier. If omitted, it's chosen
%                 automatically so the largest arrow spans roughly 15%
%                 of the vehicle's extent along c1 -- override this if
%                 the arrows look too small/large for your geometry.
%
% This only draws the geometry and pressure field -- it doesn't compute
% anything new. Run sliceSTLPlane and mntCLCD first and pass their
% outputs straight in.

if nargin < 4 || isempty(scaleFactor)
    charLen = max([segs(:,1); segs(:,3)]) - min([segs(:,1); segs(:,3)]);
    maxCp = max(abs(Cp(~isnan(Cp))));
    if isempty(maxCp) || maxCp == 0
        maxCp = 1;
    end
    scaleFactor = 0.15 * charLen / maxCp;
end

mid = [(segs(:,1)+segs(:,3))/2, (segs(:,2)+segs(:,4))/2];
valid = ~any(isnan(normals),2) & ~isnan(Cp);

vec = (Cp .* normals) * scaleFactor;   % elementwise: magnitude=Cp, direction=normal

posMask    = valid & Cp > 0;
nonPosMask = valid & Cp <= 0;

figure('Name', 'Vehicle surface Cp');
hold on; grid on; axis equal;

% Outline
X = [segs(:,1), segs(:,3), nan(size(segs,1),1)]';
Z = [segs(:,2), segs(:,4), nan(size(segs,1),1)]';
plot(X(:), Z(:), 'k-', 'LineWidth', 1.2, 'DisplayName', 'Vehicle outline');

% Windward panels: Cp*normal arrows, pointing outward, length ~ Cp
quiver(mid(posMask,1), mid(posMask,2), vec(posMask,1), vec(posMask,2), 0, ...
    'r', 'LineWidth', 1, 'DisplayName', 'Cp > 0 (windward)');

% Shadowed/non-positive panels: usually ~zero length under basic
% Newtonian shadowing, so mark the panel instead of drawing an
% invisible arrow -- still useful to see which panels are shadowed
plot(mid(nonPosMask,1), mid(nonPosMask,2), 'b.', 'MarkerSize', 6, ...
    'DisplayName', 'Cp \leq 0 (shadowed)');

xlabel('c_1'); ylabel('c_2');
title(sprintf('Surface Cp distribution (arrow length = Cp \\times %.4g)', scaleFactor));
legend('Location', 'best');

end