%% checkPlanarSlice.m
% Sanity check for sliceSTLPlane.m: slices the vehicle at y = 0 (its
% centerline plane), plots the resulting 2D outline, and overlays each
% panel's outward normal as an arrow. Checking the normals visually
% matters here -- mntCLCD.m's force signs depend entirely on them
% pointing away from the body, not into it.

clear; clc;

stlFile     = 'Bunny_Bomb - Revolve2.stl';  % <-- update path
planeAxis   = 'x';
planeOffset = 0;
rotationDeg = 90;             % CCW rotation applied to the 2D outline ("left"); 0 to disable

%% Slice
[segs, normals] = sliceSTLPlane(stlFile, planeAxis, planeOffset, rotationDeg);
fprintf('%d intersection segments found.\n', size(segs,1));

nAmbiguous = sum(any(isnan(normals),2));
if nAmbiguous > 0
    fprintf('%d of %d panels have an ambiguous normal (near-parallel to cut plane) and will be skipped by mntCLCD.\n', ...
        nAmbiguous, size(segs,1));
end

%% Closure check
% Every point on a clean cut through a solid should be shared by exactly
% 2 segments. Odd multiplicities flag a problem (open/non-watertight
% mesh, plane grazing an edge or vertex, non-manifold geometry, etc.)
allPts = [segs(:,1:2); segs(:,3:4)];
[uPts, ~, ic] = uniquetol(allPts, 1e-7, 'ByRows', true, 'DataScale', 1);
counts = accumarray(ic, 1);
nOdd = sum(mod(counts,2) ~= 0);
if nOdd > 0
    warning(['%d points have odd multiplicity -- outline may not be a clean ' ...
        'closed contour. Check mesh watertightness or try a slightly ' ...
        'different planeOffset.'], nOdd);
else
    fprintf('Closure check passed: all %d unique contour points have even multiplicity.\n', ...
        size(uPts,1));
end

%% Plot: outline + outward normal arrows
figure('Name', 'Planar slice sanity check');
hold on; grid on; axis equal;

X = [segs(:,1), segs(:,3), nan(size(segs,1),1)]';
Z = [segs(:,2), segs(:,4), nan(size(segs,1),1)]';
plot(X(:), Z(:), 'b-', 'LineWidth', 1.2, 'DisplayName', 'Sliced outline');

mid = [(segs(:,1)+segs(:,3))/2, (segs(:,2)+segs(:,4))/2];
valid = ~any(isnan(normals),2);
quiver(mid(valid,1), mid(valid,2), normals(valid,1), normals(valid,2), ...
    0.3, 'r', 'DisplayName', 'Outward normals');

xlabel('c_1 (rotated)'); ylabel('c_2 (rotated)');
title(sprintf('Vehicle cross-section at %s = %.4g, rotated %g%s CCW ("left")', ...
    planeAxis, planeOffset, rotationDeg, char(176)));
legend('Location', 'best');
