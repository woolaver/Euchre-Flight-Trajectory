function plotVehicleCp3D(stlFile, Cp, showVectors, vectorScale)
% PLOTVEHICLECP3D Plot the full 3D STL surface colored by Cp, with an
% optional overlay of Cp*normal vectors (magnitude = Cp, direction =
% outward normal) at a decimated set of face centroids.
%
% plotVehicleCp3D(stlFile, Cp, showVectors, vectorScale)
%
% Inputs:
%   stlFile     - path to .stl file (same one passed to newtonianCLCD3D)
%   Cp          - Nfaces x 1 vector from newtonianCLCD3D, in the same
%                 face order as the STL's ConnectivityList
%   showVectors - true to overlay Cp*normal arrows (default false --
%                 for a mesh with thousands of faces, the colored
%                 surface alone is usually the more readable view;
%                 arrows are decimated to ~300 for legibility when on)
%   vectorScale - arrow-length multiplier, auto-picked if omitted
%
% The colored surface is the primary view for a full 3D mesh (direction
% is implicit -- arrows point along the local outward normal -- so color
% carries the magnitude, exactly analogous to the 2D Cp*normal plot but
% adapted for a mesh too dense for per-face arrows to stay readable).

if nargin < 3 || isempty(showVectors), showVectors = false; end

TR = stlread(stlFile);
V = TR.Points;
F = TR.ConnectivityList;

figure('Name', '3D vehicle surface Cp');
patch('Faces', F, 'Vertices', V, 'FaceVertexCData', Cp, ...
    'FaceColor', 'flat', 'EdgeColor', 'none');
colormap(jet);
cb = colorbar;
cb.Label.String = 'C_p';
axis equal; view(3); grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title('3D Newtonian Cp distribution');

if showVectors
    n = faceNormal(TR);
    v1 = V(F(:,1),:); v2 = V(F(:,2),:); v3 = V(F(:,3),:);
    centroids = (v1 + v2 + v3) / 3;

    if nargin < 4 || isempty(vectorScale)
        charLen = max(V(:,1)) - min(V(:,1));
        maxCp = max(abs(Cp));
        if maxCp == 0, maxCp = 1; end
        vectorScale = 0.1 * charLen / maxCp;
    end

    nFaces = size(F,1);
    step = max(1, floor(nFaces/300));   % cap at ~300 arrows for legibility
    idx = 1:step:nFaces;
    vec = (Cp(idx) .* n(idx,:)) * vectorScale;

    hold on;
    quiver3(centroids(idx,1), centroids(idx,2), centroids(idx,3), ...
        vec(:,1), vec(:,2), vec(:,3), 0, 'k', 'LineWidth', 1);
end

end