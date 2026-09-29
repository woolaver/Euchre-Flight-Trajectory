function plotSTLRaw(stlFile, showAsSurface)
% PLOTSTLRAW Quick 3D plot of an STL file's raw geometry, to sanity-check
% that the mesh looks like the vehicle you expect before running any
% slicing or Newtonian-theory code on it.
%
% plotSTLRaw(stlFile, showAsSurface)
%
% Inputs:
%   stlFile       - path to .stl file
%   showAsSurface - true (default): render the triangulated surface
%                   (patch), which is usually the more useful check
%                   since it shows the actual faceted shape, not just a
%                   cloud of points. false: plot raw vertices only, as a
%                   3D scatter -- useful if you specifically want to see
%                   mesh density/vertex distribution, or if the surface
%                   render looks wrong and you want to check whether the
%                   underlying vertices themselves are sane.
%
% Prints basic diagnostics (vertex count, face count, bounding box) so
% you can catch obviously wrong imports (e.g. units off by 1000x,
% empty/degenerate mesh) without even looking at the plot.

if nargin < 2 || isempty(showAsSurface), showAsSurface = true; end

TR = stlread(stlFile);
V = TR.Points;
F = TR.ConnectivityList;

fprintf('%s: %d vertices, %d faces\n', stlFile, size(V,1), size(F,1));
fprintf('  x range: [%.4g, %.4g]\n', min(V(:,1)), max(V(:,1)));
fprintf('  y range: [%.4g, %.4g]\n', min(V(:,2)), max(V(:,2)));
fprintf('  z range: [%.4g, %.4g]\n', min(V(:,3)), max(V(:,3)));

figure('Name', 'Raw STL check');

if showAsSurface
    patch('Faces', F, 'Vertices', V, ...
        'FaceColor', [0.8 0.8 0.9], 'EdgeColor', [0.3 0.3 0.3], ...
        'FaceAlpha', 1.0);
else
    plot3(V(:,1), V(:,2), V(:,3), 'k.', 'MarkerSize', 2);
end

axis equal; grid on; view(3);
camlight('headlight'); lighting gouraud;
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('%s (%d verts, %d faces)', stlFile, size(V,1), size(F,1)), ...
    'Interpreter', 'none');

end
