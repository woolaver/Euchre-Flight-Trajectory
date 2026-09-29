function [segs, normals] = sliceSTLPlane(stlFile, planeAxis, planeOffset, rotationDeg)
% SLICESTLPLANE Intersect an STL mesh with a plane perpendicular to a
% global axis, returning the intersection line segments and each
% segment's outward surface normal (projected into the cut plane).
%
% [segs, normals] = sliceSTLPlane(stlFile, planeAxis, planeOffset, rotationDeg)
%
% Inputs:
%   stlFile     - path to .stl file
%   planeAxis   - 'x','y', or 'z': axis normal to the cutting plane
%                 (default 'y' -- i.e. slice at y = planeOffset, giving
%                 the vehicle's profile in the x-z plane)
%   planeOffset - offset of the plane along planeAxis (default 0)
%   rotationDeg - CCW rotation (deg) applied to the returned 2D geometry
%                 about the origin, e.g. 90 moves the point that was
%                 "up" (+c2) to the left (-c1), i.e. "rotate left"
%                 (default 0, no rotation). Applied to both segs and
%                 normals (normals are direction vectors, so they rotate
%                 cleanly with no translation issue).
%
% Outputs:
%   segs    - Ns x 4 matrix, each row = [c1a c2a c1b c2b], the two
%             endpoints of one intersection segment in the 2D plane
%             coordinates (c1, c2 are the two axes NOT equal to
%             planeAxis, kept in their original x/y/z order before
%             rotation)
%   normals - Ns x 2 matrix, each row = the STL face's outward unit
%             normal for that segment's source triangle, projected into
%             the (c1,c2) plane and re-normalized. Rows are NaN where
%             the face normal is nearly perpendicular to the cut plane
%             (i.e. the true 3-D normal points almost entirely along
%             planeAxis) -- those panels don't have a meaningful 2-D
%             outward direction and should be skipped or handled
%             separately downstream (e.g. in an MNT sweep).
%
% This does NOT stitch the segments into an ordered polyline -- each
% triangle-plane intersection is returned as its own segment. That's
% fine (in fact convenient) for panel-sum methods like Modified
% Newtonian Theory, which don't care what order the panels are summed
% in, only their length, position, and outward normal.
%
% IMPORTANT CAVEAT: this is a planar (2-D strip) cut through a 3-D body.
% Integrating pressure over these segments gives a 2-D/per-unit-span
% result, not the true 3-D surface integral -- it ignores any
% out-of-plane curvature or relieving effects. Treat CL/CD computed this
% way as an engineering approximation, not a substitute for a full 3-D
% panel method.

if nargin < 2 || isempty(planeAxis), planeAxis = 'y'; end
if nargin < 3 || isempty(planeOffset), planeOffset = 0; end
if nargin < 4 || isempty(rotationDeg), rotationDeg = 0; end

TR = stlread(stlFile);
V = TR.Points;
F = TR.ConnectivityList;

switch lower(planeAxis)
    case 'x', normalCol = 1; cols = [2 3];
    case 'y', normalCol = 2; cols = [1 3];
    case 'z', normalCol = 3; cols = [1 2];
    otherwise, error('planeAxis must be ''x'', ''y'', or ''z''');
end

d = V(:, normalCol) - planeOffset;

segs    = nan(size(F,1), 4);
normals = nan(size(F,1), 2);
nSeg = 0;

for i = 1:size(F,1)
    tri = F(i,:);
    dv = d(tri);
    s = sign(dv);

    if all(s > 0) || all(s < 0)
        continue % triangle entirely on one side, no crossing
    end

    crossPts = [];
    for e = 1:3
        a = tri(e);
        b = tri(mod(e,3)+1);
        da = d(a); db = d(b);
        if da == 0
            crossPts(end+1,:) = V(a, cols); %#ok<AGROW>
        elseif sign(da) ~= sign(db)
            t = da / (da - db);
            p = V(a,:) + t*(V(b,:) - V(a,:));
            crossPts(end+1,:) = p(cols); %#ok<AGROW>
        end
    end

    crossPts = unique(crossPts, 'rows');
    if size(crossPts,1) == 2
        nSeg = nSeg + 1;
        segs(nSeg,:) = [crossPts(1,:), crossPts(2,:)];

        fn = faceNormal(TR, i);   % 1x3 unit normal, true 3-D orientation
        n2 = fn(cols);
        nrm = norm(n2);
        if nrm < 1e-6
            normals(nSeg,:) = [NaN, NaN];
        else
            normals(nSeg,:) = n2 / nrm;
        end
    end
    % (degenerate cases -- plane passing exactly through a vertex/edge --
    % are silently skipped; they're rare and don't affect the outline)
end

segs    = segs(1:nSeg, :);
normals = normals(1:nSeg, :);

if rotationDeg ~= 0
    segs(:,1:2) = rotate2D(segs(:,1:2), rotationDeg);
    segs(:,3:4) = rotate2D(segs(:,3:4), rotationDeg);
    normals     = rotate2D(normals, rotationDeg);
end

end

function ptsOut = rotate2D(ptsIn, angleDeg)
% Standard CCW rotation about the origin. For points this rotates about
% the origin of the (c1,c2) plane coordinates (not the outline's own
% centroid); for direction vectors (normals) this is exact regardless of
% centering, since there's no translation involved.
theta = deg2rad(angleDeg);
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
ptsOut = (R * ptsIn')';
end
