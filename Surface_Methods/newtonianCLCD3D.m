function [CL, CD, CA, CN, CY, Cp, sinTheta] = newtonianCLCD3D(stlFile, alpha, CpMax, Sref)
% NEWTONIANCLCD3D Modified Newtonian Theory over a full 3D STL mesh (no
% planar slicing -- every triangle in the mesh contributes).
%
% [CL, CD, CA, CN, CY, Cp] = newtonianCLCD3D(stlFile, alphaDeg, CpMax, Sref)
%
% Inputs:
%   stlFile  - path to .stl file
%   alpha    - angle of attack, radians. ASSUMES body axes with x = axial
%              (nose-to-tail) and z = "up" -- freestream direction is
%              taken as [cos(alpha), 0, sin(alpha)]. If your STL's axes
%              don't match this (e.g. z is spanwise instead), rotate the
%              mesh first or edit the Vinf line below.
%   CpMax    - stagnation Cp behind a normal shock at the flight Mach
%              number (see rayleighPitotCpMax.m)
%   Sref     - reference AREA for nondimensionalizing the force integral
%              (e.g. planform area or a chosen cross-sectional area).
%              This replaces the reference LENGTH used in the old 2D
%              per-unit-span version -- a true 3D integral needs an
%              area, not a length. Keep it consistent with whatever
%              convention the rest of your framework uses.
%
% Outputs:
%   CL, CD     - lift and drag coefficients, wind axes
%   CA, CN, CY - axial, normal, and side force coefficients, body axes
%                (CY is a useful symmetry check: for a body with no
%                sideslip and a left-right symmetric mesh, CY should
%                come out near zero -- a nonzero CY flags either real
%                asymmetry in the geometry or a normal-orientation bug)
%   Cp         - Nfaces x 1 vector, per-triangle pressure coefficient
%                (same face order as the STL's ConnectivityList -- pass
%                straight into plotVehicleCp3D for visualization)
%
% Method: for each triangle, sin(theta) = dot(freestream_dir, outward
% normal). Windward faces (sin(theta) > 0) get Cp = CpMax*sin(theta)^2;
% shadowed faces get Cp = 0 (basic Newtonian shadowing -- same caveat as
% before: a leeside correction would improve accuracy here). Force on
% each face acts along -outward_normal, weighted by the face's true
% area; sum over all faces and rotate by alpha into wind axes.

TR = stlread(stlFile);
V = TR.Points;
F = TR.ConnectivityList;

n = faceNormal(TR);              % Nf x 3, outward unit normals

v1 = V(F(:,1),:); v2 = V(F(:,2),:); v3 = V(F(:,3),:);
areaVec = 0.5 * cross(v2 - v1, v3 - v1, 2);
dA = sqrt(sum(areaVec.^2, 2));   % Nf x 1, true triangle areas

%currently set so that nose is pointing in positive z direction because
%thats what trevor did for some reason
Vinf = [0, -sin(alpha), cos(alpha)];

sinTheta = n * Vinf';             % Nf x 1
Cp = zeros(size(sinTheta));
windward = sinTheta > 0;
Cp(windward) = CpMax * sinTheta(windward).^2;

dF = -Cp .* n .* dA;              % Nf x 3, force contribution per face
Ftotal = sum(dF, 1);              % 1 x 3

CA = -Ftotal(3) / Sref;
CY = Ftotal(1) / Sref;
CN = Ftotal(2) / Sref;

CL = CN*cos(alpha) - CA*sin(alpha);
CD = CN*sin(alpha) + CA*cos(alpha);

end
