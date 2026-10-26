function [CL, CD, CA, CN, CY, Cp, sinTheta, x_cp] = newtonianCLCD3D(stlFile, alpha, minf, Sref_lift, Sref_drag, gamma)
% NEWTONIANCLCD3D Modified Newtonian Theory over a full 3D STL mesh
%
% Inputs:
%   stlFile  - stlFile path
%   alpha    - angle of attack, radians
%  
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


TR = stlread(stlFile);

if(stlFile == "Virginia.stl" || stlFile == "Bunny_Bomb - Revolve2.stl"|| stlFile == 'N:\dotfiles\Desktop\euchre\blend_sharp.STL')
    V = TR.Points; %Use for Bunny_Bomb/Virgina
else
    V = TR.Points./1000; %Use for HARV/Ice Cream Cone
end
F = TR.ConnectivityList;

n = faceNormal(TR);              % Nf x 3, outward unit normals

v1 = V(F(:,1),:); v2 = V(F(:,2),:); v3 = V(F(:,3),:);
areaVec = 0.5 * cross(v2 - v1, v3 - v1, 2);
dA = sqrt(sum(areaVec.^2, 2));   % Nf x 1, true triangle areas

%get center point coordinate for each panel
centroid = (v1 + v2 + v3) / 3;  % Nf x 3 panel center coordinates

%currently set so that nose is pointing in positive z direction because
%thats what trevor did for some reason
Vinf = [0, -sin(alpha), cos(alpha)];

sinTheta = n * Vinf';             % Nf x 1
Cp = zeros(size(sinTheta));


p_pinf = (((gamma+1)^2*minf^2)/((4*gamma*minf^2) - 2*(gamma - 1)))^(gamma/(gamma-1))*((1-gamma+(2*gamma*minf^2))/(gamma+1)); 
CpMax = 2/(gamma*minf^2)*(p_pinf - 1);


windward = sinTheta > 0;
Cp(windward) = CpMax * sinTheta(windward).^2;

dF = -Cp .* n .* dA;              % Nf x 3, force contribution per face
Ftotal = sum(dF, 1);              % 1 x 3

%solvng for center of pressure
%multiply df for cz by x, divide by total length
x_axial = centroid(:,3);
x_cp = sum(dF(:,2).*x_axial)./sum(dF(:,2));
x_cp=(x_cp - min(x_axial))/(max(x_axial) - min(x_axial));


CA = -Ftotal(3) / Sref_drag;
CY = Ftotal(1) / Sref_drag;
CN = Ftotal(2) / Sref_lift;

CL = CN*cos(alpha) - CA*sin(alpha);
CD = CN*sin(alpha) + CA*cos(alpha);

end