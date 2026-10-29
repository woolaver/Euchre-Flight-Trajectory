function [CL, CD, Cx, Cz, Cy, Cp, theta, centroid, x_cp] = smalldist_supersonic(stlFile, alpha, minf, sref, y)
%The aim of this script is to make a low fidelity solver for supersonic 
%will miss blunt edges correct later
TR = stlread(stlFile);

if(stlFile == "csv_files/Virginia.STL" || stlFile == "csv_files/Bunny_Bomb - Revolve2.STL")
    V = TR.Points; %Use for Bunny_Bomb/Virgina
else
    V = TR.Points./1000; %Use for HARV/Ice Cream Cone
end

F = TR.ConnectivityList;

n = faceNormal(TR);              % Nf x 3, outward unit normals

v1 = V(F(:,1),:); 
v2 = V(F(:,2),:); 
v3 = V(F(:,3),:);
areaVec = 0.5 * cross(v2 - v1, v3 - v1, 2);
dA = sqrt(sum(areaVec.^2, 2));   % Nf x 1, true triangle areas

%get center point coordinate for each panel
centroid = (v1 + v2 + v3) / 3;  % Nf x 3 panel center coordinates

% Vinf = [cos(alpha), 0, sin(alpha)]; % 1x3 da vector for velocity

if(stlFile == "csv_files/blend_flat.STL" || stlFile == "csv_files/blend_round.STL" || stlFile == "csv_files/blend_sharp.STL")
    vinf = [0, -sin(alpha), -cos(alpha)];
else
    vinf = [0, -sin(alpha), cos(alpha)];
end

theta = asin(n * vinf');        %nfx1 this is dot product between normals and the velocity vector to find theta    
Cp = 2.*theta./(sqrt(minf.^2-1));

%still feel a little sus about this
%need if statement to prevent cp from getting too insane

p02_pinf=(((y+1)^2*minf^2)/(4*y*minf^2- 2*(y-1)))^(y/(y-1))*(1-y+2*y*minf^2)/(y+1);

Cp_max = (2 / (y * minf^2)) * (p02_pinf - 1);
Cp_min = -2 / (y * minf^2); %if p02/pinf was 0 basically
%disp(p02_pinf);
%disp(Cp_max);
%disp(Cp_min);
Cp = max(min(Cp, Cp_max), Cp_min);   % going through and min maxing the cp points


dF = -Cp .* n .* dA;              % Nf x 3, force contribution per face
Ftotal = sum(dF, 1);              % 1 x 3


%solvng for center of pressure
%multiply df for cz by x, divide by total length
x_axial = centroid(:,3);
x_cp = sum(dF(:,2).*x_axial)./sum(dF(:,2));
x_cp=(x_cp - min(x_axial))/(max(x_axial) - min(x_axial));

if(stlFile == "csv_files/blend_flat.STL" || stlFile == "csv_files/blend_round.STL" || stlFile == "csv_files/blend_sharp.STL")
    Cx = Ftotal(3) / sref;  %forces in x direction, nose to tail should be 1
    Cy = Ftotal(1) / sref;  %forces in y drection, out wing,  shouldnt need this should be 2
    Cz = Ftotal(2) / sref;  %forces in z, up should be 3 but bunny bomb is cooked
else
    Cx = -Ftotal(3) / sref;  %forces in x direction, nose to tail should be 1
    Cy = Ftotal(1) / sref;   %forces in y drection, out wing,  shouldnt need this should be 2
    Cz = Ftotal(2) / sref;   %forces in z, up should be 3 but bunny bomb is cooked
end

CL = Cz*cos(alpha) - Cx*sin(alpha);
CD = Cz*sin(alpha) + Cx*cos(alpha);

end