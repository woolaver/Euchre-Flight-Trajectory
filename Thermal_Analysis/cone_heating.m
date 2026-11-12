function cone_heating(stlFile,trajectory,kval,thick)
%cone heating from stl file of nosecone
%Inputs:
%   stl file 
%   trajectory 'Half_Cone_Trajectory'
%   r = .1 %nose radius (m) ask adam
%   kvalmatrix of
%Outputs:
%   Ffx: X-direction friction force (N)
%   Ffy: Y-direction friction force (N)
%   Q_tot: Total heat flux (W)

%constants
nosetip = [ 0 0 0]; %this is the cordinates of the tip
R = insulation_calc


%to get angles, straight claude
TR = stlread(stlFile);
V = TR.Points; %Use for Bunny_Bomb/Virgina
F = TR.ConnectivityList;
n = faceNormal(TR);              % Nf x 3, outward unit normals
v1 = V(F(:,1),:); 
v2 = V(F(:,2),:); 
v3 = V(F(:,3),:);
areaVec = 0.5 * cross(v2 - v1, v3 - v1, 2);
dA = sqrt(sum(areaVec.^2, 2));   % Nf x 1, true triangle areas
%get center point coordinate for each panel
centroid = (v1 + v2 + v3) / 3;  % Nf x 3 panel center coordinates

%create a matrix for the surface temp to use at each prior step
t_s_old = 298*ones(size(centroid))
%starting for loop to go though every timestep
for i =1:size(trajectory,2)
%getting data from trajectory
alpha = trajectory(i,6);
alt = trajectory(i,6);
minf = trajectory(i,6);

vinf = [0, -sin(alpha), -cos(alpha)];

theta = asin(n * vinf');        %nfx1 this is dot product between normals and the velocity vector to find theta 
q = size(theta); %this is a matrix of all the q's
distance = sqrt((centroid(1)-nosetip(1)).^2+(centroid(2)-nosetip(2)).^2+(centroid(2)-nosetip(2)).^2);

for j = 1:size(theta,2)
    if theta(j)>deg2rad(85)  %basically if the wall is flat, use fayriddel to get
    
    else


    end

end

end