function [t_n, t_i, qan,qni] = cone_heating(stlFile, trajectory,kval,thick)
%cone heating from stl file of nosecone assuming that nose cone can be
%treated as a lumped system and doesnt use any dissasoiation for the
%tangent wedge part
%Inputs:
%   stl file 
%   trajectory 'Half_Cone_Trajectory'
%   kval matrix for each layer
%   thickness matrix of layers
%Outputs:
%   t_n nose cone temp over time
%   t_i structre temp over time

%constants

trajectory = readtable(trajectory);
alphas = trajectory.AngleOfAttack;%angle of attack in radians
machs = trajectory.MachNumber;
height = trajectory.Height;




nosetip = [0 0 0]; %this is the cordinates of the tip
surfacearea = .018941; %exposed surface area of nose cone for radiation
cp_tung = 141; %j/kg*k
m_n  = 3.54; %mass of nose (kg)
cp_al = 880; %j/kg*k for 7075
m_i  = 8; %mass of internal, probly just going nto do the mass of the nearest rib i guess
R_ins = insulation_calc(kval,thick);
r=.01; %nose cone radius
backarea = .002717;



R= 287;
sb = 5.67e-8; %stefanboltzman constant
dt = .1;
y=1.4;
PR = .75;
recovery  =sqrt(PR);
cp_air = 1004;
emis = .23; %set as .23 from some table

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

%find the distance from the nose to the centroid of each panel
distance = max(vecnorm(centroid - nosetip, 2, 2), 0.01);   % Nf x 1, floored at the tip
%create a variables for the temp to use at each prior step and to store fo over time
t_n = zeros(size(trajectory,1),1); t_n(1)=298;%nose cone temp
t_i = zeros(size(trajectory,1),1); t_i(1)=298;%internal temp
qan = zeros(size(trajectory,1),1);
qni = zeros(size(trajectory,1),1);



%starting for loop to go though every timestep
for i =2:size(trajectory,1)
%getting data from trajectory
alpha = alphas(i,1);%angle of attack in radians
alt = height(i,1);
minf = machs(i,1);
[~,~,Tamb] = stdatmfull(alt);


vinf = [0, -sin(alpha), -cos(alpha)];

theta = asin(n * vinf');        %nfx1 this is dot product between normals and the velocity vector to find theta 
q = zeros(size(theta)); %this is a matrix of all the q's
qfr = fayRiddellMethod(r, minf, alt, t_n(i-1));
for j = 1:size(theta,1)

    if theta(j)>deg2rad(85)  %basically if the wall is flat, use fayriddel to get
        q(j) = qfr;
    elseif theta(j)<0
        q(j)=0;
    else
        [a,ok] = tangent_wedge(minf,alt,theta(j));
        if~ok
            q(j) = qfr*max(sin(theta(j)),0.1)^1.5;
        else
            tw_te = t_n(i-1)/a.T;
            tstar = a.T*(1 + .032*a.M^2 + .58*(tw_te - 1));
            mustar = calcSutherlands(tstar);
            rstar = a.p/R/tstar;
            rex_star = rstar*a.ue*distance(j)/mustar;
            ch = .332./sqrt(rex_star).*PR^(-2/3);
    
            haw = recovery*.5*a.ue^2+1004*a.T;
            q(j) = ch*rstar*a.ue*(haw-cp_air*t_n(i-1));
        end
    end

end

q_an = sum(q.*dA); %heat flux air to nose
q_ni = (t_n(i-1)-t_i(i-1))/R_ins*backarea;  %heat flux nose to internal
q_nrad = emis*surfacearea*sb*(t_n(i-1)^4-Tamb^4); %radation heatflux from nose to air

qan(i) = q_an;
qni(i) = q_ni;

%update temps
netq = q_an-q_ni-q_nrad;
t_n(i) = t_n(i-1)+dt*netq/(m_n*cp_tung);
t_i(i) = t_i(i-1)+dt*q_ni/(m_i*cp_al);
fprintf(' Internal Temp = %4f k,  Nose Temp = %4f k\n',t_i(i,1),t_n(i,1))
fprintf(' nose q = %4f w,\n\n',q_an)


end
fprintf('Final Internal Temp = %4f k, Final Nose Temp = %4f k\n',t_i(end,end),t_n(end,end))
fprintf('Max Internal Temp = %4f k, Max Nose Temp = %4f k\n',max(t_i),max(t_n))

end