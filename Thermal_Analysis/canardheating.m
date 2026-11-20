function [q_fin] = canardheating(vehicleProperties,h,minf,alpha_f,t_f)
%very rough estimate, using rectangular canard and eckharts 
%constants 
pr = .75;
r = sqrt(pr);
cp=1004;
y=1.4;
R=287;
[rinf,pinf,tinf] = stdatmfull(h);


%geometry
fin_length = vehicleProperties.fin_length;
fin_height = vehicleProperties.fin_thickness/2;
Sref_fin = vehicleProperties.Sref_fin;
fin_width = vehicleProperties.fin_width;
panel_length = sqrt((fin_length/2)^2 + fin_height^2);
theta = atan(fin_height/(fin_length/2)); %half angle of diamond airfoil


%%
%1 upper surface leading panel, 2 is lower surface leading panel, 3 is
%upper surface trailing panel, 4 is lower surface trailing panel

theta1 = theta - alpha_f;
theta2 = theta + alpha_f;
theta3 = 2*theta;
theta4 = 2*theta;

%compression
if(theta1 > 0)
    [p1_pinf, t1_tinf, r1_rinf, ~, M1] = tbm(minf, theta1);
    p1 = p1_pinf*pinf;
    t1 = t1_tinf*tinf;
    r1 = r1_rinf*rinf;
    u1 = M1*sqrt(y*R*t1);
%expansion
else
    nu_s = pm(minf);
    nu_1 = nu_s + abs(theta1);
    M1 = ipm(nu_1, minf);

    ps_pt = (1 + (y-1)/2 * minf^2)^(-y/(y-1));
    p1_pt = (1 + (y-1)/2 * M1^2)^(-y/(y-1));
    p1_ps = p1_pt/ps_pt;
    p1 = p1_ps*pinf;

    ts_tt = (1 + (y-1)/2 * minf^2);
    t1_tt = (1 + (y-1)/2 * M1^2);
    t1_ts = t1_tt/ts_tt;
    t1 = t1_ts*tinf;
    r1 = p1/R/t1;
    
    u1 = M1*sqrt(y*R*t1);
end

%compression
if(theta2 > 0)
    [p2_pinf, t2_tinf, r2_rinf, ~, M2] = tbm(minf, theta2);
    p2 = p2_pinf*pinf;
    t2 = t2_tinf*tinf;
    r2 = r2_rinf*rinf;
    u2 = M2*sqrt(y*R*t2);
%expansion
else
    nu_s = pm(minf);
    nu_2 = nu_s + abs(theta2);
    M2 = ipm(nu_2, minf);

    ps_pt = (1 + (y-1)/2 * minf^2)^(-y/(y-1));
    p2_pt = (1 + (y-1)/2 * M1^2)^(-y/(y-1));
    p2_ps = p2_pt/ps_pt;
    p2 = p2_ps*pinf;

    ts_tt = (1 + (y-1)/2 * minf^2);
    t2_tt = (1 + (y-1)/2 * M2^2);
    t2_ts = t2_tt/ts_tt;
    t2 = t2_ts*tinf;
    r2 = p2/R/t2;
    u2 = M2*sqrt(y*R*t2);
end

%expansion over surface 3 and 4
nu_1 = pm(M1);
nu_3 = nu_1 + abs(theta3);
M3 = ipm(nu_3, M1);

p1_pt = (1 + (y-1)/2 * M1^2)^(-y/(y-1));
p3_pt = (1 + (y-1)/2 * M3^2)^(-y/(y-1));
p3_p1 = p3_pt/p1_pt;
p3 = p3_p1*p1;

t1_tt = (1 + (y-1)/2 * M1^2);
t3_tt = (1 + (y-1)/2 * M3^2);
t3_t1 = t3_tt/t1_tt;
t3 = t3_t1*t1;
r3 = p3/R/t3;
u3 = M3*sqrt(y*R*t3);

%-------------------------------------------------------
nu_2 = pm(M2);
nu_4 = nu_2 + abs(theta4);
M4 = ipm(nu_4, M2);

p2_pt = (1 + (y-1)/2 * M2^2)^(-y/(y-1));
p4_pt = (1 + (y-1)/2 * M4^2)^(-y/(y-1));
p4_p2 = p4_pt/p2_pt;
p4 = p4_p2*p2;

t2_tt = (1 + (y-1)/2 * M2^2);
t4_tt = (1 + (y-1)/2 * M4^2);
t4_t2 = t4_tt/t2_tt;
t4 = t4_t2*t2;
r4 = p4/R/t4;
u4 = M4*sqrt(y*R*t4);
%%

%eckharts
haw1 = r*.5*u1^2+cp*t1;
ts1 = t1*(1+.032*M1^2+.58*(t_f/t1-1));
rs1 = p1/R/ts1;
mus1 = calcSutherlands(ts1);
qw1 = .332*(pr^(-2/3))/sqrt(rs1*u1/mus1)*rs1*u1*(haw1-cp*t_f)*2*sqrt(panel_length); %integrated over panel length

haw2 = r*.5*u2^2 + cp*t2;
ts2 = t2*(1 + .032*M2^2 + .58*(t_f/t2 - 1));
rs2 = p2/(R*ts2);
mus2 = calcSutherlands(ts2);
qw2 = .332*(pr^(-2/3))/sqrt(rs2*u2/mus2)*rs2*u2*(haw2 - cp*t_f)*2*sqrt(panel_length);

haw3 = r*.5*u3^2 + cp*t3;
ts3 = t3*(1 + .032*M3^2 + .58*(t_f/t3 - 1));
rs3 = p3/(R*ts3);
mus3 = calcSutherlands(ts3);
qw3 = .332*(pr^(-2/3))/sqrt(rs3*u3/mus3)*rs3*u3*(haw3 - cp*t_f)*2*sqrt(panel_length);

haw4 = r*.5*u4^2 + cp*t4;
ts4 = t4*(1 + .032*M4^2 + .58*(t_f/t4 - 1));
rs4 = p4/(R*ts4);
mus4 = calcSutherlands(ts4);
qw4 = .332*(pr^(-2/3))/sqrt(rs4*u4/mus4)*rs4*u4*(haw4 - cp*t_f)*2*sqrt(panel_length);

q_fin = (qw1 + qw2 + qw3 + qw4)*fin_width;
end