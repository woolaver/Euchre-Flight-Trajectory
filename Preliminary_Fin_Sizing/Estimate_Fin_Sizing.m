%estimate fin sizing required for 15g 
%turns out at worst possible condition we are nowhere close to meeting 15g,
%CL_req is 5.4427, this is with Virginia which is a decently high lifting
%body and the only parameter we can change is S_ref which is capped by
%sizing requirements so I'm assuming then we only have to be able to do a
%15g turn at any condition, going to use Mach 3 at 3000 m as this is where
%we actually would want to perform this kind of manuever

%capability at Mach 3 caps ~ 13000 m

%at 3000m alpha required is 6.6332, at 0m, alpha required is 5.1256 deg

%NOTE: This code does not currently take into account the loss of lift
%caused by deflecting the fin, will need an iterative process to find this
%but its 1:30 and I'm tired

clear
clc
close all

h = 3000; %m
[~, a, ~, rho] = atmoscoesa(h);

Mach = 3;
V = Mach*a; %m/s

m = 120; %kg
g = 9.8; %m/s^2

W = m*g; %N

L_b = 15*W; %N, 15 G turn

Sref_lift_b = .439; %m^2, using virginia for estimate, pretty similar for all designs

CL_req_b = L_b/(.5*rho*V^2*Sref_lift_b);

CL_vec = readmatrix('Virginia_blended_CL.csv');

alpha_req_b = findAlphaReq(CL_req_b, Mach, CL_vec);

disp("Angle of Attack required for 15G turn: " + rad2deg(alpha_req_b) + " deg")

c = 2; %m, reference length of vehicle

x_cp = .7*c; %m, estimate
x_cg = .6*c; %m, estimate

%note: estimating that x_cp is behind x_cg for now, may need to change this
%if x_cp turns out to be ahead of x_cg

x_f = .9*c; %m, fin cp location, estimate

M_b = -1*(L_b * (x_cp - x_cg)); %N*m, moment about cg due to body lifting force

L_f_req = M_b/(x_f - x_cg); %N, lifting force required from fins to balance

%shock-expansion theory over fin to get lifting force

delta_fin = deg2rad(-20); %deflection of fin

alpha_fin = alpha_req_b + delta_fin;

gamma = 1.4;

%top surface expansion

nu1 = pm(Mach);
nu2 = alpha_fin + nu1;
Mach_upper = ipm(nu2, Mach);

pu_pinf = ((1 + ((gamma-1)/2)*Mach^2)/(1 + ((gamma-1)/2)*Mach_upper^2))^(gamma/(gamma-1));
Cpu = 2/(gamma*Mach^2)*(pu_pinf - 1);

%bottom surface shock

[pl_pinf, ~, ~] = tbm(Mach, alpha_fin);
Cpl = 2/(gamma*Mach^2)*(pl_pinf - 1);

CL_fin = Cpl*cos(alpha_fin) - Cpu*cos(alpha_fin);
CD_fin = Cpl*sin(alpha_fin) + Cpu*sin(alpha_fin);

S_fin = 0; %m^2, set to 0 to get while loop to run
L_fin = .5*rho*V^2*S_fin*CL_fin;

while (L_fin > L_f_req)
    S_fin = S_fin + .005;
    L_fin = .5*rho*V^2*S_fin*CL_fin;
end

S_f_req = S_fin;

disp("Preliminary Fin Size Required at " + rad2deg(delta_fin) + " degree deflection: " + S_f_req + " m^2")