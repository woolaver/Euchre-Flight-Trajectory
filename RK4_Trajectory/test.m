clear
clc
close all

vehicleProperties = blend_round_properties();

height_init = 1;

[~, a, ~, ~] = atmoscoesa(height_init);

Mach_init = 8;
V_init = a*Mach_init;

launch_angle = vehicleProperties.launch_angle; 
x_init = 0;

state0 = [V_init; launch_angle; height_init; x_init];

dt = vehicleProperties.dt;
    
climb_end = 0;
cruise_end = 1;
state(:, 1) = state0;
t(1) = 0;

m = vehicleProperties.m; 
g = 9.8; %m/s^2

CL_vec = readmatrix(vehicleProperties.CL_vec_file);
CD_vec = readmatrix(vehicleProperties.CD_vec_file);

%populate first column of usefulData
[~, a, ~, rho] = atmoscoesa(state0(3));
    
Mach = state(1, 1)/a;

%assume CL and CD are 0 for the first input
CL = 0;
CD = 0;

Tw = vehicleProperties.Tw; %K, room temperature starting temperature
T_int = vehicleProperties.T_int; %K, room temperature starting temperature
T_nose = vehicleProperties.T_nose; %K, room temperature starting temperature

alpha_climb = vehicleProperties.alpha_climb;
Sref = vehicleProperties.Sref;
length = vehicleProperties.length;
height = vehicleProperties.height;

usefulData(:, 1) = [0; CL; CD; rho; Mach; Tw; T_int; T_nose];

%alpha_vec and Mach_vec must be the same size and steps used to
%create CL_vec_file and CD_vec_file for interpolation to work
alpha_vec = deg2rad(linspace(-30, 30, 200));
Mach_vec = linspace(2, 8, 200);

for i = 1:2

    %disp("Height: " + state(3, i) + " m")
    %disp("Velocity: " + state(1, i) + "m/s")
    
    [~, a, ~, rho] = atmoscoesa(state(3, i));
    Mach = state(1, i)/a;
    
    %for now hard set alpha constant during climb
    alpha = alpha_climb;
    
    CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha, Mach, 'linear');
    CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha, Mach, 'linear');
    
    if(isnan(CL))
        CL = 0;
        %disp("CL returned Nan, default value of 0 is used")
    end
    
    if(isnan(CD))
        CD = .1;
        %disp("CD returned NaN, default value of .1 is used")
    end
    
    D = 1/2*rho*state(1, i)^2*Sref*CD;
    L = 1/2*rho*state(1, i)^2*Sref*CL;
    
    [Fx, Fy, Q_conv_rate] = refTempMethod(length, height, state(3, i), Mach, alpha, Tw);
    Q_conv_nose_rate = fayRiddellMethod(vehicleProperties, Mach, state(3, i), T_nose);
    
    D = D + Fx;
    L = L + Fy;
    
    [Tw, T_int, T_nose] = updateTemperatures(vehicleProperties, state(3, i), Tw, T_int, T_nose, Q_conv_rate, Q_conv_nose_rate);
    
    %{
    state_dot = [-D/m - g*sin(state(2));
             L/(m*state(1)) - (g/state(1))*cos(state(2));
             state(1)*sin(state(2));
             state(1)*cos(state(2))];
    %}
    
    %RK-4 iteration
    %k1 = dt*f(t, yn)
    %k2 = dt*f(t + dt/2, yn + k1/2)
    %k3 = dt*f(t + dt/2, yn + k2/2)
    %k4 = dt*f(t + dt, yn + k3)
    %yn+1 = yn + dt/6(k1 + 2k2 + 2k3 + k4)
    
    k1 = [-D/m - g*sin(state(2, i));
             L/(m*state(1, i)) - (g/state(1, i))*cos(state(2, i));
             state(1, i)*sin(state(2, i));
             state(1, i)*cos(state(2, i))];
    
    k2 = [-D/m - g*sin(state(2, i) + k1(2)/2);
             L/(m*(state(1, i) + k1(1)/2)) - (g/(state(1, i) + k1(1)/2))*cos((state(2, i) + k1(2)/2));
             (state(1, i) + k1(1)/2)*sin((state(2, i) + k1(2)/2));
             (state(1, i) + k1(1)/2)*cos((state(2, i) + k1(2)/2))];
    
    k3 = [-D/m - g*sin(state(2, i) + k2(2)/2);
             L/(m*(state(1, i) + k2(1)/2)) - (g/(state(1, i) + k2(1)/2))*cos((state(2, i) + k2(2)/2));
             (state(1, i) + k2(1)/2)*sin((state(2, i) + k2(2)/2));
             (state(1, i) + k2(1)/2)*cos((state(2, i) + k2(2)/2))];
    
    k4 = [-D/m - g*sin(state(2, i) + k3(2));
             L/(m*(state(1, i) + k3(1))) - (g/(state(1, i) + k3(1)))*cos((state(2, i) + k3(2)));
             (state(1, i) + k3(1))*sin((state(2, i) + k3(2)));
             (state(1, i) + k3(1))*cos((state(2, i) + k3(2)))];
    
    state(:, i+1) = state(:, i) + (dt./6).*(k1 + 2.*k2 + 2.*k3 + k4)
    
    usefulData(:, i+1) = [alpha, CL, CD, rho, Mach, Tw, T_int, T_nose];
    
    t(i+1) = t(i)+dt;
end