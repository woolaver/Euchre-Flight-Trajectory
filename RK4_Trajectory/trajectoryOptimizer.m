%Script to optimize trajectory
%Parameters: launch angle, alpha climb, level height, alpha angle, dive
%speed, and alpha dive
%Constraints: Mach 3 landing speed, 30000 m height limit, max internal temp below 350K, alphas capped from -20 to 20 degrees
%Goal: maximize range

clc
clear
close all

vehicleProperties = half_cone_properties();
init_launch_angle = vehicleProperties.launch_angle;
init_alpha_climb = vehicleProperties.alpha_climb;
init_level_height = vehicleProperties.level_height;
init_alpha_level = vehicleProperties.alpha_level;
init_dive_speed = vehicleProperties.dive_speed;
init_alpha_dive = vehicleProperties.alpha_dive;

%Parameters for optimizer to update
parameters = [init_launch_angle, init_alpha_climb, init_level_height, init_alpha_level, init_dive_speed, init_alpha_dive];

%Bounds on parameters
upper_bounds = [deg2rad(90), deg2rad(20), 30000, deg2rad(20), 2700, deg2rad(20)];
lower_bounds = [deg2rad(0), deg2rad(-20), 0, deg2rad(-20), 0, deg2rad(-20)];

%constraints
max_height_constraint = 30000;
impact_mach_constraint = 3;
internal_temp_constraint = 350;

%set up surrogate model
options = optimoptions("surrogateopt", 'Display', 'iter', 'MaxFunctionEvaluations', 100, 'UseParallel', false, 'PlotFcn', 'surrogateoptplot', 'InitialPoints', parameters);

func = @(x) obj(x, vehicleProperties, max_height_constraint, impact_mach_constraint, internal_temp_constraint);

[optimal_parameters, fval, exitflag, output, trials] = surrogateopt(func, lower_bounds, upper_bounds, options);

disp("Optimal Parameters:")
disp("Launch Angle: " + rad2deg(optimal_parameters(1)) + " rad")
disp("Alpha Climb: " + rad2deg(optimal_parameters(2)) + " rad")
disp("Level Height: " + optimal_parameters(3) + " m")
disp("Alpha Level: " + rad2deg(optimal_parameters(4)) + " rad")
disp("Dive Speed: " + optimal_parameters(5) + " m/s")
disp("Alpha Dive: " + rad2deg(optimal_parameters(6)) + " deg")

function out = obj(x, vehicleProperties, max_height_constraint, impact_mach_constraint, internal_temp_constraint)
    
    %set parameters to vehicleProperties so they can be passed into the
    %trajectory
    vehicleProperties.launch_angle = x(1);
    vehicleProperties.alpha_climb = x(2);
    vehicleProperties.level_height = x(3);
    vehicleProperties.alpha_level = x(4);
    vehicleProperties.dive_speed = x(5);
    vehicleProperties.alpha_dive = x(6);

    height_init = 1;
    x_init = 0;

    [~, a, ~, ~] = atmoscoesa(height_init);

    Mach_init = 8;
    V_init = a*Mach_init;

    state0 = [V_init; vehicleProperties.launch_angle; height_init; x_init];

    tempOffset = 0;

    try
        [t, state, usefulData] = RK4_3DoF(state0, vehicleProperties, tempOffset);
    catch
        out.Fval = 10000;
        return
    end

    max_height = max(state(3, :));

    impact_mach = usefulData(6, end);

    internal_temp = max(usefulData(8, :));

    max_internal_temp = max(internal_temp);

    range = state(4, end);

    out.Fval = -range;
    out.Ineq = [max_height - max_height_constraint; impact_mach_constraint - impact_mach; max_internal_temp - internal_temp_constraint];
end