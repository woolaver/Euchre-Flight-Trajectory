function vehicleProperties = virginia_properties()
%function that contains all vehicle properties of virginia vehicle
%call this function when running trajectory for virginia round vehicle
%Outputs:
%   vehicleProperties: struct containing all necessary vehicle properties
%   for trajectory analysis

    vehicleProperties = struct();

    %aero coefficient csv files
    vehicleProperties.CL_vec_file = "csv_files/Virginia_blended_CL.csv";
    vehicleProperties.CD_vec_file = "csv_files/Virginia_blended_CD.csv";
    
    %trajectory save file
    vehicleProperties.Trajectory_Save_File = "Virginia_Trajectory.csv";
    
    %dimensions
    vehicleProperties.Sref = .439; %m^2
    vehicleProperties.length = 1.953; %m
    vehicleProperties.height = .176; %m, height is from center of vehicle to top point (radius for perfectly circular cross section)
    vehicleProperties.A = 1.478; %m^2, total surface area
    vehicleProperties.wall_thickness = .0254; %m, assuming 1 inch thickness for testing purpose, need structure for real value
    
    %mass properties
    vehicleProperties.m_skin = 60; %kg, initial estimate, will need structure for better estimate
    vehicleProperties.m_structure = 50; %kg, intitial estimate, will need structure for better estimate
    vehicleProperties.m = 120; %kg, total mass
    
    %thermal properties
    vehicleProperties.Tw = 298; %K, initial wall temperature, room temperature
    vehicleProperties.T_int = 298; %K, initial internal temperature, room temperature
    vehicleProperties.cp_w = 523; %J/kg*K, wall specific heat capacity for TI6ALV4
    vehicleProperties.cp_int = 870; %J/kg*K, internal specific heat capacity for 7075-T6 Al
    vehicleProperties.emissivity = .4; % wall emissivity, estimate for TI6ALV4
    vehicleProperties.k_w = 7; %W/m*K, estimate for TI6ALV4, thermal conductivity of wall
    vehicleProperties.k_int = 130; %W/m*K, estimate for 7075-T6 Al, thermal conductivity of internal structure

    %trajectory properties
    vehicleProperties.alpha_climb = deg2rad(0);
    vehicleProperties.dt = .1; %sec
    vehicleProperties.launch_angle = deg2rad(65);
    vehicleProperties.alpha_glide = deg2rad(0);
    vehicleProperties.alpha_dive = deg2rad(-5);
    vehicleProperties.alpha_level = deg2rad(-13);
    vehicleProperties.dive_speed = 1410;

    %nose properties (nose properties are the same as blend flat)
    vehicleProperties.R_nose = .01; %m
    vehicleProperties.A_nose = 1.578*10^(-4); %m^2
    vehicleProperties.k_nose = 30; %W/m*K
    vehicleProperties.cp_nose = 1200; %J/kg*K, estimate for high temperatures
    vehicleProperties.m_nose = .0059; %kg, estiamte of 2800 kg/m^3, volume calculated as hemisphere with radius R_nose
    vehicleProperties.emissivity_nose = .8; %SiC/SiC CMC estimate
    vehicleProperties.D_nose = .02; %m^2, circular esimtate from R_nose
    vehicleProperties.T_nose = 298; %K, initial nose temperature, room temperature
end