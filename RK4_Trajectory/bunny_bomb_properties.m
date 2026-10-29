function vehicleProperties = bunny_bomb_properties()
%function that contains all vehicle properties of Bunny Bomb vehicle
%call this function when running trajectory for Bunny_Bomb round vehicle
%Outputs:
%   vehicleProperties: struct containing all necessary vehicle properties
%   for trajectory analysis

    vehicleProperties = struct();

    %aero coefficient csv files
    vehicleProperties.CL_vec_file = "csv_files/Bunny_Bomb_blended_CL.csv";
    vehicleProperties.CD_vec_file = "csv_files/Bunny_Bomb_blended_CD.csv";
    
    %trajectory save file
    vehicleProperties.Trajectory_Save_File = "Bunny_Bomb_Trajectory.csv";
    
    %dimensions
    vehicleProperties.Sref = .467; %m^2
    vehicleProperties.length = 1.950; %m
    vehicleProperties.height = .160; %m, height is from center of vehicle to top point (radius for perfectly circular cross section)
    vehicleProperties.A = 1.094; %m^2, total surface area
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
    vehicleProperties.alpha_climb = deg2rad(3);
    vehicleProperties.dt = .01; %sec
    vehicleProperties.launch_angle = deg2rad(18.5);
end