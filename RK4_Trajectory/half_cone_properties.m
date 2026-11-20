function vehicleProperties = half_cone_properties()
%function that contains all vehicle properties of blend sharp vehicle
%call this function when running trajectory for blend sharp vehicle
%Outputs:
%   vehicleProperties: struct containing all necessary vehicle properties
%   for trajectory analysis
    
    vehicleProperties = struct();

    %aero coefficient csv files
    vehicleProperties.CL_vec_file = "csv_files/Half_Cone_blended_CL.csv";
    vehicleProperties.CD_vec_file = "csv_files/Half_Cone_blended_CD.csv";
    
    %trajectory save file
    vehicleProperties.Trajectory_Save_File = "Half_Cone_Trajectory.csv";
    
    %dimensions
    vehicleProperties.Sref = .52; %m^2
    vehicleProperties.length = 1.919; %m
    vehicleProperties.height = .189; %m, height is from center of vehicle to top point (radius for perfectly circular cross section)
    vehicleProperties.A = 1.252; %m^2, total surface area
    vehicleProperties.wall_thickness = .01; %m, assuming 1 inch thickness for testing purpose, need structure for real value
    
    %mass properties
    vehicleProperties.m_skin = 75; %kg, initial estimate, will need structure for better estimate
    vehicleProperties.m_structure = 35; %kg, intitial estimate, will need structure for better estimate
    vehicleProperties.m = 120; %kg, total mass
    
    %thermal properties
    vehicleProperties.Tw = 298; %K, initial wall temperature, room temperature
    vehicleProperties.T_int = 298; %K, initial internal temperature, room temperature
    vehicleProperties.cp_w = 523; %J/kg*K, wall specific heat capacity for TI6ALV4
    vehicleProperties.cp_int = 870; %J/kg*K, internal specific heat capacity for 7075-T6 Al
    vehicleProperties.emissivity = .4; % wall emissivity, estimate for TI6ALV4
    vehicleProperties.k_w = 7; %W/m*K, estimate for TI6ALV4, thermal conductivity of wall
    vehicleProperties.k_int = 130; %W/m*K, estimate for 7075-T6 Al, thermal conductivity of internal structure
    vehicleProperties.k_insulation = .04; %J/kg*K, constant for now but should turn into a linear interpolation
    vehicleProperties.insulation_thickness = .01; %m

    %trajectory properties
    vehicleProperties.alpha_climb = deg2rad(3);
    vehicleProperties.dt = .1; %sec
    vehicleProperties.launch_angle = deg2rad(50);
    vehicleProperties.alpha_glide = deg2rad(20);
    vehicleProperties.alpha_dive = deg2rad(-3);
    vehicleProperties.alpha_level = deg2rad(-13.3);
    vehicleProperties.dive_speed = 1225;
    vehicleProperties.level_height = 15000;

    %nose properties
    vehicleProperties.R_nose = .01; %m
    vehicleProperties.A_nose = .002717; %m^2
    vehicleProperties.k_nose = 172; %W/m*K
    vehicleProperties.cp_nose = 141; %J/kg*K, estimate for high temperatures
    vehicleProperties.m_nose = 3.54; %kg
    vehicleProperties.emissivity_nose = .23; %Tungsten estimate
    vehicleProperties.D_nose = .02; %m^2, circular esimtate from R_nose
    vehicleProperties.T_nose = 298; %K, initial nose temperature, room temperature

    %fin properties
    vehicleProperties.fin_x_location = .4820; %1/4 body length from the nose
    vehicleProperties.fin_y_location = .094; %distance from vehicle edge to the vehicle centerline
    vehicleProperties.Sref_fin = .042; %m^2
    vehicleProperties.fin_thickness = .02; %m, estimated fin thickness for diamond airfoil, results in half angle of 5.7124
    vehicleProperties.fin_length = .2; %m, assuming rectangular planform
    vehicleProperties.fin_width = .1; %m, assuming rectangular planform
    vehicleProperties.cg_x_location = vehicleProperties.length*.5; %estimating cg at half length
    vehicleProperties.cp_x_location = vehicleProperties.length*.55; %estimating cp at .55 length, later on this should change with angle of attack but fine for testing rn
    vehicleProperties.fin_volume = .0001565;
end