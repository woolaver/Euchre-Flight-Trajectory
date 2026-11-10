function [Lift, Drag, alpha_req, fin_alpha] = findAlphaReq(vehicleProperties, L_req, Mach, alt)
%finds the required alpha to meet a certain lift value (note: function averages the CL
%immediatley above and below the CL requried, if this value is too far off
%then the fineness of CL_vec_file will need to be increased
%Inputs:
%   vehicleProperties: struct containing all necessary vehicle properties
%   L_req: required lift to meet condition
%   Mach: current Mach number
%   alt: current altitude
%Outputs:
%   Lift: final lift value (N)
%   Drag: final drag value (N)
%   alpha_req: required body angle of attack to meet condition (rad)
%   fin_alpha: angle of attack of fins to trim vehicle (rad)
    
    %find proper row to use in CL_vec_file (note: takes row above the
    %required value as odds are the exact Mach number will not line up

    Mach_vec = linspace(1.5, 9, 200); 
    alpha_vec = deg2rad(linspace(-30, 30, 200));

    Mach_row = find(Mach_vec >= Mach, 1, "first");

    CL_vec_at_mach = CL_vec(Mach_row, :);

    alpha_req_index = find(CL_vec_at_mach > CL_req, 1, "first");

    alpha_req = (alpha_vec(alpha_req_index) + alpha_vec(alpha_req_index - 1))/2;
end