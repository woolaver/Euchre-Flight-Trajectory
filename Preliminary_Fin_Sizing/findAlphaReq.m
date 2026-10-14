function alpha_req = findAlphaReq(CL_req, Mach, CL_vec)
%finds the required alpha to meet a certain CL (note: function averages the CL
%immediatley above and below the CL requried, if this value is too far off
%then the fineness of CL_vec_file will need to be increased
%Inputs:
%   Cl_req: required CL to meet condition
%   Mach: current Mach number
%   CL_vec_file: file to use for CL values
%Outputs:
%   alpha_req: required angle of attack to meet condition (rad)
    
    %find proper row to use in CL_vec_file (note: takes row above the
    %required value as odds are the exact Mach number will not line up

    Mach_vec = linspace(2, 8, 200);
    alpha_vec = deg2rad(linspace(-30, 30, 200));

    Mach_row = find(Mach_vec >= Mach, 1, "first");

    CL_vec_at_mach = CL_vec(Mach_row, :);

    alpha_req_index = find(CL_vec_at_mach > CL_req, 1, "first");

    alpha_req = (alpha_vec(alpha_req_index) + alpha_vec(alpha_req_index - 1))/2;
end