function [Lift, Drag, alpha_req, alpha_fin] = findAlphaReq(vehicleProperties, L_req, Mach, alt)
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
%   alpha_fin: angle of attack of fins to trim vehicle (rad)
    
    %find proper row to use in CL_vec_file (note: takes row above the
    %required value as odds are the exact Mach number will not line up\

    [~, a, ~, rho] = atmoscoesa(alt);
    V = Mach*a;
    Sref = vehicleProperties.Sref;

    CL_vec = readmatrix(vehicleProperties.CL_vec_file);
    CD_vec = readmatrix(vehicleProperties.CD_vec_file);

    Mach_vec = linspace(1.5, 9, 200); 
    alpha_vec = deg2rad(linspace(-30, 30, 200));

    Mach_row = find(Mach_vec >= Mach, 1, "first");

    CL_vec_at_mach = CL_vec(Mach_row, :);

    CL_req = L_req/((1/2)*rho*V^2*Sref);

    alpha_req_index = find(CL_vec_at_mach > CL_req, 1, "first");

    alpha_guess = (alpha_vec(alpha_req_index) + alpha_vec(alpha_req_index - 1))/2;
    
    %loop through angles of attack near alpha guess to find one where total
    %lift (trim lift + body lift) equals L_req

    diff = 1000;
    tol = 100;

    delta_alpha_guess = .01;
    prev_jump = "none";

    while (diff > tol)       
        CL_body = interp2(alpha_vec, Mach_vec, CL_vec, alpha_guess, Mach, 'linear');
        CD_body = interp2(alpha_vec, Mach_vec, CD_vec, alpha_guess, Mach, 'linear');

        L_body = (1/2)*rho*V^2*Sref*CL_body;
        D_body = (1/2)*rho*V^2*Sref*CD_body;

        [D_fin, L_fin, alpha_fin] = calcTrimDrag(vehicleProperties, alpha_guess, L_body, D_body, alt, Mach);
        
        L_calc = L_body + L_fin;

        diff = abs(L_req - L_calc);

        if(L_req > L_calc)
            alpha_guess = alpha_guess + delta_alpha_guess;
            if(prev_jump == "decrease")
                delta_alpha_guess = delta_alpha_guess/10;
                prev_jump = "none";
            else
                prev_jump = "increase";
            end
            
        else
            %if switching between increasing and decreasing then decrease
            %delta_alpha_fin
            alpha_guess = alpha_guess - delta_alpha_guess;

            if(prev_jump == "increase")
                delta_alpha_guess = delta_alpha_guess/10;
                prev_jump = "none";
            else
                prev_jump = "decrease";
            end
        end

        if(delta_alpha_guess < 1e-10)
            disp("findAlphaReq tolerance not reached")
            break
        end
    end

    Lift = L_calc;
    Drag = D_body + D_fin;
    alpha_req = alpha_guess;
    alpha_fin = alpha_fin;
end