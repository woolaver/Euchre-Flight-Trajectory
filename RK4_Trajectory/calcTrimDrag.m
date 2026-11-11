function [trimDrag, trimLift, alpha_fin] = calcTrimDrag(vehicleProperties, alpha_body, L_body, D_body, alt, M_inf)
%Uses diamond airfoil assumption to calculate total drag and lift created
%by canards to trim the vehicle
%Inputs:
%   vehicleProperties: Struct containing all vehicle properties
%   alpha_body: current angle of attack of the body (rad)
%   L_body: Lift produced by body (N)
%   D_body: Drag produced by body (N)
%   alt: current altitude (m)
%   M_inf: freestream Mach number
%Outputs:
%   trimDrag: additional drag from the fins to trim the vehicle (N)
%   trimLift: additional lift from the fins to trim the vehicle (N)
%   alpha_fin: angle of attack required from the fins to trim the vehicle
%   (rad)

    %calculate post leading edge shock properties, assuming constant
    %properties over the fins

    %initial assumption is simple oblique shock at the nose, should be far
    %enough downstream of nose that bow shock does not have major affect

    [~, a_inf, p_inf, rho_inf] = atmoscoesa(alt);

    V_inf = M_inf*a_inf;

    %use fin location to calculate angle to assume oblique shock at
    theta_OS = tan(vehicleProperties.fin_y_location/vehicleProperties.fin_x_location);

    [ps_pinf, M_s, ~] = tbm(M_inf, theta_OS);

    p_s = ps_pinf*p_inf;

    %find alpha_fin required to trim vehicle

    Sref_fin = vehicleProperties.Sref_fin;

    M_cg = 1000; %dummy value to get while loop started
    alpha_fin = deg2rad(0);

    fin_marm = vehicleProperties.cg_x_location - vehicleProperties.fin_x_location;
    body_marm = vehicleProperties.cp_x_location - vehicleProperties.cg_x_location;
    
    tol = 1e-3;
    delta_alpha_fin = deg2rad(1);
    prev_jump = "none";

    while abs(M_cg) > tol
        [CL_fin, CD_fin] = calcFinForces(vehicleProperties, p_s, M_s, alpha_fin, M_inf, p_inf);
        L_fin = 2*(.5*rho_inf*V_inf^2*Sref_fin*CL_fin);
        D_fin = 2*(.5*rho_inf*V_inf^2*Sref_fin*CD_fin);

        M_cg = body_marm*(cos(alpha_body)*L_body - sin(alpha_body)*D_body) + fin_marm*(-cos(alpha_body + alpha_fin)*L_fin + sin(alpha_body + alpha_fin)*D_fin);

        %if moment is positive -> vehicle pitch up -> need lower fin_alpha
        if(M_cg > tol)
            alpha_fin = alpha_fin + delta_alpha_fin;
            %if switching between increasing and decreasing then decrease
            %delta_alpha_fin
            if(prev_jump == "decrease")
                delta_alpha_fin = delta_alpha_fin/10;
                prev_jump = "none";
            else
                prev_jump = "increase";
            end
        else
            alpha_fin = alpha_fin - delta_alpha_fin;

            %if switching between increasing and decreasing then decrease
            %delta_alpha_fin
            if(prev_jump == "increase")
                delta_alpha_fin = delta_alpha_fin/10;
                prev_jump = "none";
            else
                prev_jump = "decrease";
            end
        end

        if(delta_alpha_fin < 1e-10)
            disp("Some shit happened")
            break;
        end
    end
    trimDrag = D_fin;
    trimLift = L_fin;
end