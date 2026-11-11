function [CL_fin, CD_fin] = calcFinForces(vehicleProperties, ps, Ms, alpha_fin, M_inf, p_inf)
%Using diamond airfoil shock-expansion theory calculates fin forces
%Inputs:
%   vehicleProperties: struct containing all necessary vehicleProperties
%   ps: post leading edge oblique shock estimation pressure (Pa)
%   Ms: post leading edge oblique shock estimation Mach Number
%   alpha_fin: local angle of attack of the fins (rad)
%   M_inf: freestream Mach Number
%   p_inf: freestream pressure (Pa)
%Outputs:
%   L_fin: fin CL
%   D_fin: fin CD

    gamma = 1.4;

    fin_length = vehicleProperties.fin_length;
    fin_height = vehicleProperties.fin_thickness/2;
    Sref_fin = vehicleProperties.Sref_fin;
    fin_width = vehicleProperties.fin_width;

    panel_length = sqrt((fin_length/2)^2 + fin_height^2);

    theta = atan(fin_height/(fin_length/2)); %half angle of diamond airfoil

    %1 upper surface leading panel, 2 is lower surface leading panel, 3 is
    %upper surface trailing panel, 4 is lower surface trailing panel

    theta1 = theta - alpha_fin;
    theta2 = theta + alpha_fin;
    theta3 = 2*theta;
    theta4 = 2*theta;

    %compression
    if(theta1 > 0)
        [p1_ps, M1, ~] = tbm(Ms, theta1);
        p1 = p1_ps*ps;
    %expansion
    else
        nu_s = pm(Ms);
        nu_1 = nu_s + abs(theta1);
        M1 = ipm(nu_1, Ms);

        ps_pt = (1 + (gamma-1)/2 * Ms^2)^(-gamma/(gamma-1));
        p1_pt = (1 + (gamma-1)/2 * M1^2)^(-gamma/(gamma-1));
        p1_ps = p1_pt/ps_pt;
        p1 = p1_ps*ps;
    end

    %compression
    if(theta2 > 0)
        [p2_ps, M2, ~] = tbm(Ms, theta2);
        p2 = p2_ps*ps;
    %expansion
    else
        nu_s = pm(Ms);
        nu_2 = nu_s + abs(theta2);
        M2 = ipm(nu_2, Ms);

        ps_pt = (1 + (gamma-1)/2 * Ms^2)^(-gamma/(gamma-1));
        p2_pt = (1 + (gamma-1)/2 * M2^2)^(-gamma/(gamma-1));
        p2_ps = p2_pt/ps_pt;
        p2 = p2_ps*ps;
    end

    %expansion over surface 3 and 4
    nu_1 = pm(M1);
    nu_3 = nu_1 + abs(theta3);
    M3 = ipm(nu_3, M1);

    p1_pt = (1 + (gamma-1)/2 * M1^2)^(-gamma/(gamma-1));
    p3_pt = (1 + (gamma-1)/2 * M3^2)^(-gamma/(gamma-1));
    p3_p1 = p3_pt/p1_pt;
    p3 = p3_p1*p1;

    nu_2 = pm(M2);
    nu_4 = nu_2 + abs(theta4);
    M4 = ipm(nu_4, M2);

    p2_pt = (1 + (gamma-1)/2 * M2^2)^(-gamma/(gamma-1));
    p4_pt = (1 + (gamma-1)/2 * M4^2)^(-gamma/(gamma-1));
    p4_p2 = p4_pt/p2_pt;
    p4 = p4_p2*p2;

    %reference to freestream to make solving for Lift and Drag easier
    p1_pinf = p1/p_inf;
    p2_pinf = p2/p_inf;
    p3_pinf = p3/p_inf;
    p4_pinf = p4/p_inf;

    Cp1 = (p1_pinf - 1)/(.5*gamma*M_inf^2);
    Cp2 = (p2_pinf - 1)/(.5*gamma*M_inf^2);
    Cp3 = (p3_pinf - 1)/(.5*gamma*M_inf^2);
    Cp4 = (p4_pinf - 1)/(.5*gamma*M_inf^2);

    CL_fin = (-Cp1*cos(theta1) + Cp2*cos(theta2) - Cp3*cos(theta2) + Cp4*cos(theta1))*(panel_length*fin_width)/(Sref_fin/2);
    CD_fin = (Cp1*sin(theta1) + Cp2*sin(theta2) + Cp3*sin(theta2) + Cp4*sin(theta1))*(panel_length*fin_width)/(Sref_fin/2);    
end