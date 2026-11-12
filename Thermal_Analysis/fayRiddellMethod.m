function [Q_stag_rate] = fayRiddellMethod(R, M_inf, alt, T_nose)
%Fay-Riddell Method to get stagnation point heating at nose of vehicle
%Assuming equilibrium flow with a spherical nose 
%All heat transfer will be assumed to convectivally transfer into the nose
%which will then transfer heat conductivly to the internal
%Inputs:
%   vehicleProperties: struct containing all necessary vehicle properties
%   M_inf: Freestream Mach number
%   alt: altitude (m)
%   T_nose: current nose wall temperature (K)
%Outputs:
%   Q_stag_rate: Total heat transfer rate from stagnation heating

    gamma = 1.4;


    [T_inf, ~, p_inf, rho_inf] = atmoscoesa(alt);

    %Normal Shock to get shock properties
    
    M_s = ((gamma-1)*M_inf^2 + 2)/(2*gamma*M_inf^2 - (gamma-1));
    rhos_rhoinf = ((gamma+1)*M_inf^2)/((gamma-1)*M_inf^2 + 2);
    Ts_Tinf = ((2*gamma*M_inf^2 - (gamma-1))*((gamma-1)*M_inf^2 + 2))/((gamma+1)^2*M_inf^2);
    ps_pinf = (2*gamma*M_inf^2 - (gamma-1))/(gamma+1);
    
    rho_s = rhos_rhoinf*rho_inf;
    T_s = Ts_Tinf*T_inf;
    p_s = ps_pinf*p_inf;

    %isentropic stagnation to get to BL edge state

    rhow_rhos = (1 + (gamma-1)/2 * M_s^2)^(1/(gamma-1));
    Te_Ts = (1 + (gamma-1)/2 * M_s^2);
    pe_ps = (1 + (gamma-1)/2 * M_s^2);

    rho_e = rhow_rhos*rho_s;
    T_e = Te_Ts*T_s;
    p_e = pe_ps*p_s;

    mu_e = calcSutherlands(T_e);

    due_dx = (1/R)*sqrt((2*(p_e - p_inf))/rho_e);

    %calculate wall state from current wall temperature
    cp_air = 1005; %J/kg*K
    h_w = T_nose*cp_air;

    p_w = p_e; %assumption that wall pressure is constant in boundary layer

    R = 287; %J/kg*K
    rho_w = p_w/(R*T_nose);
    mu_w = calcSutherlands(T_nose);

    h_e = cp_air*T_e;

    %Ignoring dissociation effects for now, Claude estimates ~2-4%
    %difference but need to verify this
    Le = 1.4; %air assumption
    Pr = 0.75;

    %Q_stag_rate = 0.763*Pr^(-0.6)*(rho_e*mu_e)^(.4)*(rho_w*mu_w)^(.1)*sqrt(due_dx)*(h_e - h_w)*(1 + (Le^(.52) - 1)*(hD/h_e));
    qw_stag_rate = 0.763*Pr^(-0.6)*(rho_e*mu_e)^(.4)*(rho_w*mu_w)^(.1)*sqrt(due_dx)*(h_e - h_w);

    Q_stag_rate = qw_stag_rate;
end


% 
% function [Q_stag_rate] = fayRiddellMethod(vehicleProperties, M_inf, alt, T_nose)
% %Fay-Riddell Method to get stagnation point heating at nose of vehicle
% %Assuming equilibrium flow with a spherical nose 
% %All heat transfer will be assumed to convectivally transfer into the nose
% %which will then transfer heat conductivly to the internal
% %Inputs:
% %   vehicleProperties: struct containing all necessary vehicle properties
% %   M_inf: Freestream Mach number
% %   alt: altitude (m)
% %   T_nose: current nose wall temperature (K)
% %Outputs:
% %   Q_stag_rate: Total heat transfer rate from stagnation heating
% 
%     gamma = 1.4;
% 
%     R_nose = vehicleProperties.R_nose;
% 
%     [T_inf, ~, p_inf, rho_inf] = atmoscoesa(alt);
% 
%     %Normal Shock to get shock properties
% 
%     M_s = ((gamma-1)*M_inf^2 + 2)/(2*gamma*M_inf^2 - (gamma-1));
%     rhos_rhoinf = ((gamma+1)*M_inf^2)/((gamma-1)*M_inf^2 + 2);
%     Ts_Tinf = ((2*gamma*M_inf^2 - (gamma-1))*((gamma-1)*M_inf^2 + 2))/((gamma+1)^2*M_inf^2);
%     ps_pinf = (2*gamma*M_inf^2 - (gamma-1))/(gamma+1);
% 
%     rho_s = rhos_rhoinf*rho_inf;
%     T_s = Ts_Tinf*T_inf;
%     p_s = ps_pinf*p_inf;
% 
%     %isentropic stagnation to get to BL edge state
% 
%     rhow_rhos = (1 + (gamma-1)/2 * M_s^2)^(1/(gamma-1));
%     Te_Ts = (1 + (gamma-1)/2 * M_s^2);
%     pe_ps = (1 + (gamma-1)/2 * M_s^2);
% 
%     rho_e = rhow_rhos*rho_s;
%     T_e = Te_Ts*T_s;
%     p_e = pe_ps*p_s;
% 
%     mu_e = calcSutherlands(T_e);
% 
%     due_dx = (1/R_nose)*sqrt((2*(p_e - p_inf))/rho_e);
% 
%     %calculate wall state from current wall temperature
%     cp_air = 1005; %J/kg*K
%     h_w = T_nose*cp_air;
% 
%     p_w = p_e; %assumption that wall pressure is constant in boundary layer
% 
%     R = 287; %J/kg*K
%     rho_w = p_w/(R*T_nose);
%     mu_w = calcSutherlands(T_nose);
% 
%     h_e = cp_air*T_e;
% 
%     %Ignoring dissociation effects for now, Claude estimates ~2-4%
%     %difference but need to verify this
%     Le = 1.4; %air assumption
%     Pr = 0.75;
% 
%     %Q_stag_rate = 0.763*Pr^(-0.6)*(rho_e*mu_e)^(.4)*(rho_w*mu_w)^(.1)*sqrt(due_dx)*(h_e - h_w)*(1 + (Le^(.52) - 1)*(hD/h_e));
%     qw_stag_rate = 0.763*Pr^(-0.6)*(rho_e*mu_e)^(.4)*(rho_w*mu_w)^(.1)*sqrt(due_dx)*(h_e - h_w);
% 
%     Q_stag_rate = qw_stag_rate*vehicleProperties.A_nose;
% end