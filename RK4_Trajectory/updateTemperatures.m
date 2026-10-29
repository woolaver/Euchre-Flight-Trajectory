function [Tw_updated, T_int_updated] = updateTemperatures(vehicleProperties, alt, Tw_curr, T_int_curr, Q_conv_rate)
%Updates wall and internal temperatures given wall temperature, internal
%temperature, and aero heating heat flux
%Inputs:
%   vehicleProperties: struct containing vehicle information
%   alt: altitude (m)
%   Tw_curr: Current wall temperature (K)
%   T_int_curr: Current internal temperature (K)
%   Q_conv_rate: Convective heat transfer rate from refTempMethod (W)
%Outputs:
%   Tw_updated: new wall temperature (K)
%   T_int_updated: new internal temperature (K)

    [T_env, ~, ~, ~] = atmoscoesa(alt);

    emissivity = vehicleProperties.emissivity;
    sigma = 5.67*10^(-8); %W/m^2*K^4

    k_w = vehicleProperties.k_w;
    A = vehicleProperties.A;
    wall_thickness = vehicleProperties.wall_thickness;
    m_skin = vehicleProperties.m_skin;
    m_int = vehicleProperties.m_structure;
    dt = vehicleProperties.dt;
    cp_w = vehicleProperties.cp_w;
    cp_int = vehicleProperties.cp_int;

    %radiative heat transfer
    Q_rad_rate = emissivity*sigma*A*(T_env);
    
    %conductive heat transfer using steady state assumption for each time point
    Q_cond_rate = (k_w*A*(Tw_curr - T_int_curr))/(wall_thickness);
    
    Qw_total_rate = Q_rad_rate + Q_conv_rate - Q_cond_rate;
    
    Qw_total = Qw_total_rate * dt;
    Q_int_total = Q_cond_rate*dt;
    
    delta_Tw = Qw_total/(m_skin * cp_w);
    delta_T_int = Q_int_total/(m_int * cp_int);
    
    Tw_updated = Tw_curr + delta_Tw;
    T_int_updated = T_int_curr + delta_T_int;
end