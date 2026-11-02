function [Tw_updated, T_int_updated, T_nose_updated] = updateTemperatures(vehicleProperties, alt, Tw_curr, T_int_curr, T_nose_curr, Q_conv_rate, Q_conv_nose_rate)
%Updates wall and internal temperatures given wall temperature, internal
%temperature, and aero heating heat flux
%Inputs:
%   vehicleProperties: struct containing vehicle information
%   alt: altitude (m)
%   Tw_curr: Current wall temperature (K)
%   T_int_curr: Current internal temperature (K)
%   T_nose_curr: Current nose temperature (K)
%   Q_conv_rate: Convective heat transfer rate from refTempMethod (W)
%   Q_conv_nose_rate: Convective heat transfer rate for the nose from
%   fayRiddellMethod (W)
%Outputs:
%   Tw_updated: new wall temperature (K)
%   T_int_updated: new internal temperature (K)
%   T_nose_updated: new nose tempearture (K)

    [T_env, ~, ~, ~] = atmoscoesa(alt);

    emissivity = vehicleProperties.emissivity;
    sigma = 5.67*10^(-8); %W/m^2*K^4

    dt = vehicleProperties.dt;
    
    %Wall Properties
    k_w = vehicleProperties.k_w;
    A = vehicleProperties.A;
    wall_thickness = vehicleProperties.wall_thickness;
    m_skin = vehicleProperties.m_skin;
    cp_w = vehicleProperties.cp_w;

    %Internal Properties
    m_int = vehicleProperties.m_structure;
    cp_int = vehicleProperties.cp_int;
    
    %Nose Properties
    A_nose = vehicleProperties.A_nose;
    m_nose = vehicleProperties.m_nose;
    R_nose = vehicleProperties.R_nose;
    k_nose = vehicleProperties.k_nose;
    D_nose = vehicleProperties.D_nose;
    cp_nose = vehicleProperties.cp_nose;
    emissivity_nose = vehicleProperties.emissivity_nose;

    %radiative heat transfer
    Q_rad_rate = emissivity*sigma*A*(Tw_curr^4 - T_env^4);
    
    %conductive heat transfer using steady state assumption for each time point
    Q_cond_rate = (k_w*A*(Tw_curr - T_int_curr))/(wall_thickness);

    %nose conductive heat transfer
    Q_rad_nose_rate = emissivity_nose*sigma*A_nose*(T_nose_curr^4 - T_env^4);

    %conductive heat transfer between internal structure and nose
    Q_cond_nose_rate = (k_nose*D_nose*(T_nose_curr - T_int_curr))/(R_nose);
    
    Qw_total_rate = Q_conv_rate - Q_cond_rate - Q_rad_rate;
    Q_total_nose_rate = Q_conv_nose_rate - Q_cond_nose_rate - Q_rad_nose_rate;
    Q_int_total_rate = Q_cond_rate + Q_cond_nose_rate;
    
    Qw_total = Qw_total_rate * dt;
    Q_int_total = Q_int_total_rate*dt;
    Q_nose_total = Q_total_nose_rate*dt;
    
    delta_Tw = Qw_total/(m_skin * cp_w);
    delta_T_int = Q_int_total/(m_int * cp_int);
    delta_T_nose = Q_nose_total/(m_nose*cp_nose);
    
    Tw_updated = Tw_curr + delta_Tw;
    T_int_updated = T_int_curr + delta_T_int;
    T_nose_updated = T_nose_curr + delta_T_nose;
end