function [Ffx, Ffy, Q_total] = refTempMethod(length, height, alt, M_inf, alpha, Tw)
%Eckerts Reference Temperature Method, using an assumption of a cone with
%half angle equal to that of the vehicle nose with the same characteristic
%length and total surface area
%Inputs:
%   length: reference length of vehicle (m)
%   height: height of back edge of vehicle (m)
%   alt: altitude (m)
%   M_inf: Freestream mach number
%   alpha: angle of attack (rad)
%   Tw: wall temperature (K)
%Outputs:
%   Ffx: X-direction friction force (N)
%   Ffy: Y-direction friction force (N)
%   Q_tot: Total heat flux (W)
    
    %assming calorically perfect gas
    gamma = 1.4;
    R = 287; %J/kg*K
    cp = 1005; %J/kg*K
    Pr = .75;
    r = sqrt(Pr);

    %read in cone tables
    cone_p_pinf = readmatrix('coneTables\cone_pw_pinf.csv');
    cone_T_Tinf = readmatrix('coneTables\cone_Tw_Tinf.csv');
    cone_M2 = readmatrix('coneTables\cone_M2.csv');
    cone_rho_rhoinf = readmatrix('coneTables\cone_rhow_rhoinf.csv');

    %ignore first row and column to remove headers
    cone_p_pinf = cone_p_pinf(2:end, 2:end);
    cone_T_Tinf = cone_T_Tinf(2:end, 2:end);
    cone_rho_rhoinf = cone_rho_rhoinf(2:end, 2:end);
    cone_M2 = cone_M2(2:end, 2:end);

    %vectors for interpolation
    cone_mach_vec = [2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10];
    cone_angle_vec = [2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30]; %deg

    % Estimate cone geometry and freestream atmospheric properties
    half_angle = atan(height/length);
    [T_inf, a_inf, p_inf, rho_inf] = atmoscoesa(alt);
    V_inf = M_inf * a_inf;
    mu_inf = calcSutherlands(T_inf);
    
    length_cone = length/cos(half_angle);

    x = linspace(0, length_cone, 1000);
    r_x = sin(half_angle).*x;
    
    upper_angle = rad2deg(half_angle - alpha);
    lower_angle = rad2deg(half_angle + alpha);

    %if upper angle is positive then calcualte values, if negative then
    %assume friction force and heat flux is zero
    if(upper_angle > 0)
        %get edge states
        upper_p_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_p_pinf, upper_angle, M_inf, 'linear');
        upper_T_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_T_Tinf, upper_angle, M_inf, 'linear');
        upper_rho_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_rho_rhoinf, upper_angle, M_inf, 'linear');
        M2_u = interp2(cone_angle_vec, cone_mach_vec, cone_M2, upper_angle, M_inf, 'linear');

        pu_e = upper_p_ratio*p_inf;
        Tu_e = upper_T_ratio*T_inf;
        rhou_e = upper_rho_ratio*rho_inf;
        mu_e = calcSutherlands(Tu_e);

        Tw_Tu = Tw/Tu_e;
        Tstar_u = Tu_e*(1 + .032*M2_u^2 + .58*(Tw_Tu - 1));
    
        rhostar_u = pu_e/(R*Tstar_u);
        mustar_u = calcSutherlands(Tstar_u);
    
        a_u = sqrt(gamma*R*Tu_e);
        uu_e = M2_u*a_u;
    
        Rex_star_u = (rhostar_u.*uu_e.*x)./mustar_u;
    
        %assuming laminar flow for now, will try later with turbulent flow to
        %see if there are major differences
        cf_star_u = .664./(sqrt(Rex_star_u));
    
        q_inf_star_u = .5*rhostar_u*uu_e^2;
    
        dx = x(2) - x(1);
        disk_area_x = r_x.*pi.*dx; %Note: half the disk area because upper and lower surfaces add
    
        friction_u = sqrt(3)*sum(cf_star_u(2:end).*disk_area_x(2:end).*q_inf_star_u);
    
        %heat transfer
        he_u = cp*Tu_e;
        h_u = cp*Tw;
        haw_u = he_u + r*(uu_e^2/2);
    
        CH_star_u = (.332./sqrt(Rex_star_u)).*(Pr^(-2/3));
    
        qw_star_u = CH_star_u.*rhostar_u.*uu_e.*(haw_u - h_u);
    
        Q_star_u = sqrt(3)*sum(qw_star_u(2:end).*disk_area_x(2:end));

    else
        friction_u = 0;
        Q_star_u = 0;
    end

    %lower surface
    %get edge states
    lower_p_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_p_pinf, lower_angle, M_inf, 'linear');
    lower_T_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_T_Tinf, lower_angle, M_inf, 'linear');
    lower_rho_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_rho_rhoinf, lower_angle, M_inf, 'linear');
    M2_l = interp2(cone_angle_vec, cone_mach_vec, cone_M2, lower_angle, M_inf, 'linear');

    pl_e = lower_p_ratio*p_inf;
    Tl_e = lower_T_ratio*T_inf;
    rhol_e = lower_rho_ratio*rho_inf;
    ml_e = calcSutherlands(Tl_e);

    Tw_Tl = Tw/Tl_e;
    Tstar_l = Tl_e*(1 + .032*M2_l^2 + .58*(Tw_Tl - 1));

    rhostar_l = pl_e/(R*Tstar_l);
    mustar_l = calcSutherlands(Tstar_l);

    a_l = sqrt(gamma*R*Tl_e);
    ul_e = M2_l*a_l;

    Rex_star_l = (rhostar_l.*ul_e.*x)./mustar_l;

    %assuming laminar flow for now, will try later with turbulent flow to
    %see if there are major differences
    cf_star_l = .664./(sqrt(Rex_star_l));

    q_inf_star_l = .5*rhostar_l*ul_e^2;

    dx = x(2) - x(1);
    disk_area_x = r_x.*pi.*dx; %Note: half the disk area because upper and lower surfaces added

    friction_l = sqrt(3)*sum(cf_star_l(2:end).*disk_area_x(2:end).*q_inf_star_l);

    Ffx = friction_l*cos(deg2rad(lower_angle)) + friction_u*cos(deg2rad(upper_angle));
    Ffy = -friction_l*sin(deg2rad(lower_angle)) + friction_u*sin(deg2rad(upper_angle));

    %heat transfer
    he_l = cp*Tl_e;
    h_l = cp*Tw;
    haw_l = he_l + r*(ul_e^2/2);

    CH_star_l = (.332./sqrt(Rex_star_l)).*(Pr^(-2/3));

    qw_star_l = CH_star_l.*rhostar_l.*ul_e.*(haw_l - h_l);

    Q_star_l = sqrt(3)*sum(qw_star_l(2:end).*disk_area_x(2:end));

    Q_total = Q_star_l + Q_star_u;
end