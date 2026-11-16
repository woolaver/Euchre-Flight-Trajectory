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

%LE sweep ~ 6 deg
%LE blunt = .01 m ~ .39 inch

%% Laminar and Turbulent Assumption
%{
    
    %assming calorically perfect gas
    gamma = 1.4;
    R = 287; %J/kg*K
    cp = 1005; %J/kg*K
    Pr = .75;
    r = sqrt(Pr);

    Re_T = 1e6; %dummy value for now, need to find way to get this

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

    %if upper angle is larger than 2.5 then calculate values, if not then
    %assume friction force and heat flux is zero
    %need to set to 2.5 for interp2 to handle cone table values
    if(upper_angle > 2.5)
        %get edge states
        upper_p_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_p_pinf, upper_angle, M_inf, 'linear');
        upper_T_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_T_Tinf, upper_angle, M_inf, 'linear');
        upper_rho_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_rho_rhoinf, upper_angle, M_inf, 'linear');
        M2_u = interp2(cone_angle_vec, cone_mach_vec, cone_M2, upper_angle, M_inf, 'linear');

        pu_e = upper_p_ratio*p_inf;
        Tu_e = upper_T_ratio*T_inf;
        rhou_e = upper_rho_ratio*rho_inf;
        mu_e = calcSutherlands(Tu_e);
        a_u = sqrt(gamma*R*Tu_e);
        uu_e = M2_u*a_u;        

        Rex_u = (rhou_e.*uu_e.*x)./mu_e;

        x_T_u = find(Rex_u > Re_T, 1, 'first');

        Tw_Tu = Tw/Tu_e;

        Tstar_u_lam = Tu_e*(.45 + .55*Tw_Tu + .16*r*((gamma-1)/2)*M2_u^2);
        Tstar_u_turb = Tu_e*(.5*(1 + Tw_Tu) + .16*r*((gamma-1)/2)*M2_u^2);

        rhostar_u_lam = pu_e/(R*Tstar_u_lam);
        rhostar_u_turb = pu_e/(R*Tstar_u_turb);

        mustar_u_lam = calcSutherlands(Tstar_u_lam);
        mustar_u_turb = calcSutherlands(Tstar_u_turb);
    
        Rex_star_u_lam = (rhostar_u_lam.*uu_e.*x(1:x_T_u - 1))./mustar_u_lam;
        Rex_star_u_turb = (rhostar_u_turb.*uu_e.*x(x_T_u:end))./mustar_u_turb;

        Rex_star_u = [Rex_star_u_lam, Rex_star_u_turb];
    
        cf_star_u_lam  = .664./sqrt(Rex_star_u_lam);
        cf_star_u_turb = .02296./(Rex_star_u_turb.^.139);
        cf_star_u = [cf_star_u_lam, cf_star_u_turb];
    
        q_inf_star_u_lam = .5*rhostar_u_lam*uu_e^2;
        q_inf_star_u_turb = .5*rhostar_u_turb*uu_e^2;
    
        dx = x(2) - x(1);
        disk_area_x = r_x.*pi.*dx; %Note: half the disk area because upper and lower surfaces add
        
        friction_u_lam = sqrt(3)*sum(cf_star_u(2:x_T_u - 1).*disk_area_x(2:x_T_u-1).*q_inf_star_u_lam);
        friction_u_turb = sqrt(3)*sum(cf_star_u(x_T_u:end).*disk_area_x(x_T_u:end).*q_inf_star_u_turb);

        friction_u = friction_u_lam + friction_u_turb;
    
        %heat transfer
        he_u = cp*Tu_e;
        h_u = cp*Tw;
        haw_u = he_u + r*(uu_e^2/2);
    
        CH_star_u = (.332./sqrt(Rex_star_u)).*(Pr^(-2/3));
    
        qw_star_u_lam = CH_star_u.*rhostar_u_lam.*uu_e.*(haw_u - h_u);
        qw_star_u_turb = CH_star_u.*rhostar_u_turb.*uu_e.*(haw_u - h_u);
    
        Q_star_u_lam = sqrt(3)*sum(qw_star_u_lam(2:x_T_u-1).*disk_area_x(2:x_T_u-1));
        Q_star_u_turb = sqrt(3)*sum(qw_star_u_turb(x_T_u:end).*disk_area_x(x_T_u:end));

        Q_star_u = Q_star_u_lam + Q_star_u_turb;
    else
        friction_u = 0;
        Q_star_u = 0;
    end

    %lower surface
    %get edge states
    if (lower_angle > 2.5)
        lower_p_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_p_pinf, lower_angle, M_inf, 'linear');
        lower_T_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_T_Tinf, lower_angle, M_inf, 'linear');
        lower_rho_ratio = interp2(cone_angle_vec, cone_mach_vec, cone_rho_rhoinf, lower_angle, M_inf, 'linear');
        M2_l = interp2(cone_angle_vec, cone_mach_vec, cone_M2, lower_angle, M_inf, 'linear');
    
        pl_e = lower_p_ratio*p_inf;
        Tl_e = lower_T_ratio*T_inf;
        rhol_e = lower_rho_ratio*rho_inf;
        mul_e = calcSutherlands(Tl_e);
        a_l = sqrt(gamma*R*Tl_e);
        ul_e = M2_l*a_l;

        Rex_l = (rhol_e.*ul_e.*x)./mul_e;

        x_T_l = find(Rex_l > Re_T, 1, 'first');

        Tw_Tl = Tw/Tl_e;

        Tstar_l_lam = Tl_e*(.45 + .55*Tw_Tl + .16*r*((gamma-1)/2)*M2_l^2);
        Tstar_l_turb = Tl_e*(.5*(1 + Tw_Tl) + .16*r*((gamma-1)/2)*M2_l^2);

        rhostar_l_lam = pl_e/(R*Tstar_l_lam);
        rhostar_l_turb = pl_e/(R*Tstar_l_turb);

        mustar_l_lam = calcSutherlands(Tstar_l_lam);
        mustar_l_turb = calcSutherlands(Tstar_l_turb);
    
        Rex_star_l_lam = (rhostar_l_lam.*ul_e.*x(1:x_T_l - 1))./mustar_l_lam;
        Rex_star_l_turb = (rhostar_l_turb.*ul_e.*x(x_T_l:end))./mustar_l_turb;

        Rex_star_l = [Rex_star_l_lam, Rex_star_l_turb];
    
        cf_star_l_lam  = .664./sqrt(Rex_star_l_lam);
        cf_star_l_turb = .02296./(Rex_star_l_turb.^.139);
        cf_star_l = [cf_star_l_lam, cf_star_l_turb];
    
        q_inf_star_l_lam = .5*rhostar_l_lam*ul_e^2;
        q_inf_star_l_turb = .5*rhostar_l_turb*ul_e^2;
    
        dx = x(2) - x(1);
        disk_area_x = r_x.*pi.*dx; %Note: half the disk area because upper and lower surfaces add
        
        friction_l_lam = sqrt(3)*sum(cf_star_l(2:x_T_l - 1).*disk_area_x(2:x_T_l-1).*q_inf_star_l_lam);
        friction_l_turb = sqrt(3)*sum(cf_star_l(x_T_l:end).*disk_area_x(x_T_l:end).*q_inf_star_l_turb);

        friction_l = friction_l_lam + friction_l_turb;
    
        %heat transfer
        he_l = cp*Tl_e;
        h_l = cp*Tw;
        haw_l = he_l + r*(ul_e^2/2);
    
        CH_star_l = (.332./sqrt(Rex_star_l)).*(Pr^(-2/3));
    
        qw_star_l_lam = CH_star_l.*rhostar_l_lam.*ul_e.*(haw_l - h_l);
        qw_star_l_turb = CH_star_l.*rhostar_l_turb.*ul_e.*(haw_l - h_l);
    
        Q_star_l_lam = sqrt(3)*sum(qw_star_l_lam(2:x_T_l-1).*disk_area_x(2:x_T_l-1));
        Q_star_l_turb = sqrt(3)*sum(qw_star_l_turb(x_T_l:end).*disk_area_x(x_T_l:end));

        Q_star_l = Q_star_l_lam + Q_star_l_turb;
    else
        friction_l = 0;
        Q_star_l = 0;
    end

    Ffx = friction_l*cos(deg2rad(lower_angle)) + friction_u*cos(deg2rad(upper_angle));
    Ffy = -friction_l*sin(deg2rad(lower_angle)) + friction_u*sin(deg2rad(upper_angle));

    Q_total = Q_star_u + Q_star_l;

%}
    %% Pure Laminar Assumption
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
    if(upper_angle > 2.5)
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
    if (lower_angle > 2.5)
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
    
        %heat transfer
        he_l = cp*Tl_e;
        h_l = cp*Tw;
        haw_l = he_l + r*(ul_e^2/2);
    
        CH_star_l = (.332./sqrt(Rex_star_l)).*(Pr^(-2/3));
    
        qw_star_l = CH_star_l.*rhostar_l.*ul_e.*(haw_l - h_l);
    
        Q_star_l = sqrt(3)*sum(qw_star_l(2:end).*disk_area_x(2:end));
    else
        friction_l = 0;
        Q_star_l = 0;
    end

    Ffx = friction_l*cos(deg2rad(lower_angle)) + friction_u*cos(deg2rad(upper_angle));
    Ffy = -friction_l*sin(deg2rad(lower_angle)) + friction_u*sin(deg2rad(upper_angle));

    Q_total = Q_star_u + Q_star_l;

end
