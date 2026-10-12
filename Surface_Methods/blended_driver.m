function [Cl,Cd,Cp, x_cp]  = blended_driver(stlfile, alpha, minf, sref_lift, sref_drag,y)
    
    [Cl_s, Cd_s, ~, ~, ~, Cp_s, theta_s, ~, x_cp] = smalldist_supersonic(stlfile, alpha, minf, sref_lift, sref_drag,y);
    [Cl_n, Cd_n, ~, ~, ~, Cp_n, sinTheta] = newtonianCLCD3D(stlfile, alpha,minf, sref_lift, y);

    %interpolation
    M1 = 3;   % below small dist
    M2 = 5;   % above newtonian

    t = (minf - M1) / (M2 - M1);
    t = min(max(t, 0), 1);          % clamp

    w = 3*t.^2 - 2*t.^3;             % smoothstep:

    Cp = (1 - w).*Cp_s + w.*Cp_n;
    Cl = (1 - w).*Cl_s + w.*Cl_n;
    Cd = (1 - w).*Cd_s + w.*Cd_n;

    fprintf('Small Dist: Cl = %4f, Cd = %4f \n',Cl_s, Cd_s);
    fprintf('Newtonian: Cl = %4f, Cd = %4f \n',Cl_n, Cd_n);
    fprintf('Blended: Cl = %4f, Cd = %4f \n',Cl, Cd);

    plotVehicleCp3D(stlfile, Cp, rad2deg(theta_s) , false, [])

end