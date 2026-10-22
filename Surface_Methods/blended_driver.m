function [Cl,Cd,Cp, x_cp]  = blended_driver(stlfile, alpha, minf, sref_lift, sref_drag,y)
    
    [Cl_s, Cd_s, ~, ~, ~, Cp_s, theta_s, ~, x_cp_s] = smalldist_supersonic(stlfile, alpha, minf, sref_lift, sref_drag,y);
    [Cl_n, Cd_n, ~, ~, ~, Cp_n, ~, x_cp_n] = newtonianCLCD3D(stlfile, alpha, minf, sref_lift,sref_drag, y);

    %interpolation
    M1 = 3;   % below small dist
    M2 = 5;   % above newtonian

    t = (minf - M1) / (M2 - M1);
    t = min(max(t, 0), 1);          % clamp

    w = 3*t.^2 - 2*t.^3;             % smoothstep:

    Cp = (1 - w).*Cp_s + w.*Cp_n;
    Cl = (1 - w).*Cl_s + w.*Cl_n;
    Cd = (1 - w).*Cd_s + w.*Cd_n;
    x_cp = (1 - w).*x_cp_s + w.*x_cp_n;
    x_cp = 1- x_cp;  %this there becuase for some reason its gone from the back to the front tismallah

    %{
    fprintf('Small Dist: Cl = %4f, Cd = %4f \n',Cl_s, Cd_s);
    fprintf('Newtonian: Cl = %4f, Cd = %4f \n',Cl_n, Cd_n);
    fprintf('Blended: Cl = %4f, Cd = %4f \n',Cl, Cd);
    %}

    %plotVehicleCp3D(stlfile, Cp, rad2deg(theta_s) , false, [])

end