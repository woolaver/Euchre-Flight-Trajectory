function [Cf] = friction_drag(total_vel, alt, ref_length)
%Calculates skin friction coefficient
    
    [T, ~, ~, rho] = atmoscoesa(alt);

    %calculate kinematic viscocity with Sutherland's formula
    S = 110.4; %K
    mu_0 = 1.716*10^(-5); %Pa*s
    T_0 = 273.15; %K
    mu = mu_0*(T/T_0)^(1.5)*((T_0+S)/(T+S));

    %discretize the length of the projectile into 10^6 parts
    
    dx = ref_length/10^6;
    x_pos = 0;
    Cf = 0;
    
    Re_crit = 5*10^5;

    %integrate over length of projectile to calculate Cf
    for i = 1:10^6
        x_pos = x_pos + dx;
        Re_x = rho*total_vel*x_pos/mu;

        if Re_x < Re_crit
            Cf_local = (.664/(sqrt(Re_x)))*dx;
        else
            Cf_local = (.0576/(Re_x^(1/5)))*dx;
        end
        Cf = Cf + Cf_local;
    end
end