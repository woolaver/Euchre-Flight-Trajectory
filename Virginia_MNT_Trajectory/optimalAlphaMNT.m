function optimal_alpha = optimalAlphaMNT(Mach)
    %Finds the optimal alpha for best L/D for modified newtonian theory
    
    N = 100;
    alpha = deg2rad(linspace(-30, 30, N));
    
    %alpha_vec and Mach_vec must be same as the ones used to create CL_vec
    %and CD_vec
    alpha_vec = deg2rad(linspace(-30, 30, 100));
    Mach_vec = linspace(2, 8, 100);

    CL_vec = readmatrix("Virginia_CL.csv");
    CD_vec = readmatrix("Virginia_CD.csv");

    CL_at_alpha = zeros([1, N]);
    CD_at_alpha = zeros([1, N]);

    for i = 1:N
        CL_at_alpha(i) = interp2(alpha_vec, Mach_vec, CL_vec, alpha(i), Mach, 'linear');
        CD_at_alpha(i) = interp2(alpha_vec, Mach_vec, CD_vec, alpha(i), Mach, 'linear');    
    end

    L_D = CL_at_alpha ./ CD_at_alpha; % Calculate lift-to-drag ratio
    [~, index] = max(L_D); % Find maximum L/D and its index
    optimal_alpha = alpha(index); % Get the angle of attack for maximum L/D
end