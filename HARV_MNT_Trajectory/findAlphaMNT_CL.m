function alpha_required = findAlphaMNT_CL(CL_req, Mach)
    %Finds the optimal alpha for best L/D for modified newtonian theory
        
    N = 1000;
    alpha = deg2rad(linspace(-30, 30, N));
    tol = 1e-2;

    Mach_vec = linspace(2, 8, 50);
    alpha_vec = deg2rad(linspace(-30, 30, 500));

    CL_vec = readmatrix("HARV_CL_50_500.csv");

    for i = 1:N
        CL = interp2(alpha_vec, Mach_vec, CL_vec, alpha(i), Mach, 'linear');
        if (abs(CL - CL_req) < tol)
            alpha_required = alpha(i);
            return
        end
    end
end