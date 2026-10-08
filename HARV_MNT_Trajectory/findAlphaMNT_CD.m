function alpha_required = findAlphaMNT_CD(CD_req, Mach)
    %Finds the optimal alpha for best L/D for modified newtonian theory
        
    N = 100;
    alpha = deg2rad(linspace(-30, 30, N));
    tol = 1e-2;

    Mach_vec = linspace(2, 8, 50);
    alpha_vec = deg2rad(linspace(-30, 30, 500));

    CD_vec = readmatrix("HARV_CD_50_500.csv");

    alpha_best = optimalAlphaMNT(Mach);

    CD_best = interp2(alpha_vec, Mach_vec, CD_vec, alpha_best, Mach, 'linear');

    if (CD_best < CD_req)
        alpha_required = alpha_best;
    else
        for i = 1:N
            CD = interp2(alpha_vec, Mach_vec, CD_vec, alpha(i), Mach, 'linear');
            if (abs(CD - CD_req) < tol)
                alpha_required = alpha(i);
                return
            end
        end
    end
end