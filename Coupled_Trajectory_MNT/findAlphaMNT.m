function alpha_required = findAlphaMNT(CD_req, Cp_max)
    %Finds the optimal alpha for best L/D for modified newtonian theory
        
    N = 100;
    alpha = deg2rad(linspace(0, 40, N));

    geometry = "HARV.STL";
    S_ref = .029; 

    for i = 1:N
        [~, CD, ~, ~, ~, ~, ~] = newtonianCLCD3D(geometry, alpha(i), Cp_max, S_ref);
        if (abs(CD - CD_req) < tol)
            alpha_required = alpha(i);
            return
        end
    end
end