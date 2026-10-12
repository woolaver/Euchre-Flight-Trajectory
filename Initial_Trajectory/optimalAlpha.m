function optimal_alpha = optimalAlpha(Mach)
%Finds the optimal alpha for best L/D for empirical CL/CD formula given
%Mach number
    alpha = linspace(0, 40, 1000);

    CD0 = .05;
    gamma = 1.4;
    Cp_max = (2/(gamma*Mach^2))*((((gamma+1)/2)*Mach^2)^((gamma)/(gamma-1))*((gamma + 1)/(2*gamma*Mach^2 - (gamma - 1)))^(1/(gamma-1))-1);
    CL = Cp_max.*sind(alpha).^2.*cosd(alpha);
    CD = Cp_max.*sind(alpha).^3 + CD0;

    L_D = CL ./ CD; % Calculate lift-to-drag ratio
    [~, index] = max(L_D); % Find maximum L/D and its index
    optimal_alpha = alpha(index); % Get the angle of attack for maximum L/D
end