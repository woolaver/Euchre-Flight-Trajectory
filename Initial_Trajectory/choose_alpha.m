function [alpha] = choose_alpha(V, flight_path_angle, h_dot, Mach, height)
% Function to control angle of attack logic during trajectory analysis
% NOTE: THIS FUNCTION IS NOT CURRENTLY USED BY ANYTHING BUT I DONT WANT TO
% DELETE IT IN CASE I GO BACK TO IT
% Inputs: V: projectile velocity
%         h_dot: projectile rate of change of altitude
%         flight_phase: phase of flight (0 = initial climb, ends once h_dot
%         first becomes negative, 1 = cruise, hold alpha at optimal L/D
%         value, 2 = hold Mach 3 to ensure collision speed is in
%         requirements (might look into doing a full nose dive later if
%         this comes out to having a greater range)
% Outputs: alpha: angle of attack
%          flight_phase: output flight_phase to ode for storage
 

    if (h_dot < 0 && flight_phase == 0)
        %empircal formula for CL/CD for testing purposes
        CD0 = .05;
        gamma = 1.4;
        alpha = deg2rad(linspace(0, 40, 1000));
        Cp_max = (2/(gamma*Mach^2))*((((gamma+1)/2)*Mach^2)^((gamma)/(gamma-1))*((gamma + 1)/(2*gamma*Mach^2 - (gamma - 1)))^(1/(gamma-1))-1);
        CL_test = Cp_max.*sin(alpha).^2.*cos(alpha);
        CD_test = Cp_max.*sin(alpha).^3 + CD0;

        L_D = CL_test ./ CD_test; % Calculate lift-to-drag ratio
        [~, index] = max(L_D); % Find maximum L/D and its index
        optimalAlpha = alpha(index); % Get the angle of attack for maximum L/D

        alpha = optimalAlpha;
    else
        alpha = 0;
    end

    [~, ~, ~, rho] = atmoscoesa(height);
    ground_speed = 3*340;

    %{
    if (flight_phase == 1 && V < ground_speed)
        num = -2*m*g*sin(flight_path_angle);
        den = rho*V^2*S;
        alpha = asin(nthroot((((num/den) - CD0)/Cp_max), 3));
    end
    %}
end