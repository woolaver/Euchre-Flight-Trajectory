function [alpha, Mach, CL, CD] = getTrajectoryData(trajectoryFile, climb_end, cruise_end)
% Uses trajectory data csv file to recreate local angles of attack, Mach
% number, CL, CD
    %   Inputs: trajectoryFile: csv file with trajectory states
    %           climb_end: timepoint where climb ends
    %           cruise_end: timepoint where cruise ends

    %   Because ODE45 sucks it cannot output any values other than the states   
    %   given, therefore we have to use the trajectory data to recreate each
    %   local state 

    trajectory = readmatrix(trajectoryFile);

    alpha_vec = deg2rad(linspace(-30, 30, 100));
    Mach_vec = linspace(2, 8, 100);

    CL_vec = readmatrix("HARV_CL_50_500.csv");
    CD_vec = readmatrix("HARV_CD_50_500.csv");

    m = 120; %kg
    g = 9.8; %m/s^2
    S_ref = (.2^2)*pi;

    time = trajectory(1, :);
    velocity = trajectory(2, :);
    flightPathAngle = trajectory(3, :);
    height = trajectory(4, :);

    numPoints = length(time);

    alpha = zeros([numPoints, 1]);
    Mach = zeros([numPoints, 1]);
    CL = zeros([numPoints, 1]);
    CD = zeros([numPoints, 1]);

    for i = 1:numPoints
        [~, a, ~, rho] = atmoscoesa(height(i));
        Mach(i) = velocity(i)/a;

        if (time(i) <= climb_end)
            alpha(i) = deg2rad(0);

        elseif (time(i) <= cruise_end)
            if(height(i) < 27500)
                alpha(i) = optimalAlphaMNT(Mach(i));
            else
       
                alpha(i) = deg2rad(-3);
            end

        else
            % CD required to make V_dot = 0
            CD_req = -2*m*g*sin(flightPathAngle(i))/(rho*velocity(i)^2*S_ref);
        
            % Minimum achievable drag coefficient
            alpha_zero = 0;
            CD_min = interp2(alpha_vec, Mach_vec, CD_vec, alpha_zero, Mach(i), 'linear');
            if(height(i) < 27500)
                if CD_req >= CD_min
                    % Constant velocity is physically achievable
                    alpha(i) = findAlphaMNT(CD_req, Mach(i));
                else
                    alpha(i) = 0;
                end
            else
                alpha(i) = deg2rad(-3);
            end
        end

        CL(i) = interp2(alpha_vec, Mach_vec, CL_vec, alpha(i), Mach(i), 'linear');
        CD(i) = interp2(alpha_vec, Mach_vec, CD_vec, alpha(i), Mach(i), 'linear');
    end
end
