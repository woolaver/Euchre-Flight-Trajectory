function alpha_optimal = optimal_LD(Mach, LD_vec)
%Finds the optimal angle of attack to fly to maximize L/D at a given Mach
%number
%Inputs:
%   Mach: current Mach number
%   LD_vec_file: csv file containing L/D values varying with Mach in each
%   row and angle of attack for each column
%Outputs:
%   alpha_optimal: angle of attack that maximizes L/D at given Mach (rad)

    Mach_vec = linspace(2, 8, 200);
    alpha_vec = deg2rad(linspace(-30, 30, 200));

    Mach_row = find(Mach_vec >= Mach, 1, "first");

    LD_vec_at_mach = LD_vec(Mach_row, :);

    [~, LD_index] = max(LD_vec_at_mach);

    alpha_optimal = alpha_vec(LD_index);
end