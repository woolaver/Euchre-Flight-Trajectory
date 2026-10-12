% Energy accounting for the existing ODE_Model_2D trajectory.
% Run ODE_Model_2D first, then run trajectory_energy_losses.
% This script reads phase_names, phase_times, and phase_states without
% rerunning or changing the trajectory. Results are stored in energy_report.
%
% Specific mechanical energy: e = V^2/2 + g*h [J/kg].
% Total mechanical energy:    E = m*e [J].
% For this constant-mass, unpowered model, de/dt = -D*V/m.
% Thus start energy minus end energy measures modeled drag loss, including
% numerical error. A speed decrease alone is not an energy loss: it may
% represent conversion of kinetic energy into altitude.

if ~exist('phase_names', 'var') || ~exist('phase_times', 'var') || ...
        ~exist('phase_states', 'var')
    error('Energy:MissingTrajectory', ...
          'Run ODE_Model_2D first to create the phase histories.');
end

energy_report = struct();
energy_report.mass_kg = 120; % Must match the ODE functions.
energy_report.gravity_mps2 = 9.8;
energy_report.phase_count = numel(phase_names);
assert(numel(phase_times) == energy_report.phase_count && ...
       numel(phase_states) == energy_report.phase_count, ...
       'Phase names, times, and states must have matching lengths.');

% Columns: duration, range, starting energy, ending energy, energy loss,
% specific energy loss, and energy loss per meter of horizontal travel.
energy_report.values = nan(energy_report.phase_count, 7);
energy_report.history = cell(energy_report.phase_count, 1);

for energy_phase = 1:energy_report.phase_count
    energy_state = phase_states{energy_phase};
    energy_time = phase_times{energy_phase}(:);
    if isempty(energy_state) && isempty(energy_time)
        continue % An absent phase has no measurable endpoints.
    end
    assert(size(energy_state,2) == 4 && ...
           size(energy_state,1) == numel(energy_time) && ...
           isreal(energy_state) && all(isfinite(energy_state), 'all') && ...
           isreal(energy_time) && all(isfinite(energy_time)) && ...
           all(diff(energy_time) > 0), ...
           'Each phase needs finite real [V, gamma, h, x] states and increasing times.');

    energy_specific = 0.5*energy_state(:,1).^2 + ...
                      energy_report.gravity_mps2*energy_state(:,3);
    energy_total = energy_report.mass_kg*energy_specific;
    energy_range = energy_state(end,4) - energy_state(1,4);
    energy_loss = energy_total(1) - energy_total(end);
    energy_per_meter = NaN;
    if energy_range > 1e-6
        energy_per_meter = energy_loss/energy_range;
    end
    energy_report.values(energy_phase,:) = [ ...
        energy_time(end)-energy_time(1), energy_range/1000, ...
        energy_total(1)/1e6, energy_total(end)/1e6, energy_loss/1e6, ...
        (energy_specific(1)-energy_specific(end))/1000, energy_per_meter];
    energy_report.history{energy_phase} = table(energy_time, ...
        energy_specific, energy_total, ...
        'VariableNames', {'Time_s', 'SpecificEnergy_J_per_kg', 'Energy_J'});
end

energy_report.summary = array2table(energy_report.values, ...
    'VariableNames', {'Duration_s', 'Range_km', 'StartEnergy_MJ', ...
    'EndEnergy_MJ', 'Loss_MJ', 'SpecificLoss_kJ_per_kg', 'Loss_J_per_m'}, ...
    'RowNames', cellstr(string(phase_names(:))));
disp(energy_report.summary);
fprintf('Summed phase energy loss: %.3f MJ\n', ...
    sum(energy_report.summary.Loss_MJ, 'omitnan'));
fprintf('Summed phase horizontal distance: %.3f km\n', ...
    sum(energy_report.summary.Range_km, 'omitnan'));
fprintf(['Loss = starting minus ending mechanical energy; positive means loss.\n' ...
    'NaN loss/m means no forward distance; a one-sample skipped hold has zero loss.\n' ...
    'Phase sums exclude any energy or position jumps between phases.\n']);

figure('Name', 'Trajectory Energy Losses');
tiledlayout(2,1);
nexttile;
hold on
for energy_phase = 1:energy_report.phase_count
    energy_history = energy_report.history{energy_phase};
    if ~isempty(energy_history)
        plot(energy_history.Time_s, energy_history.Energy_J/1e6, ...
            '-o', 'MarkerIndices', size(energy_history,1), ...
            'DisplayName', phase_names{energy_phase});
    end
end
xlabel('Time (s)');
ylabel('Mechanical energy (MJ)');
title('Kinetic plus gravitational potential energy');
legend('show', 'Location', 'best');
grid on
hold off

nexttile;
bar(1:energy_report.phase_count, energy_report.summary.Loss_MJ);
xticks(1:energy_report.phase_count);
xticklabels(phase_names);
ylabel('Energy lost (MJ)');
title('Mechanical energy loss by phase');
grid on

% Locate losses within climb using the saved trajectory, not another ODE run.
energy_climb_index = find(strcmpi(string(phase_names), 'Climb'), 1);
if ~isempty(energy_climb_index) && ...
        size(phase_states{energy_climb_index},1) >= 2
    energy_state = phase_states{energy_climb_index};
    energy_time = phase_times{energy_climb_index}(:);
    energy_total = energy_report.history{energy_climb_index}.Energy_J;
    energy_cumulative = (energy_total(1)-energy_total)/1e6;
    % Interval-average power lost [MJ/s = MW]; avoids dividing by vertical
    % speed, which approaches zero at the apex. Negative values are retained.
    energy_rate = -diff(energy_total)./diff(energy_time)/1e6;
    energy_mid_altitude = (energy_state(1:end-1,3)+energy_state(2:end,3))/2000;
    energy_report.climb = table(energy_time, energy_state(:,3)/1000, ...
        energy_cumulative, 'VariableNames', ...
        {'Time_s', 'Altitude_km', 'CumulativeLoss_MJ'});
    energy_report.climb_intervals = table( ...
        (energy_time(1:end-1)+energy_time(2:end))/2, ...
        energy_mid_altitude, energy_rate, 'VariableNames', ...
        {'MidTime_s', 'MidAltitude_km', 'AverageLossRate_MW'});

    figure('Name', 'Where Climb Energy Is Lost');
    tiledlayout(2,1);
    nexttile;
    plot(energy_state(:,3)/1000, energy_cumulative, 'LineWidth', 1.5);
    xlabel('Altitude (km)');
    ylabel('Cumulative energy lost (MJ)');
    title('Energy lost since launch during climb');
    grid on
    hold on
    % Mark the first crossing of 50% and 90% of total climb loss.
    if energy_cumulative(end) > 0
        for energy_fraction = [0.5 0.9]
            energy_target = energy_fraction*energy_cumulative(end);
            energy_cross = find(energy_cumulative >= energy_target, 1);
            energy_weight = (energy_target-energy_cumulative(energy_cross-1))/ ...
                (energy_cumulative(energy_cross)-energy_cumulative(energy_cross-1));
            energy_altitude = energy_state(energy_cross-1,3) + energy_weight* ...
                (energy_state(energy_cross,3)-energy_state(energy_cross-1,3));
            energy_when = energy_time(energy_cross-1) + energy_weight* ...
                (energy_time(energy_cross)-energy_time(energy_cross-1));
            plot(energy_altitude/1000, energy_target, 'o', 'LineWidth', 2);
            text(energy_altitude/1000, energy_target, ...
                sprintf('  %.0f%% at %.2f km', 100*energy_fraction, energy_altitude/1000), ...
                'VerticalAlignment', 'bottom');
            fprintf('%.0f%% of climb energy loss accumulated by %.2f km at %.2f s.\n', ...
                100*energy_fraction, energy_altitude/1000, energy_when);
        end
    end
    hold off

    nexttile;
    plot(energy_mid_altitude, energy_rate, 'LineWidth', 1.5);
    xlabel('Altitude (km)');
    ylabel('Average energy loss rate (MW)');
    title('Energy loss per second over each saved time interval');
    grid on
end

clear energy_climb_index energy_cumulative energy_rate energy_mid_altitude ...
      energy_fraction energy_target energy_cross energy_weight energy_altitude energy_when
clear energy_phase energy_state energy_time energy_specific energy_total ...
      energy_range energy_loss energy_per_meter energy_history
