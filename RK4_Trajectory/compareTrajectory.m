%Script to compare multiple trajectories using pre-saved trajectory csv
%files

clear
clc
close all

names = ["Blend_Round", "Blend_Flat", "Blend_Sharp", "Half_Cone"];
labels = ["Blend Round", "Blend Flat", "Blend Sharp", "Half Cone"];
files = "csv_files/" + names + "_Trajectory.csv";

traj = cell(1, numel(files));
for i = 1:numel(files)
    traj{i} = readtable(files(i));
end

% Each row: {column name, title, y label, convert to degrees?}
plots = {
    "Velocity",       "Velocity vs. Time",        "Velocity (m/s)",        false
    "FlightPathAngle","Flight Path Angle vs. Time","Flight Path Angle (deg)", true
    "Height",         "Height vs. Time",          "Height (m)",            false
    "Distance",       "Distance vs. Time",        "Distance (m)",          false
    "AngleOfAttack",  "Angle of Attack vs. Time", "Angle of Attack (deg)", true
    "CL",             "CL vs. Time",              "CL",                    false
    "CD",             "CD vs. Time",              "CD",                    false
    "Density",        "Density vs. Time",         "Density (kg/m^3)",      false
    "MachNumber",     "Mach Number vs. Time",     "Mach Number",           false
    };

for p = 1:size(plots, 1)
    figure()
    hold on
    for i = 1:numel(traj)
        y = traj{i}.(plots{p, 1});
        if plots{p, 4}
            y = rad2deg(y);
        end
        plot(traj{i}.Time, y, 'LineWidth', 2)
    end
    legend(labels)
    title(plots{p, 2})
    xlabel("Time (s)")
    ylabel(plots{p, 3})
    hold off
end

figure()
hold on
for i = 1:numel(traj)
    plot(traj{i}.Distance, traj{i}.Height, 'LineWidth', 2)
end
legend(labels)
title('Trajectory')
xlabel("Distance (m)")
ylabel("Height (m)")
hold off