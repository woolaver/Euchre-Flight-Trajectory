%Script to compare multiple trajectories using pre-saved trajectory csv
%files

clear
clc
close all

virginia_trajectory = readtable("Virginia_Trajectory.csv");
ice_cream_cone_trajectory = readtable("Ice_Cream_Cone_Trajectory.csv");
harv_trajectory = readtable("HARV_Trajectory.csv");
bunny_bomb_trajectory = readtable("Bunny_Bomb_Trajectory.csv");

virginia_time = virginia_trajectory.Time;
ice_cream_cone_time = ice_cream_cone_trajectory.Time;
harv_time = harv_trajectory.Time;
bunny_bomb_time = bunny_bomb_trajectory.Time;

figure()
hold on
plot(virginia_time, virginia_trajectory.Velocity, 'LineWidth', 1)
plot(ice_cream_cone_time, ice_cream_cone_trajectory.Velocity, 'LineWidth', 1)
plot(harv_time, harv_trajectory.Velocity, 'LineWidth', 1)
plot(bunny_bomb_time, bunny_bomb_trajectory.Velocity, 'LineWidth', 1)
legend("Virginia", "Ice Cream Cone", "HARV", "Bunny Bomb")
title('Velocity vs. Time')
xlabel("Time (s)")
ylabel("Velocity (m/s)")
hold off

figure()
hold on
plot(virginia_time, rad2deg(virginia_trajectory.FlightPathAngle), 'LineWidth', 1)
plot(ice_cream_cone_time, rad2deg(ice_cream_cone_trajectory.FlightPathAngle), 'LineWidth', 1)
plot(harv_time, rad2deg(harv_trajectory.FlightPathAngle), 'LineWidth', 1)
plot(bunny_bomb_time, rad2deg(bunny_bomb_trajectory.FlightPathAngle), 'LineWidth', 1)
legend("Virginia", "Ice Cream Cone", "HARV", "Bunny Bomb")
title('Flight Path Angle vs. Time')
xlabel("Time (s)")
ylabel("Flight Path Angle (deg)")
hold off

figure()
hold on
plot(virginia_time, virginia_trajectory.Height, 'LineWidth', 1)
plot(ice_cream_cone_time, ice_cream_cone_trajectory.Height, 'LineWidth', 1)
plot(harv_time, harv_trajectory.Height, 'LineWidth', 1)
plot(bunny_bomb_time, bunny_bomb_trajectory.Height, 'LineWidth', 1)
legend("Virginia", "Ice Cream Cone", "HARV", "Bunny Bomb")
title('Height vs. Time')
xlabel("Time (s)")
ylabel("Height (m)")
hold off

figure()
hold on
plot(virginia_time, virginia_trajectory.Distance, 'LineWidth', 1)
plot(ice_cream_cone_time, ice_cream_cone_trajectory.Distance, 'LineWidth', 1)
plot(harv_time, harv_trajectory.Distance, 'LineWidth', 1)
plot(bunny_bomb_time, bunny_bomb_trajectory.Distance, 'LineWidth', 1)
legend("Virginia", "Ice Cream Cone", "HARV", "Bunny Bomb")
title('Distance vs. Time')
xlabel("Time (s)")
ylabel("Distance (m)")
hold off

figure()
hold on
plot(virginia_time, rad2deg(virginia_trajectory.AngleOfAttack), 'LineWidth', 1)
plot(ice_cream_cone_time, rad2deg(ice_cream_cone_trajectory.AngleOfAttack), 'LineWidth', 1)
plot(harv_time, rad2deg(harv_trajectory.AngleOfAttack), 'LineWidth', 1)
plot(bunny_bomb_time, rad2deg(bunny_bomb_trajectory.AngleOfAttack), 'LineWidth', 1)
legend("Virginia", "Ice Cream Cone", "HARV", "Bunny Bomb")
title('Angle of Attack vs. Time')
xlabel("Time (s)")
ylabel("Angle of Attack (deg)")
hold off

figure()
hold on
plot(virginia_time, virginia_trajectory.CL, 'LineWidth', 1)
plot(ice_cream_cone_time, ice_cream_cone_trajectory.CL, 'LineWidth', 1)
plot(harv_time, harv_trajectory.CL, 'LineWidth', 1)
plot(bunny_bomb_time, bunny_bomb_trajectory.CL, 'LineWidth', 1)
legend("Virginia", "Ice Cream Cone", "HARV", "Bunny Bomb")
title('CL vs. Time')
xlabel("Time (s)")
ylabel("CL")
hold off

figure()
hold on
plot(virginia_time, virginia_trajectory.CD, 'LineWidth', 1)
plot(ice_cream_cone_time, ice_cream_cone_trajectory.CD, 'LineWidth', 1)
plot(harv_time, harv_trajectory.CD, 'LineWidth', 1)
plot(bunny_bomb_time, bunny_bomb_trajectory.CD, 'LineWidth', 1)
legend("Virginia", "Ice Cream Cone", "HARV", "Bunny Bomb")
title('CD vs. Time')
xlabel("Time (s)")
ylabel("CD")
hold off

figure()
hold on
plot(virginia_time, virginia_trajectory.Density, 'LineWidth', 1)
plot(ice_cream_cone_time, ice_cream_cone_trajectory.Density, 'LineWidth', 1)
plot(harv_time, harv_trajectory.Density, 'LineWidth', 1)
plot(bunny_bomb_time, bunny_bomb_trajectory.Density, 'LineWidth', 1)
legend("Virginia", "Ice Cream Cone", "HARV", "Bunny Bomb")
title('Density vs. Time')
xlabel("Time (s)")
ylabel("Density (kg/m^3)")
hold off

figure()
hold on
plot(virginia_time, virginia_trajectory.MachNumber, 'LineWidth', 1)
plot(ice_cream_cone_time, ice_cream_cone_trajectory.MachNumber, 'LineWidth', 1)
plot(harv_time, harv_trajectory.MachNumber, 'LineWidth', 1)
plot(bunny_bomb_time, bunny_bomb_trajectory.MachNumber, 'LineWidth', 1)
legend("Virginia", "Ice Cream Cone", "HARV", "Bunny Bomb")
title('Mach Number vs. Time')
xlabel("Time (s)")
ylabel("Mach Number")
hold off

figure()
hold on
plot(virginia_trajectory.Distance, virginia_trajectory.Height, 'LineWidth', 1)
plot(ice_cream_cone_trajectory.Distance, ice_cream_cone_trajectory.Height, 'LineWidth', 1)
plot(harv_trajectory.Distance, harv_trajectory.Height, 'LineWidth', 1)
plot(bunny_bomb_trajectory.Distance, bunny_bomb_trajectory.Height, 'LineWidth', 1)
legend("Virginia", "Ice Cream Cone", "HARV", "Bunny Bomb")
title('Trajectory')
xlabel("Distance (m)")
ylabel("Height (m)")
hold off