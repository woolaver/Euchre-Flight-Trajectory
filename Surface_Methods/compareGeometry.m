%Compare MNT CL and CD values for preliminary geometries

clear
clc
close all

HARV_CL = readmatrix("csv_files\HARV_blended_CL.csv");
HARV_CD = readmatrix("csv_files\HARV_blended_CD.csv");

Bunny_Bomb_CL = readmatrix("csv_files\Bunny_Bomb_blended_CL.csv");
Bunny_Bomb_CD = readmatrix("csv_files\Bunny_Bomb_blended_CD.csv");

Ice_Cream_Cone_CL = readmatrix("csv_files\Ice_Cream_Cone_blended_CL.csv");
Ice_Cream_Cone_CD = readmatrix("csv_files\Ice_Cream_Cone_blended_CD.csv");

Virginia_CL = readmatrix("csv_files\Virginia_blended_CL.csv");
Virginia_CD = readmatrix("csv_files\Virginia_blended_CD.csv");

Blend_Round_CL = readmatrix("csv_files\Blend_Round_blended_CL.csv");
Blend_Round_CD = readmatrix("csv_files\Blend_Round_blended_CD.csv");

Blend_Sharp_CL = readmatrix("csv_files\Blend_Sharp_blended_CL.csv");
Blend_Sharp_CD = readmatrix("csv_files\Blend_Sharp_blended_CD.csv");

Blend_Flat_CL = readmatrix("csv_files\Blend_Flat_blended_CL.csv");
Blend_Flat_CD = readmatrix("csv_files\Blend_Flat_blended_CD.csv");

Half_Cone_CL = readmatrix("csv_files\Half_Cone_blended_CL.csv");
Half_Cone_CD = readmatrix("csv_files\Half_Cone_blended_CD.csv");

%plot drag polar at Mach 5 for each
figure()
hold on
plot(Blend_Round_CD(100, :), Blend_Round_CL(100, :), 'LineWidth', 1)
plot(Blend_Sharp_CD(100, :), Blend_Sharp_CL(100, :), 'LineWidth', 1)
plot(Blend_Flat_CD(100, :), Blend_Flat_CL(100, :), 'LineWidth', 1)
plot(Half_Cone_CD(100, :), Half_Cone_CL(100, :), 'LineWidth', 1)
legend('Blend Round', 'Blend Sharp', 'Blend Flat', 'Half Cone')
title('Drag Polar at Mach 5')
xlabel('CD')
ylabel('CL')
hold off

%plot contours of L/D vs. alpha and Mach number

HARV_L_D = HARV_CL./HARV_CD;
Virginia_L_D = Virginia_CL./Virginia_CD;
Ice_Cream_Cone_L_D = Ice_Cream_Cone_CL./Ice_Cream_Cone_CD;
Bunny_Bomb_L_D = Bunny_Bomb_CL./Bunny_Bomb_CD;
Blend_Round_L_D = Blend_Round_CL./Blend_Round_CD;
Blend_Sharp_L_D = Blend_Sharp_CL ./ Blend_Sharp_CD;
Blend_Flat_L_D = Blend_Flat_CL ./ Blend_Flat_CD;
Half_Cone_L_D = Half_Cone_CL ./ Half_Cone_CD;
%{
figure()
hold on
contourf(linspace(-30, 30, 200), linspace(2, 8, 200), HARV_L_D)
colorbar
title('HARV L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off

figure()
hold on
contourf(linspace(-30, 30, 200), linspace(2, 8, 200), Virginia_L_D)
colorbar
title('Virginia L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off

figure()
hold on
contourf(linspace(-30, 30, 200), linspace(2, 8, 200), Ice_Cream_Cone_L_D)
colorbar
title('Ice Cream Cone L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off

figure()
hold on
contourf(linspace(-30, 30, 200), linspace(2, 8, 200), Bunny_Bomb_L_D)
colorbar
title('Bunny Bomb L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off

figure()
hold on
contourf(linspace(-30, 30, 200), linspace(2, 8, 200), Blend_Round_L_D)
colorbar
title('Blend Round L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off

figure()
hold on
contourf(linspace(-30, 30, 200), linspace(2, 8, 200), Blend_Sharp_L_D)
colorbar
title('Blend Sharp L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off

figure()
hold on
contourf(linspace(-30, 30, 200), linspace(2, 8, 200), Blend_Flat_L_D)
colorbar
title('Blend Flat L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off
%}

%plot CD vs. alpha at Mach 5
figure()
hold on
plot(linspace(-30, 30, 200), Blend_Round_CD(100,:), 'LineWidth', 1)
plot(linspace(-30, 30, 200), Blend_Sharp_CD(100,:), 'LineWidth', 1)
plot(linspace(-30, 30, 200), Blend_Flat_CD(100,:), 'LineWidth', 1)
plot(linspace(-30, 30, 200), Half_Cone_CD(100, :), 'LineWidth', 1)
legend('Blend Round', 'Blend Sharp', 'Blend Flat', 'Half Cone')
title('C_D vs. Angle of Attack at Mach 5')
xlabel('Alpha (deg)')
ylabel('C_D')
hold off

%plot CL vs. alpha at Mach 5
figure()
hold on
plot(linspace(-30, 30, 200), Blend_Round_CL(100,:), 'LineWidth', 1)
plot(linspace(-30, 30, 200), Blend_Sharp_CL(100,:), 'LineWidth', 1)
plot(linspace(-30, 30, 200), Blend_Flat_CL(100,:), 'LineWidth', 1)
plot(linspace(-30, 30, 200), Half_Cone_CL(100, :), 'LineWidth', 1)
legend('Blend Round', 'Blend Sharp', 'Blend Flat', 'Half Cone')
title('C_L vs. Angle of Attack at Mach 5')
xlabel('Alpha (deg)')
ylabel('C_L')
hold off

%plot L/D vs. alpha at Mach 5
figure()
hold on
plot(linspace(-30, 30, 200), Blend_Round_L_D(100,:), 'LineWidth', 1)
plot(linspace(-30, 30, 200), Blend_Sharp_L_D(100,:), 'LineWidth', 1)
plot(linspace(-30, 30, 200), Blend_Flat_L_D(100,:), 'LineWidth', 1)
plot(linspace(-30, 30, 200), Half_Cone_L_D(100,:), 'LineWidth', 1)
legend('Blend Round', 'Blend Sharp', 'Blend Flat', 'Half Cone')
title('L/D vs. Angle of Attack at Mach 5')
xlabel('Alpha (deg)')
ylabel('L/D')
hold off