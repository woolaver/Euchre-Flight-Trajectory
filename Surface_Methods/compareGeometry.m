%Compare MNT CL and CD values for preliminary geometries

clear
clc
close all

HARV_CL = readmatrix("HARV_CL_1000.csv");
HARV_CD = readmatrix("HARV_CD_1000.csv");

Virginia_CL = readmatrix("Virginia_CL_1000.csv");
Virginia_CD = readmatrix("Virginia_CD_1000.csv");

Ice_Cream_Cone_CL = readmatrix("Ice_Cream_Cone_CL_1000.csv");
Ice_Cream_Cone_CD = readmatrix("Ice_Cream_Cone_CD_1000.csv");

Bunny_Bomb_CL = readmatrix("Bunny_Bomb_CL.csv");
Bunny_Bomb_CD = readmatrix("Bunny_Bomb_CD.csv");

%plot drag polar at Mach 5 for each
figure()
hold on
plot(HARV_CD(501, :), HARV_CL(501,:), 'LineWidth', 1)
plot(Virginia_CD(501, :), Virginia_CL(501, :), 'LineWidth', 1)
plot(Ice_Cream_Cone_CD(501, :), Ice_Cream_Cone_CL(501, :), 'LineWidth', 1)
plot(Bunny_Bomb_CD(51, :), Bunny_Bomb_CL(51, :), 'LineWidth', 1)
legend('HARV', 'Virginia', 'Ice Cream Cone', 'Bunny Bomb')
title('Drag Polar at Mach 5')
xlabel('CD')
ylabel('CL')
hold off

%plot contours of L/D vs. alpha and Mach number

HARV_L_D = HARV_CL./HARV_CD;
Virginia_L_D = Virginia_CL./Virginia_CD;
Ice_Cream_Cone_L_D = Ice_Cream_Cone_CL./Ice_Cream_Cone_CD;
Bunny_Bomb_L_D = Bunny_Bomb_CL./Bunny_Bomb_CD;

figure()
hold on
contourf(linspace(-30, 30, 1000), linspace(2, 8, 1000), HARV_L_D)
colorbar
title('HARV L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off

figure()
hold on
contourf(linspace(-30, 30, 1000), linspace(2, 8, 1000), Virginia_L_D)
colorbar
title('Virginia L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off

figure()
hold on
contourf(linspace(-30, 30, 1000), linspace(2, 8, 1000), Ice_Cream_Cone_L_D)
colorbar
title('Ice Cream Cone L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off

figure()
hold on
contourf(linspace(-30, 30, 100), linspace(2, 8, 100), Bunny_Bomb_L_D)
colorbar
title('Bunny Bomb L/D over Alpha and Mach Number')
xlabel('Alpha (deg)')
ylabel('Mach')
hold off

%plot CD vs. alpha at Mach 5
figure()
hold on
plot(linspace(-30, 30, 1000), HARV_CD(501,:), 'LineWidth', 1)
plot(linspace(-30, 30, 1000), Virginia_CD(501,:), 'LineWidth', 1)
plot(linspace(-30, 30, 1000), Ice_Cream_Cone_CD(501,:), 'LineWidth', 1)
plot(linspace(-30, 30, 100), Bunny_Bomb_CD(51,:), 'LineWidth', 1)
legend('HARV', 'Virginia', 'Ice Cream Cone', 'Bunny Bomb')
title('C_D vs. Angle of Attack at Mach 5')
xlabel('Alpha (deg)')
ylabel('C_D')
hold off

%plot CL vs. alpha at Mach 5
figure()
hold on
plot(linspace(-30, 30, 1000), HARV_CL(501,:), 'LineWidth', 1)
plot(linspace(-30, 30, 1000), Virginia_CL(501,:), 'LineWidth', 1)
plot(linspace(-30, 30, 1000), Ice_Cream_Cone_CL(501,:), 'LineWidth', 1)
plot(linspace(-30, 30, 100), Bunny_Bomb_CL(51,:), 'LineWidth', 1)
legend('HARV', 'Virginia', 'Ice Cream Cone', 'Bunny Bomb')
title('C_L vs. Angle of Attack at Mach 5')
xlabel('Alpha (deg)')
ylabel('C_L')
hold off

%plot L/D vs. alpha at Mach 5
figure()
hold on
plot(linspace(-30, 30, 1000), HARV_L_D(501,:), 'LineWidth', 1)
plot(linspace(-30, 30, 1000), Virginia_L_D(501,:), 'LineWidth', 1)
plot(linspace(-30, 30, 1000), Ice_Cream_Cone_L_D(501,:), 'LineWidth', 1)
plot(linspace(-30, 30, 100), Bunny_Bomb_L_D(51,:), 'LineWidth', 1)
legend('HARV', 'Virginia', 'Ice Cream Cone', 'Bunny Bomb')
title('L/D vs. Angle of Attack at Mach 5')
xlabel('Alpha (deg)')
ylabel('L/D')
hold off