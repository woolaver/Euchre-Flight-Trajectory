%Using reference values from https://arc.aiaa.org/doi/pdfplus/10.2514/6.2025-1338
%compare HARV CL and CD values from blended driver to determine accuracy

clear
clc
close all

%Reference data digitized from plots 
alpha = deg2rad(linspace(0, 28, 15));

% (a) Mach 4, 20 km
ref_M4_20km_CL = [0 0.22 0.48 0.74 0.95 1.18 1.43 1.68 1.94 2.24 2.55 2.84 3.14 3.42 3.72];
ref_M4_20km_CD = [0.115 0.13 0.16 0.21 0.27 0.35 0.46 0.59 0.73 0.91 1.13 1.36 1.63 1.90 2.20];

% (b) Mach 12, 20 km
ref_M12_20km_CL = [0 0.12 0.25 0.40 0.58 0.75 0.96 1.20 1.46 1.72 2.04 2.32 2.66 3.00 3.32];
ref_M12_20km_CD = [0.085 0.10 0.115 0.14 0.19 0.25 0.34 0.45 0.59 0.74 0.94 1.15 1.40 1.70 2.00];

% (c) Mach 4, 40 km
ref_M4_40km_CL = [0 0.22 0.47 0.67 0.90 1.15 1.40 1.68 1.97 2.26 2.56 2.84 3.16 3.43 3.73];
ref_M4_40km_CD = [0.125 0.13 0.165 0.20 0.25 0.34 0.46 0.59 0.75 0.93 1.13 1.36 1.64 1.93 2.23];

% (d) Mach 12, 40 km
ref_M12_40km_CL = [0 0.12 0.25 0.40 0.57 0.75 0.96 1.20 1.45 1.71 2.03 2.32 2.66 2.98 3.32];
ref_M12_40km_CD = [0.085 0.10 0.115 0.14 0.185 0.245 0.34 0.445 0.59 0.745 0.945 1.15 1.40 1.72 2.03];

% (a) Mach 4, 20 km
ref_M4_20km_LD = [0 1.74 2.93 3.45 3.53 3.38 3.13 2.88 2.66 2.45 2.25 2.08 1.93 1.80 1.67];

% (b) Mach 12, 20 km
ref_M12_20km_LD = [0 1.20 2.20 2.82 3.08 3.10 2.83 2.68 2.50 2.33 2.17 2.02 1.88 1.76 1.65];

% (c) Mach 4, 40 km
ref_M4_40km_LD = [0 1.73 2.88 3.38 3.50 3.35 3.10 2.88 2.66 2.46 2.25 2.08 1.92 1.77 1.65];

% (d) Mach 12, 40 km
ref_M12_40km_LD = [0 1.25 2.22 2.82 3.10 3.12 2.85 2.70 2.50 2.32 2.18 2.03 1.88 1.74 1.63];

%% Blended Driver Data

geometry = "csv_files/HARV.STL";
Sref = .359; %m^2
gamma = 1.4;
Tw = 500; %not needed cause just need skin friction but just using to run refTempMethod
length = 1.778;
height = .089;

blended_M4_20km_CL = zeros(size(alpha));
blended_M4_40km_CL = zeros(size(alpha));
blended_M12_20km_CL = zeros(size(alpha));
blended_M12_40km_CL = zeros(size(alpha));

blended_M4_20km_CD = zeros(size(alpha));
blended_M4_40km_CD = zeros(size(alpha));
blended_M12_20km_CD = zeros(size(alpha));
blended_M12_40km_CD = zeros(size(alpha));

for i = 1:numel(alpha)
    M_inf = 4;
    [blended_M4_20km_CL(i), blended_M4_20km_CD(i)] = blended_driver(geometry, alpha(i), M_inf, Sref, gamma);
    blended_M4_40km_CL(i) = blended_M4_20km_CL(i);
    blended_M4_40km_CD(i) = blended_M4_20km_CD(i);

    alt = 20000;
    [~, a, ~, rho] = atmoscoesa(alt);
    V = M_inf*a;
    [Fx, Fy, ~] = refTempMethod(length, height, alt, M_inf, alpha(i), Tw);

    if(isnan(Fx))
        [Fx, Fy, ~] = refTempMethod(length, height, alt, M_inf, deg2rad(3), Tw);
    end
    
    Lf = Fy*cos(alpha(i)) - Fx*sin(alpha(i));
    Df = Fy*sin(alpha(i)) + Fx*cos(alpha(i));

    CL_f = 2*Lf/(rho*V^2*Sref);
    CD_f = 2*Df/(rho*V^2*Sref);

    blended_M4_20km_CL(i) = blended_M4_20km_CL(i) + CL_f;
    blended_M4_20km_CD(i) = blended_M4_20km_CD(i) + CD_f;

    alt = 40000;
    [~, a, ~, rho] = atmoscoesa(alt);
    V = M_inf*a;
    [Fx, Fy, ~] = refTempMethod(length, height, alt, M_inf, alpha(i), Tw);

    if(isnan(Fx))
        [Fx, Fy, ~] = refTempMethod(length, height, alt, M_inf, deg2rad(3), Tw);
    end
    
    Lf = Fy*cos(alpha(i)) - Fx*sin(alpha(i));
    Df = Fy*sin(alpha(i)) + Fx*cos(alpha(i));

    CL_f = 2*Lf/(rho*V^2*Sref);
    CD_f = 2*Df/(rho*V^2*Sref);

    blended_M4_40km_CL(i) = blended_M4_40km_CL(i) + CL_f;
    blended_M4_40km_CD(i) = blended_M4_40km_CD(i) + CD_f;

    %code does not take altitude into account so same for 20 and 40 km
    M_inf = 12;
    [blended_M12_20km_CL(i), blended_M12_20km_CD(i)] = blended_driver(geometry, alpha(i), M_inf, Sref, gamma);
    blended_M12_40km_CL(i) = blended_M12_20km_CL(i);
    blended_M12_40km_CD(i) = blended_M12_20km_CD(i);
    
    alt = 20000;
    M_inf = 10; %need to use 10 because cone tables max at 10
    [~, a, ~, rho] = atmoscoesa(alt);
    V = M_inf*a;
    [Fx, Fy, ~] = refTempMethod(length, height, alt, M_inf, alpha(i), Tw);

    if(isnan(Fx))
        [Fx, Fy, ~] = refTempMethod(length, height, alt, M_inf, deg2rad(3), Tw);
    end
    
    Lf = Fy*cos(alpha(i)) - Fx*sin(alpha(i));
    Df = Fy*sin(alpha(i)) + Fx*cos(alpha(i));

    CL_f = 2*Lf/(rho*V^2*Sref);
    CD_f = 2*Df/(rho*V^2*Sref);

    blended_M12_20km_CL(i) = blended_M12_20km_CL(i) + CL_f;
    blended_M12_20km_CD(i) = blended_M12_20km_CD(i) + CD_f;   

    alt = 40000;
    [~, a, ~, rho] = atmoscoesa(alt);
    V = M_inf*a;
    [Fx, Fy, ~] = refTempMethod(length, height, alt, M_inf, alpha(i), Tw);

    if(isnan(Fx))
        [Fx, Fy, ~] = refTempMethod(length, height, alt, M_inf, deg2rad(3), Tw);
    end
    
    Lf = Fy*cos(alpha(i)) - Fx*sin(alpha(i));
    Df = Fy*sin(alpha(i)) + Fx*cos(alpha(i));

    CL_f = 2*Lf/(rho*V^2*Sref);
    CD_f = 2*Df/(rho*V^2*Sref);

    blended_M12_40km_CL(i) = blended_M12_40km_CL(i) + CL_f;
    blended_M12_40km_CD(i) = blended_M12_40km_CD(i) + CD_f;
end

blended_M4_20km_LD = blended_M4_20km_CL./blended_M4_20km_CD;
blended_M4_40km_LD = blended_M4_40km_CL./blended_M4_40km_CD;
blended_M12_20km_LD = blended_M12_20km_CL./blended_M12_20km_CD;
blended_M12_40km_LD = blended_M12_40km_CL./blended_M12_40km_CD;

percent_error_M4_20km = (abs(blended_M4_20km_LD - ref_M4_20km_LD)./ref_M4_20km_LD).*100;
percent_error_M4_40km = (abs(blended_M4_40km_LD - ref_M4_40km_LD)./ref_M4_40km_LD).*100;
percent_error_M12_20km = (abs(blended_M12_20km_LD - ref_M12_20km_LD)./ref_M12_20km_LD).*100;
percent_error_M12_40km = (abs(blended_M12_40km_LD - ref_M12_40km_LD)./ref_M12_40km_LD).*100;


figure(); 
subplot(2,2,1); 
hold on
plot(ref_M4_20km_CD, ref_M4_20km_CL, '-o', 'LineWidth', 1)
plot(blended_M4_20km_CD, blended_M4_20km_CL, '-o', 'LineWidth', 1)
xlabel('C_D') 
ylabel('C_L')
legend("Reference Values", "Blended Driver Values")
title('Mach 4 20km Drag Polar')
hold off

subplot(2,2,2); 
hold on
plot(ref_M12_20km_CD, ref_M12_20km_CL, '-o', 'LineWidth', 1)
plot(blended_M12_20km_CD, blended_M12_20km_CL, '-o', 'LineWidth', 1)
xlabel('C_D') 
ylabel('C_L')
legend("Reference Values", "Blended Driver Values")
title('Mach 12 20km Drag Polar')
hold off

subplot(2,2,3); 
hold on
plot(ref_M4_40km_CD, ref_M4_40km_CL, '-o', 'LineWidth', 1)
plot(blended_M4_40km_CD, blended_M4_40km_CL, '-o', 'LineWidth', 1)
xlabel('C_D') 
ylabel('C_L')
legend("Reference Values", "Blended Driver Values")
title('Mach 4 40km Drag Polar')
hold off

subplot(2,2,4); 
hold on
plot(ref_M12_40km_CD, ref_M12_40km_CL, '-o', 'LineWidth', 1)
plot(blended_M12_40km_CD, blended_M12_40km_CL, '-o', 'LineWidth', 1)
xlabel('C_D') 
ylabel('C_L')
legend("Reference Values", "Blended Driver Values")
title('Mach 12 40km Drag Polar')
hold off


figure(); 
subplot(2,2,1); 
hold on
plot(rad2deg(alpha), ref_M4_20km_LD, '-o', 'LineWidth', 1)
plot(rad2deg(alpha), blended_M4_20km_LD, '-o', 'LineWidth', 1)
xlabel('\alpha (deg)')
ylabel('L/D')
legend("Reference Values", "Blended Driver Values")
title('Mach 4 20km L/D vs. \alpha')
hold off

subplot(2,2,2); 
hold on
plot(rad2deg(alpha), ref_M12_20km_LD, '-o', 'LineWidth', 1)
plot(rad2deg(alpha), blended_M12_20km_LD, '-o', 'LineWidth', 1)
xlabel('\alpha (deg)')
ylabel('L/D')
legend("Reference Values", "Blended Driver Values")
title('Mach 12 20km L/D vs. \alpha')
hold off

subplot(2,2,3); 
hold on
plot(rad2deg(alpha), ref_M4_40km_LD, '-o', 'LineWidth', 1)
plot(rad2deg(alpha), blended_M4_40km_LD, '-o', 'LineWidth', 1)
xlabel('\alpha (deg)')
ylabel('L/D')
legend("Reference Values", "Blended Driver Values")
title('Mach 4 40km L/D vs. \alpha')
hold off

subplot(2,2,4); 
hold on
plot(rad2deg(alpha), ref_M12_40km_LD, '-o', 'LineWidth', 1)
plot(rad2deg(alpha), blended_M12_40km_LD, '-o', 'LineWidth', 1)
xlabel('\alpha (deg)')
ylabel('L/D')
legend("Reference Values", "Blended Driver Values")
title('Mach 12 40km L/D vs. \alpha')
hold off

figure(); 
subplot(2,2,1); 
plot(rad2deg(alpha), percent_error_M4_20km, '-o', 'LineWidth', 1)
xlabel('\alpha (deg)')
ylabel('Percent Error L/D')
title('Percent Error Mach 4 20km of L/D')

subplot(2,2,2); 
plot(rad2deg(alpha), percent_error_M12_20km, '-o', 'LineWidth', 1)
xlabel('\alpha (deg)')
ylabel('Percent Error L/D')
title('Percent Error Mach 12 20km of L/D')

subplot(2,2,3); 
plot(rad2deg(alpha), percent_error_M4_40km, '-o', 'LineWidth', 1)
xlabel('\alpha (deg)')
ylabel('Percent Error L/D')
title('Percent Error Mach 4 40km of L/D')

subplot(2,2,4); 
plot(rad2deg(alpha), percent_error_M12_40km, '-o', 'LineWidth', 1)
xlabel('\alpha (deg)')
ylabel('Percent Error L/D')
title('Percent Error Mach 12 40km of L/D')