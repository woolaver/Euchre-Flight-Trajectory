CD_vec = readmatrix("HARV_CD_50_500.csv");
CL_vec = readmatrix("HARV_CL_50_500.csv");

alpha_vec = deg2rad(linspace(-30, 30, 500));
mach_vec = linspace(2, 8, 50);

alpha_test = deg2rad(linspace(-30, 30, 5000));
mach_test = linspace(2, 8, 500);

for i = 1:5000
    for j = 1:500
        CL = interp2(alpha_vec, mach_vec, CL_vec, alpha_test(i), mach_test(j));
        CD = interp2(alpha_vec, mach_vec, CD_vec, alpha_test(i), mach_test(j));

        if(isnan(CL))
            disp("Fuck justin")
        end
        if(isnan(CD))
            disp("FUCK JUSTIN")
        end
    end
end