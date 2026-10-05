function gamma = gamma_empirical(T)
%Use the emperical formula for cp to calculate gamma
%   Detailed explanation goes here
    if (T > 1000)
        a1 = 2.42144865*10^4;
        a2 = -1.18915152*10^2;
        a3 = 4.00427357*10^0;
        a4 = 3.31011382*10^(-4);
        a5 = -1.40223124*10^(-7);
        a6 = 2.73033502*10^(-11);
        a7 = -2.03055745*10^(-15);
    else
        a1 = 4.18431872*10^3;
        a2 = -4.30514102*10^1;
        a3 = 3.56839620*10^0;
        a4 = -6.78872905*10^(-4);
        a5 = 1.55370228*10^(-6);
        a6 = -3.29937060*10^(-12);
        a7 = -4.66395342*10^(-15);
    end
    cp_R = a1*T^(-2) + a2*T^(-1) + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4;

    gamma = (cp_R)/(cp_R - 1);
end