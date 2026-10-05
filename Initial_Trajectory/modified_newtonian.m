function [CL, CD] = modified_newtonian(alpha, total_vel, alt, T_prev)
%Simple Modified Newtonian Theory for caluculation of CL and CD
    [~, a, P1, ~] = atmoscoesa(alt);
    gamma = gamma_empirical(T_prev);
    
    Mach = total_vel/a;
    

end