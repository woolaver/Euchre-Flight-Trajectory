function [a,ok] = tangent_wedge(minf,h,theta)
    ok=true;
    y=1.4;
    R=287;
    [rho,p,T] = stdatmfull(h);
    try
        [p2op1, t2op1, rho2op1, ~, M2] = tbm(minf, theta);
    catch
        ok=false;
    end
    a.p = p*p2op1;
    a.rho = rho*rho2op1;
    a.T = T*t2op1;
    a.M = M2;
    a.ue = M2/sqrt(y*R*a.T);



end