function a = tangent_wedge(minf,h,theta)
    y=1.4;
    R=287;
    [rho,p,T] = stdatmfull(h);
    [p2op1, t2op1, rho2op1, ~, M2] = tbm(minf, theta);
    if
    a.p = p*p2op1;
    a.rho = rho*rho2op1;
    a.T = T*t2op1;
    a.M = M2;
    a.ue = M2/sqrt(y*R*a.T);



end