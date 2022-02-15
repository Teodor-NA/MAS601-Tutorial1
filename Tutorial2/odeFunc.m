function qDot = odeFunc(t, q)
    global L1 L2 Fg1 Fg2 L0 k b Minv
    s1 = q(1:2);
    phi1 = q(3);
    s2 = q(4:5);
    phi2 = q(6);
    sDot1 = q(7:8);
    phiDot1 = q(9);
    sDot2 = q(10:11);
    phiDot2 = q(12);

    aux1 = L1/2*(phiDot1^2);
    aux2 = L2/2*(phiDot2^2);

    gamma = [
        -aux1*cos(phi1);
        -aux1*sin(phi1);
        aux1*cos(phi1)+aux2*cos(phi2);
        aux1*sin(phi1)+aux2*sin(phi2)
    ];

    d = s2 - s1;
    L = sqrt(d.'*d);
    u = d/L;

    Fk = k*(L - L0)*u;
    
    dDot = sDot2 - sDot1;
    LDot = u.'*dDot;
    Fb = b*LDot*u;

    ha = [
        Fk + Fb + Fg1; 
        0; 
        -Fk - Fb + Fg2; 
        0
    ];

    D = [
        1, 0, L1/2*sin(phi1), zeros(1,3); 
        0, 1, -L1/2*cos(phi1), zeros(1,3);
        1, 0, -L1/2*sin(phi1), -1, 0, -L2/2*sin(phi2);
        0, 1, L1/2*cos(phi1), 0, -1, L2/2*cos(phi2)
    ];

    cDotDot = Minv*(D.')*inv(D*Minv*(D.'))*(gamma - D*Minv*ha) + Minv*ha;
    
    qDot = [q(7:12); cDotDot];

end