function qDot = ODE_Funct(t, q)
    global m1 m2 m3 k1 k2 k3 b1 b2 b3 a d c L0 t1 t2 t3 F1 F2;

    x1 = q(1);
    x2 = q(2);
    x3 = q(3);
    xDot1 = q(4);
    xDot2 = q(5);
    xDot3 = q(6);

    Fk1 = (x1 - L0)*k1;
    Fk2 = (x2 - x1 - d - L0)*k2;
    Fk3 = (x3 - L0)*k3;

    if (x3 - x1 > d) || (x1 > x3 + 8*d)
        Fb1 = 0;
    else
        Fb1 = (xDot1 - xDot3)*b1;
    end
    Fb2 = xDot1*b2;
    if (x3 - x2 > d) || (x2 > x3 + 8*d)
        Fb3 = 0;
    else
        Fb3 = (xDot2 - xDot3)*b3;
    end
    
%    Fb1 = 0;
%    Fb2 = 0;
%    Fb3 = 0;

    if t < t1
        f = F1*t/t1;
    elseif t < t2
        f = F1 + (F2 - F1)/(t2 - t1)*(t - t1);
    elseif t < t3
        f = F2 + (0 - F2)/(t3 - t2)*(t - t2);
    else
        f = 0;
    end
    
%      f = 0;

%     f = 10;
    
    xDotDot1 = (-Fb2 - Fb1 - Fk1 + Fk2)/m1;
    xDotDot2 = (-Fk2 - Fb3 + f)/m2;
    xDotDot3 = (-Fk3 + Fb1 + Fb3)/m3;
    
    % qDot = zeros(6, 1);
    qDot = [xDot1; xDot2; xDot3; xDotDot1; xDotDot2; xDotDot3];
end