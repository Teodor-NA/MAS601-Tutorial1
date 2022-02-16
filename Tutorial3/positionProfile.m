function q2ref = positionProfile(t)
    a = 0.1;
    q2ref = zeros(3,1);
    q2ref(1) = a*sin(pi/2*t);
    q2ref(2) = a*pi/2*cos(pi/2*t);
    q2ref(3) = -a*(pi^2)/4*sin(pi/2*t);
    
end