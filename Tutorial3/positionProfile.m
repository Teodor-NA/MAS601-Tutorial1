function q2ref = positionProfile(t)
    q2ref = zeros(3,1);
    q2ref(1) = 0.1*sin(pi/2*t);
    q2ref(2) = 0.1*pi/2*cos(pi/2*t);
    q2ref(3) = -0.1*(pi^2)/4*sin(pi/2*t);
    
end