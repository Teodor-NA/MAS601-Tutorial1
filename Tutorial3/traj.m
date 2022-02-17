function q2ref = traj(t, amplitude)
    q2ref = zeros(3,1);
    q2ref(1) = amplitude*sin(pi/2*t);
    q2ref(2) = amplitude*pi/2*cos(pi/2*t);
    q2ref(3) = -amplitude*(pi^2)/4*sin(pi/2*t);   
end