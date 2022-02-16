function qDot = odeFunct(t, q)

global L1 L2 L3 L4 Theta1 MInv T tSpan dt m2 m3 m4 g
    
    s2 = q(1:2);
    theta2 = q(3);
    s3 = q(4:5);
    theta3 = q(6);
    s4 = q(7:8);
    theta4 = q(9);
    sDot2 = q(10:11);
    thetaDot2 = q(12);
    sDot3 = q(13:14);
    thetaDot3 = q(15);
    sDot4 = q(16:17);
    thetaDot4 = q(18);

    D = zeros(8,9);
    D(1:2, 1:3) = [eye(2), L2/2*[sin(theta2); -cos(theta2)]];
    D(3:4, 1:3) = [-eye(2), L2/2*[sin(theta2); -cos(theta2)]];
    D(3:4, 4:6) = [eye(2), L3/2*[sin(theta3); -cos(theta3)]];
    D(5:6, 4:6) = [-eye(2), L3/2*[sin(theta3); -cos(theta3)]];
    D(5:6, 7:9) = [eye(2), L4/2*[sin(theta4); -cos(theta4)]];
    D(7:8, 7:9) = [eye(2), L4/2*[sin(theta4-pi); -cos(theta4-pi)]];
    DTrans = D';
        
    gamma = [
        -L2*(thetaDot2^2)/2*[cos(theta2); sin(theta2)];
        -L2*(thetaDot2^2)/2*[cos(theta2); sin(theta2)] - L3*(thetaDot3^2)/2*[cos(theta3); sin(theta3)];
        -L3*(thetaDot3^2)/2*[cos(theta3); sin(theta3)] - L4*(thetaDot4^2)/2*[cos(theta4); sin(theta4)];
        -L4*(thetaDot4^2)/2*[cos(theta4-pi); sin(theta4-pi)];        
    ];
    
    index = find(abs(tSpan - t) < dt/2);
    
    if length(index) ~= 1
        index
    end
    
    ha = [0; -m2*g; T(index); 0; -m3*g; 0; 0; -m4*g; 0;];
    
    cDotDot = MInv*(DTrans)*inv(D*MInv*(DTrans))*(gamma - D*MInv*ha) + MInv*ha;
    
    qDot = [q(10:18); cDotDot];

end