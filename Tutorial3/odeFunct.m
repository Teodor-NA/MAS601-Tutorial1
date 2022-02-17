function qDot = odeFunct(t, q)

global L2 L3 L4 M MInv Theta2Max h
    % Extract what we need from q
    theta2 = q(3);
    theta3 = q(6);
    theta4 = q(9);
    dtheta2 = q(12);
    dtheta3 = q(15);
    dtheta4 = q(18);

    % Pre-compute sin and cos values
    c2 = cos(theta2);
    s2 = sin(theta2);
    c3 = cos(theta3);
    s3 = sin(theta3);
    c4 = cos(theta4);
    s4 = sin(theta4);
    c4_pi = cos(theta4 - pi);
    s4_pi = sin(theta4 - pi);

    % Find trajectory angle, angle vel and angle acc for current timestep
    qRef = traj(t, Theta2Max);
    
    % Find appropriate Torque for current timestep
    [pos, ~, acc] = invKin(qRef, [theta3; theta4], [L2; L3; L4]);
    lambdaT = invDyn(pos, acc, M, h, [L2; L3; L4]);
    T = lambdaT(9);
    
    % Calculate forward dynamics
    D = zeros(8,9);
    D(1:2, 1:3) = [eye(2), L2/2*[s2; -c2]];
    D(3:4, 1:3) = [-eye(2), L2/2*[s2; -c2]];
    D(3:4, 4:6) = [eye(2), L3/2*[s3; -c3]];
    D(5:6, 4:6) = [-eye(2), L3/2*[s3; -c3]];
    D(5:6, 7:9) = [eye(2), L4/2*[s4; -c4]];
    D(7:8, 7:9) = [eye(2), L4/2*[s4_pi; -c4_pi]];
    DTrans = D';
        
    gamma = [
        -L2*(dtheta2^2)/2*[c2; s2];
        -L2*(dtheta2^2)/2*[c2; s2] - L3*(dtheta3^2)/2*[c3; s3];
        -L3*(dtheta3^2)/2*[c3; s3] - L4*(dtheta4^2)/2*[c4; s4];
        -L4*(dtheta4^2)/2*[c4_pi; s4_pi];        
    ];
    
    ha = h + [zeros(2,1); T; zeros(6,1)];
    
    tmpA = MInv*DTrans;
    tmpB = D*tmpA;
    tmpC = MInv*ha;
    
    cDotDot = tmpA/tmpB*(gamma - D*tmpC) + tmpC;
 
    qDot = [q(10:18); cDotDot];

end