% Inputs
% qRef = [theta2; d/dt(theta2); d^2/dt^2(theta2)] (reference trajectory)
% theta34 = [theta3; theta4]
% L = [L2; L3; L4]
% options = fsolve options
% Outputs
% pos = [x2; y2; theta2; x3; y3; theta3; x4; y4; theta4]
% vel = d/dt(pos)
% acc = d/dt(vel)
function [pos, vel, acc] = invKin(qRef, theta34, L)
    % Angular kinematics
    % Position
    theta2 = qRef(1);
    dtheta2 = qRef(2);
    ddtheta2 = qRef(3);

    L2 = L(1);
    L3 = L(2);
    L4 = L(3);

    theta3 = theta34(1);
    theta4 = theta34(2);
    
    % Pre-compute sin and cos values
    c2 = cos(theta2);
    s2 = sin(theta2);
    c3 = cos(theta3);
    s3 = sin(theta3);
    c4 = cos(theta4);
    s4 = sin(theta4);
    
    % Velocity    
    A = [
        -L3*s3, -L4*s4;
        L3*c3, L4*c4
    ];
    B = [
        L2*s2*dtheta2; 
        -L2*c2*dtheta2
    ];
    dtheta34 = A\B;

    dtheta3 = dtheta34(1);
    dtheta4 = dtheta34(2);

    % Acceleration
    C = [
        L2*(s2*ddtheta2 + c2*dtheta2^2) + L3*c3*dtheta3^2 + L4*c4*dtheta4^2;
        -L2*(c2*ddtheta2 - s2*dtheta2^2) + L3*s3*dtheta3^2 + L4*s4*dtheta4^2
    ];

    ddtheta34 = A\C;

    ddtheta3 = ddtheta34(1);
    ddtheta4 = ddtheta34(2);

    % Linear Kinematics    
    Rg2 = L2/2*[c2; s2];
    Rg3 = Rg2 + L2/2*[c2; s2] + L3/2*[c3; s3];
    Rg4 = Rg3 + L3/2*[c3; s3] + L4/2*[c4; s4];

    Vg2 = L2/2*dtheta2*[-s2; c2];
    Vg3 = Vg2 + L2/2*dtheta2*[-s2; c2] + L3/2*dtheta3*[-s3; c3];
    Vg4 = Vg3 + L3/2*dtheta3*[-s3; c3] + L4/2*dtheta4*[-s4; c4];

    Ag2 = L2/2*ddtheta2*[-s2; c2] - L2/2*dtheta2^2*[c2; s2];
    Ag3 = Ag2 + L2/2*ddtheta2*[-s2; c2] - L2/2*dtheta2^2*[c2; s2] + L3/2*ddtheta3*[-s3; c3] - L3/2*dtheta3^2*[c3; s3];
    Ag4 = Ag3 + L3/2*ddtheta3*[-s3; c3] - L3/2*dtheta3^2*[c3; s3] + L4/2*ddtheta4*[-s4; c4] - L4/2*dtheta4^2*[c4; s4];

    % Return values
    pos = [Rg2; theta2; Rg3; theta3; Rg4; theta4];
    vel = [Vg2; dtheta2; Vg3; dtheta3; Vg4; dtheta4];
    acc = [Ag2; ddtheta2; Ag3; ddtheta3; Ag4; ddtheta4];    
end