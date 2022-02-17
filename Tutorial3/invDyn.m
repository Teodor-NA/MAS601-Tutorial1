% pos = {x2; y2; theta2; x3; y3; theta3; x4; y4; theta4]
% acc = {ddx2; ddy2; ddtheta2; ddx3; ddy3; ddtheta3; ddx4; ddy4; ddtheta4]
% L = [L2; L3; L4]
function  lambdaT = invDyn(pos, acc, M, h, L)
    theta2 = pos(3);
    theta3 = pos(6);
    theta4 = pos(9);
    L2 = L(1);
    L3 = L(2);
    L4 = L(3);

    D = [
        [eye(2), L2/2*[sin(theta2); -cos(theta2)], zeros(2, 6)];
        [-eye(2), L2/2*[sin(theta2); -cos(theta2)], eye(2), L3/2*[sin(theta3); -cos(theta3)], zeros(2, 3)]; 
        [zeros(2, 3), -eye(2), L3/2*[sin(theta3); -cos(theta3)], eye(2), L4/2*[sin(theta4); -cos(theta4)]];
        [zeros(2, 6), eye(2), L4/2*[sin(theta4-pi); -cos(theta4-pi)]];
        [zeros(1, 2), 1, zeros(1, 6)]
    ];
    
    DTrans = D';
    
    lambdaT = DTrans\(M*acc - h);
end