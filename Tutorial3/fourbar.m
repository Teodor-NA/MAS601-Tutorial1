function Phi = fourbar(x, theta2)
% Array (x contains theta3 and theta4)
global L1 L2 L3 L4 Theta1

Phi = [L2*cos(theta2) + L3*cos(x(1)) + L4*cos(x(2)) - L1*cos(Theta1)
       L2*sin(theta2) + L3*sin(x(1)) + L4*sin(x(2)) - L1*sin(Theta1)];
end
