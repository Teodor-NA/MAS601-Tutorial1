clear; clc;

global L1 L2 L3 L4 Theta1
Theta1 = 0; 
L1 = 5;  % m 
 
m2 = 2; % kg 
L2 = 0.5; %m 
Jbar2 = 1/12*m2*L2^2; %kg.m^2 
 
m3 = 1; % kg 
L3 = 8; % m 
Jbar3 = 1/12*m3*L3^2; %kg.m^2 
 
m4 = 3; % kg 
L4 = 4; % m 
Jbar4 = 1/12*m4*L2^4; %kg.m^2 
theta2 = pi/4;

dt = 0.01;
tEnd = 10;
tSpan = 0:dt:tEnd;
N = length(tSpan);

guesses = [15; 300]*pi/180;% estimates

qRef2 = zeros(3,N);
theta34 = zeros(2, N+1);
theta34(:, 1) = guesses;
thetaDot34 = zeros(2, N);
options = optimset('display', 'off');
for i = 1:N
    qRef2(:,i) = positionProfile(tSpan(i));
    theta34(:,i+1) = fsolve(@fourbar, theta34(:,i), options, qRef2(1,i));

    A = [-L3*sin(theta34(1,i+1)), -L4*sin(theta34(2,i));
        L3*cos(theta34(1,i+1)), L4*cos(theta34(2,i))];
    B = [L2*sin( qRef2(1,i)*qRef2(2,i)); -L2*cos( qRef2(1,i)*qRef2(2,i))];
    thetaDot34(:,i) = A\B;
end

figure;
plot(tSpan, 180/pi*theta34(:, 2:N+1));
hold on
plot(tSpan, 180/pi*qRef2(1,:));
legend('\theta_3','\theta_4', '\theta_2')

% qInit2 = positionProfile(0);
 

% options = optimset('display', 'off');
% angles = fsolve(@fourbar, theta34, options, theta2)
% theta3 = angles(1)
% theta4 = angles(2)

