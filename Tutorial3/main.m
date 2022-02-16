clear; clc;

global L1 L2 L3 L4 Theta1 MInv T tSpan dt m2 m3 m4 g

animate = true;
animateSim = true;

%% Init
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

g = 9.81;
M = diag([m2, m2, Jbar2, m3, m3, Jbar3, m4, m4, Jbar4]);
MInv = M^-1;
ha = [0; -m2*g; 0; 0; -m3*g; 0; 0; -m4*g; 0;];

dt = 0.01;
tEnd = 20;
tSpan = 0:dt:tEnd;
N = length(tSpan);

theta34Estimate = [15; 300]*pi/180;% estimates

%% Solve
qRef2 = zeros(3,N);
theta34 = zeros(2, N);
thetaDot34 = zeros(2, N);
thetaDotDot34 = zeros(2, N);
Rg2 = zeros(2, N);
Rg3 = zeros(2, N);
Rg4 = zeros(2, N);
Vg2 = zeros(2, N);
Vg3 = zeros(2, N);
Vg4 = zeros(2, N);
Ag2 = zeros(2, N);
Ag3 = zeros(2, N);
Ag4 = zeros(2, N);

lambdaT = zeros(9, N);

options = optimset('display', 'off');
for i = 1:N
    % Position
    qRef2(:,i) = positionProfile(tSpan(i));
    theta2 = qRef2(1,i);
    thetaDot2 = qRef2(2,i);
    thetaDotDot2 = qRef2(3,i);

    if (i == 1)
        theta34(:,i) = fsolve(@fourbar, theta34Estimate, options, theta2);
    else
        theta34(:,i) = fsolve(@fourbar, theta34(:,i - 1), options, theta2);
    end

    theta3 = theta34(1,i);
    theta4 = theta34(2,i);

    % Velocity
    A = [
        -L3*sin(theta3), -L4*sin(theta4);
        L3*cos(theta3), L4*cos(theta4)
    ];
    B = [
        L2*sin(theta2)*thetaDot2; 
        -L2*cos(theta2)*thetaDot2
    ];
    thetaDot34(:,i) = A\B;

    thetaDot3 = thetaDot34(1,i);
    thetaDot4 = thetaDot34(2,i);

    % Acceleration
    C = [
        L2*(sin(theta2)*thetaDotDot2 + cos(theta2)*thetaDot2^2) + L3*cos(theta3)*thetaDot3^2 + L4*cos(theta4)*thetaDot4^2;
        -L2*(cos(theta2)*thetaDotDot2 - sin(theta2)*thetaDot2^2) + L3*sin(theta3)*thetaDot3^2 + L4*sin(theta4)*thetaDot4^2
    ];

    thetaDotDot34(:,i) = A\C;

    thetaDotDot3 = thetaDotDot34(1,i);
    thetaDotDot4 = thetaDotDot34(2,i);

    Rg2(:,i) = L2/2*[cos(theta2); sin(theta2)];
    Rg3(:,i) = Rg2(:,i) + L2/2*[cos(theta2); sin(theta2)] + L3/2*[cos(theta3); sin(theta3)];
    Rg4(:,i) = Rg3(:,i) + L3/2*[cos(theta3); sin(theta3)] + L4/2*[cos(theta4); sin(theta4)];

    Vg2(:,i) = L2/2*thetaDot2*[-sin(theta2); cos(theta2)];
    Vg3(:,i) = Vg2(:,i) + L2/2*thetaDot2*[-sin(theta2); cos(theta2)] + L3/2*thetaDot3*[-sin(theta3); cos(theta3)];
    Vg4(:,i) = Vg3(:,i) + L3/2*thetaDot3*[-sin(theta3); cos(theta3)] + L4/2*thetaDot4*[-sin(theta4); cos(theta4)];

    Ag2(:,i) = L2/2*thetaDotDot2*[-sin(theta2); cos(theta2)] - L2/2*thetaDot2^2*[cos(theta2); sin(theta2)];
    Ag3(:,i) = Ag2(:,i) + L2/2*thetaDotDot2*[-sin(theta2); cos(theta2)] - L2/2*thetaDot2^2*[cos(theta2); sin(theta2)] + L3/2*thetaDotDot3*[-sin(theta3); cos(theta3)] - L3/2*thetaDot3^2*[cos(theta3); sin(theta3)];
    Ag4(:,i) = Ag3(:,i) + L3/2*thetaDotDot3*[-sin(theta3); cos(theta3)] - L3/2*thetaDot3^2*[cos(theta3); sin(theta3)] + L4/2*thetaDotDot4*[-sin(theta4); cos(theta4)] - L4/2*thetaDot4^2*[cos(theta4); sin(theta4)];

    D = zeros(9);
    D(1:2, 1:3) = [eye(2), L2/2*[sin(theta2); -cos(theta2)]];
    D(3:4, 1:3) = [-eye(2), L2/2*[sin(theta2); -cos(theta2)]];
    D(3:4, 4:6) = [eye(2), L3/2*[sin(theta3); -cos(theta3)]];
    D(5:6, 4:6) = [-eye(2), L3/2*[sin(theta3); -cos(theta3)]];
    D(5:6, 7:9) = [eye(2), L4/2*[sin(theta4); -cos(theta4)]];
    D(7:8, 7:9) = [eye(2), L4/2*[sin(theta4-pi); -cos(theta4-pi)]];
    D(9,3) = 1;
    DTrans = D';
    
    cDotDot = [
        Ag2(:,i);
        thetaDotDot2;
        Ag3(:,i);
        thetaDotDot3;
        Ag4(:,i);
        thetaDotDot4
    ];
    
    lambdaT(:, i) = DTrans\(M*cDotDot - ha);
%     gamma = [
%         -L2*(thetaDot2^2)/2*[cos(theta2); sin(theta2)];
%         -L2*(thetaDot2^2)/2*[cos(theta2); sin(theta2)] - L3*(thetaDot3^2)/2*[cos(theta3); sin(theta3)];
%         -L3*(thetaDot3^2)/2*[cos(theta3); sin(theta3)] - L4*(thetaDot4^2)/2*[cos(theta4); sin(theta4)];
%         -L4*(thetaDot4^2)/2*[cos(theta4-pi); sin(theta4-pi)];        
%     ];


end


%% Invoking ode45
q = [
    Rg2(:, 1);
    qRef2(1, 1);
    Rg3(:, 1);
    theta34(1, 1);
    Rg4(:, 1);
    theta34(2, 1);
    Vg2(:, 1);
    qRef2(2, 1);
    Vg3(:, 1);
    thetaDot34(1, 1);
    Vg4(:, 1);
    thetaDot34(2, 1);
];

T = lambdaT(9, :);
options = odeset('RelTol', 1e-15, 'AbsTol', 1e-15);
[Time, qT] = ode45(@odeFunct, tSpan, q);

%% Plot
close all;

% Driven
plotRows = 2;
plotCols = 3;
plotIndex = 0;

fig = figure;
plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, 180/pi*qRef2(1,:), 'r');
plot(tSpan, 180/pi*theta34(1,:), 'b');
plot(tSpan, 180/pi*theta34(2,:), 'g');
legend('\theta_2','\theta_3', '\theta_4')
xlabel('Time [s]');
ylabel('Angle [\circ]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, qRef2(2,:), 'r');
plot(tSpan, thetaDot34(1,:), 'b');
plot(tSpan, thetaDot34(2,:), 'g');
hl = legend('$\dot{\theta_2}$','$\dot{\theta_3}$', '$\dot{\theta_4}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time [s]');
ylabel('Angular velocity [rad/s]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, qRef2(3,:), 'r');
plot(tSpan, thetaDotDot34(1,:), 'b');
plot(tSpan, thetaDotDot34(2,:), 'g');
hl = legend('$\ddot{\theta_2}$','$\ddot{\theta_3}$', '$\ddot{\theta_4}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time [s]');
ylabel('Angular acceleration [rad/s^2]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, Rg2);
plot(tSpan, Rg3);
plot(tSpan, Rg4);
legend('x_2', 'y_2', 'x_3', 'y_3', 'x_4', 'y_4');
xlabel('Time [s]');
ylabel('Position [m]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, Vg2);
plot(tSpan, Vg3);
plot(tSpan, Vg4);
hl = legend('$\dot{x_2}$', '$\dot{x_2}$', '$\dot{x_3}$', '$\dot{y_3}$', '$\dot{x_4}$', '$\dot{y_4}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time [s]');
ylabel('Velocity [m/s]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, Ag2);
plot(tSpan, Ag3);
plot(tSpan, Ag4);
hl = legend('$\ddot{x_2}$', '$\ddot{x_2}$', '$\ddot{x_3}$', '$\ddot{y_3}$', '$\ddot{x_4}$', '$\ddot{y_4}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');

sgtitle('Driven');

% Torque

figure;
plot(tSpan, lambdaT(9,:), 'r');
legend('T')
xlabel('Time [s]');
ylabel('Torque [Nm]');
title('Torque');


% Simulated

plotIndex = 0;
plotRows = 2;
plotCols = 2;
figure;
plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, 180/pi*qT(:,3), 'r');
plot(tSpan, 180/pi*qT(:,6), 'b');
plot(tSpan, 180/pi*qT(:,9), 'g');
legend('\theta_2','\theta_3', '\theta_4')
xlabel('Time [s]');
ylabel('Angle [\circ]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, qT(:,12), 'r');
plot(tSpan, qT(:,15), 'b');
plot(tSpan, qT(:,18), 'g');
hl = legend('$\dot{\theta_2}$','$\dot{\theta_3}$', '$\dot{\theta_4}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time [s]');
ylabel('Angular velocity [rad/s]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, qT(:,1:2));
plot(tSpan, qT(:,4:5));
plot(tSpan, qT(:,7:8));
legend('x_2', 'y_2', 'x_3', 'y_3', 'x_4', 'y_4')
xlabel('Time [s]');
ylabel('Position [m]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, qT(:,10:11));
plot(tSpan, qT(:,13:14));
plot(tSpan, qT(:,16:17));
hl = legend('$\dot{x_2}$', '$\dot{y_2}$', '$\dot{x_3}$', '$\dot{y_3}$', '$\dot{x_4}$', '$\dot{y_4}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time [s]');
ylabel('Velocity [m/s]');

sgtitle('Simulated');

%% Animate driven
if animate
    tEnd = inf;
    stp = 4;
    spd = 1;

    figHandle = figure;
    xlabel('X [m]')
    ylabel('Y [m]')
    figHandle.WindowState = 'maximized';
    axHandle = gca;
    axis equal;
    axis([-1, 10, -1, 4]);
    hold on;
    pause(2);
    for i = 1:stp:N
        tic;
        cla(axHandle);
        title(axHandle, ['Driven - t: ', num2str(tSpan(i), '%.1f'), '[s] - Torque: ', num2str(lambdaT(9, i), '%.3f'), ' [Nm]']);
        
        drawBox(axHandle, [0; 0], 0, 0.2, 0.2, 'k', 4);
        drawBox(axHandle, [L1; 0], 0, 0.2, 0.2, 'k', 4);
        drawBox(axHandle, Rg2(:, i), qRef2(1,i), 0.2, L2, 'r', 4);
        drawBox(axHandle, Rg3(:, i), theta34(1,i), 0.2, L3, 'b', 4);
        drawBox(axHandle, Rg4(:, i), theta34(2,i), 0.2, L4, 'g', 4);

        if tSpan(i) >= tEnd
            break;
        end
        delta = toc;
        pause(dt*stp/spd - delta);
    end
end

%% Animate sim
if animateSim
    tEnd = inf;
    stp = 4;
    spd = 1;

    figHandle = figure;
    title('Simulated');
    xlabel('X [m]')
    ylabel('Y [m]')
    figHandle.WindowState = 'maximized';
    axHandle = gca;
    axis equal;
    axis([-1, 10, -1, 4]);
    hold on;
    pause(2);
    for i = 1:stp:N
        tic;
        cla(axHandle);
        title(axHandle, ['Simulated - t: ', num2str(tSpan(i), '%.1f'), '[s] - Torque: ', num2str(lambdaT(9, i), '%.3f'), ' [Nm]']);
        
        drawBox(axHandle, [0; 0], 0, 0.2, 0.2, 'k', 4);
        drawBox(axHandle, [L1; 0], 0, 0.2, 0.2, 'k', 4);
        drawBox(axHandle, qT(i,1:2)', qT(i,3), 0.2, L2, 'r', 4);
        drawBox(axHandle, qT(i,4:5)', qT(i,6), 0.2, L3, 'b', 4);
        drawBox(axHandle, qT(i,7:8)', qT(i,9), 0.2, L4, 'g', 4);

        if tSpan(i) >= tEnd
            break;
        end
        delta = toc;
        pause(dt*stp/spd - delta);
    end
end

