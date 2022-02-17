clear; clc;

global L1 L2 L3 L4 Theta1 M MInv tSpan dt m2 m3 m4 g Theta2Max h

animate = false;
animateSim = false;

%% Init
Theta1 = 0; 
L1 = 5;  % m 

Theta2Max = 1;

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
h = [0; -m2*g; 0; 0; -m3*g; 0; 0; -m4*g; 0;];

dt = 0.01;
tEnd = 20;
tSpan = 0:dt:tEnd;
N = length(tSpan);

theta34Estimate = [15; 300]*pi/180;% estimates

%% Solve driven
pos = zeros(9, N);
vel = zeros(9, N);
acc = zeros(9, N);
lambdaT = zeros(9, N);

L = [L2; L3; L4];

fsolveOpt = optimset('display', 'off');
for i = 1:N
    %% Single loop
    % Inverse Kinematics
    qRef = traj(tSpan(i), Theta2Max);

    if (i == 1)
        theta34 = fsolve(@fourbar, theta34Estimate, fsolveOpt, qRef(1));
    else
        theta34 = fsolve(@fourbar, [pos(6, i - 1); pos(9, i - 1)], fsolveOpt, qRef(1));
    end
    
    [pos(:, i), vel(:, i), acc(:, i)] = invKin(qRef, theta34, L);
    
    % Inverse Dynamics
    lambdaT(:, i) = invDyn(pos(:, i), acc(:, i), M, h, L);
    
end


%% Solve simulated
% Use the first result from driven inverse kinematics as initial conditions
q = [pos(:,1); vel(:,1)];

% Custom tolerance
% odeOpt = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
% [Time, qT] = ode45(@odeFunct, tSpan, q, odeOpt);

% Default tolerance
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
plot(tSpan, 180/pi*pos(3, :), 'r');
plot(tSpan, 180/pi*pos(6,:), 'b');
plot(tSpan, 180/pi*pos(9,:), 'g');
legend('\theta_2','\theta_3', '\theta_4')
xlabel('Time [s]');
ylabel('Angle [\circ]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, vel(3,:), 'r');
plot(tSpan, vel(6,:), 'b');
plot(tSpan, vel(9,:), 'g');
hl = legend('$\dot{\theta_2}$','$\dot{\theta_3}$', '$\dot{\theta_4}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time [s]');
ylabel('Angular velocity [rad/s]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, acc(3,:), 'r');
plot(tSpan, acc(6,:), 'b');
plot(tSpan, acc(9,:), 'g');
hl = legend('$\ddot{\theta_2}$','$\ddot{\theta_3}$', '$\ddot{\theta_4}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time [s]');
ylabel('Angular acceleration [rad/s^2]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, pos(1:2, :));
plot(tSpan, pos(4:5, :));
plot(tSpan, pos(7:8, :));
legend('x_2', 'y_2', 'x_3', 'y_3', 'x_4', 'y_4');
xlabel('Time [s]');
ylabel('Position [m]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, vel(1:2, :));
plot(tSpan, vel(4:5, :));
plot(tSpan, vel(7:8, :));
hl = legend('$\dot{x_2}$', '$\dot{x_2}$', '$\dot{x_3}$', '$\dot{y_3}$', '$\dot{x_4}$', '$\dot{y_4}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time [s]');
ylabel('Velocity [m/s]');

plotIndex = plotIndex + 1;
subplot(plotRows, plotCols, plotIndex);
hold on
plot(tSpan, acc(1:2, :));
plot(tSpan, acc(4:5, :));
plot(tSpan, acc(7:8, :));
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
        drawBox(axHandle, pos(1:2, i), pos(3, i), 0.2, L2, 'r', 4);
        drawBox(axHandle, pos(4:5, i), pos(6, i), 0.2, L3, 'b', 4);
        drawBox(axHandle, pos(7:8, i), pos(9, i), 0.2, L4, 'g', 4);

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

