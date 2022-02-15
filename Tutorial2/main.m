clear; clc;

global L1 L2 Fg1 Fg2 L0 k b Minv

animate = true;

%N = 1000;
tLim = 50;

L0 = 0.3;
k = 10;
b = 20;

g = [0; -9.81]; %m/s^2

m1 = 2; % kg 
L1 = 2; %m 
Jbar1 = 1/12*m1*L1^2; 
m2 = 1; % kg 
L2 = 1; %m 
Jbar2 = 1/12*m2*L2^2; 
 
%% Initial condition: 
phi1 = 15*pi/180; % rad 
phi2 = 30*pi/180; % rad
phiDot1 = 0; % rad/s 
phiDot2 = 0; % rad/s 


s1g = L1/2*[cos(phi1); sin(phi1)];
s2g = L1*[cos(phi1); sin(phi1)] +  L2/2*[cos(phi2); sin(phi2)];

M = diag([m1 m1 Jbar1 m2 m2 Jbar2]);
Minv = inv(M);

Fg1 = m1*g;
Fg2 = m2*g;

q = [s1g; phi1; s2g; phi2; zeros(6,1)];

% Tspan = linspace(0,tLim,N);
dt = 0.01;
Tspan = 0:dt:tLim;
N = length(Tspan);

%% Invoking ode45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[T, qT] = ode45(@odeFunc, Tspan, q);
% dt = T(2) - T(1);

%% Plot
close all;
figure;
subplot(2, 2, 1)
plot(T, qT(:, 1:2))
hold on
plot(T, qT(:, 4:5))
legend('x_1', 'y_1', 'x_2', 'y_2')
xlabel('Time (s)')
ylabel('Position (m)')

subplot(2, 2, 2)
plot(T, [qT(:, 3), qT(:, 6)]*180/pi)
legend('\phi_1', '\phi_2')
xlabel('Time (s)')
ylabel('Angle (\circ)')

subplot(2, 2, 3)
plot(T, qT(:, 7:8))
hold on
plot(T, qT(:, 10:11))
hl = legend('$\dot{x_1}$', '$\dot{y_2}$', '$\dot{x_1}$', '$\dot{y_2}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time (s)')
ylabel('Velocity (m/s)')

subplot(2, 2, 4)
plot(T, [qT(:, 9), qT(:, 12)]*180/pi)
hl = legend('$\dot{\phi_1}$', '$\dot{\phi_2}$');
set(hl, 'Interpreter', 'latex');
xlabel('Time (s)')
ylabel('Angular velocity (\circ/s)')

%% Animate
if animate
    tEnd = inf;
    stp = 4;
    spd = 1;

    figHandle = figure;
    xlabel('X [m]')
    xlabel('Y [m]')
    figHandle.WindowState = 'maximized';
    axHandle = gca;
    axis equal;
    axis([-8, 8, -4, 2]);
    hold on;
    pause(1);
    for i = 1:stp:N
        tic;
        cla(axHandle);
        title(axHandle, ['t: ', num2str(T(i), '%0.3f'), '[s]']);
        
        drawBox(axHandle, [0, 0], 0, 0.2, 0.2, 'k', 4);
        drawBox(axHandle, qT(i, 1:2).', qT(i, 3).', 0.2, L1, 'r', 4);
        drawBox(axHandle, qT(i, 4:5).', qT(i, 6).', 0.2, L2, 'b', 4);
        if T(i) >= tEnd
            break;
        end
        delta = toc;
        pause(dt*stp/spd - delta);
    end
end
