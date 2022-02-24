clc; clear;

%% Parameters
% First order
K1 = 3;
tau = 0.3;

% Second order
K2 = 250;
wc = 10;
z = 0.3;

% Simulation time step
dt = 1e-5;
% Simulation end time
tEnd = 2;
% Time to start the step/ramp
t0 = 0.2;
% Time to end the ramp
t1 = 0.7;
% Step/ramp amplitude
A = 1;

%% Setup
t = 0:dt:tEnd;
N = length(t);

us = zeros(1, N);
ur = zeros(1, N);
ys1 = zeros(1, N);
dys1 = zeros(1, N);
yr1 = zeros(1, N);
dyr1 = zeros(1, N);
ys2 = zeros(1, N);
dys2 = zeros(1, N);
ddys2 = zeros(1, N);
yr2 = zeros(1, N);
dyr2 = zeros(1, N);
ddyr2 = zeros(1, N);

%% Simulate
for i = 2:N
    % Current time
    tc = t(i);
    
    if (tc < t0)
        us(i) = 0;
        ur(i) = 0;
    elseif (tc < t1)
        us(i) = A;
        ur(i) = (tc - t0)/(t1 - t0)*A;
    else
        us(i) = A;
        ur(i) = A;
    end
    
    [ys1(i), dys1(i)] = firstOrder(ys1(i - 1), us(i), K1, tau, dt);
    [yr1(i), dyr1(i)] = firstOrder(yr1(i - 1), ur(i), K1, tau, dt);
    
    [ys2(i), dys2(i), ddys2(i)] = secondOrder(ys2(i - 1), dys2(i - 1), us(i), K2, wc, z, dt);
    [yr2(i), dyr2(i), ddyr2(i)] = secondOrder(yr2(i - 1), dyr2(i - 1), ur(i), K2, wc, z, dt);
    
end


%% Plot
close all;
yMaxRight = max([ys1, yr1])*1.1;
yFinalRight = ys1(end);
yMaxLeft = yMaxRight/yFinalRight;

figure;
title(['First order, K = ', num2str(K1), ', \tau = ', num2str(tau)]);
hold on

yyaxis left
plot(t, [us; ur]);
xlabel('Time [s]');
ylabel('Reference [-]');
axis([0, tEnd, 0, yMaxLeft]);

yyaxis right
plot(t, [ys1; yr1]);
ylabel('Response [-]');
legend('Step reference', 'Ramp reference', 'Step response', 'Ramp response');
axis([0, tEnd, 0, yMaxRight]);


yMaxRight = max([ys2, yr2])*1.1;
yFinalRight = ys2(end);
yMaxLeft = yMaxRight/yFinalRight;

figure;
title(['Second order, K = ', num2str(K2), ', \omega_c = ', ...
    num2str(wc), ', \zeta = ', num2str(z)]);
hold on;

yyaxis left
plot(t, [us; ur]);
xlabel('Time [s]');
ylabel('Reference [-]');
axis([0, tEnd, 0, yMaxLeft]);

yyaxis right
plot(t, [ys2; yr2]);
ylabel('Response [-]');
legend('Step reference', 'Ramp reference', 'Step response', 'Ramp response');
axis([0, tEnd, 0, yMaxRight]);

function [y, dy] = firstOrder(y_in, u, K, tau, dt)
    dy = (K*u - y_in)/tau;
    y = y_in + dy*dt;
end

function [y, dy, ddy] = secondOrder(y_in, dy_in, u, K, wc, z, dt)
    ddy = wc^2*u*K - 2*z*wc*dy_in - wc^2*y_in;
    dy = dy_in + ddy*dt;
    y = y_in + dy*dt;
end

