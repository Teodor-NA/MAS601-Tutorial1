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
figure;
plot(t, [us; ur; ys1; yr1]);
xlabel('Time [s]');
ylabel('Amplitude [-]');
legend('Step reference', 'Ramp reference', 'Step response', 'Ramp response');
title(['First order, K = ', num2str(K1), ', \tau = ', num2str(tau)]);

figure;
plot(t, [us; ur; ys2; yr2]);
xlabel('Time [s]');
ylabel('Amplitude [-]');
legend('Step reference', 'Ramp reference', 'Step response', 'Ramp response');
title(['Second order, K = ', num2str(K2), ', \omega_c = ', ...
    num2str(wc), ', \zeta = ', num2str(z)]);


function [y, dy] = firstOrder(y_in, u, K, tau, dt)
    dy = (K*u - y_in)/tau;
    y = y_in + dy*dt;
end

function [y, dy, ddy] = secondOrder(y_in, dy_in, u, K, wc, z, dt)
    ddy = wc^2*u*K - 2*z*wc*dy_in - wc^2*y_in;
    dy = dy_in + ddy*dt;
    y = y_in + dy*dt;
end

