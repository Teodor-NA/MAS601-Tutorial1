clear; clc;

%% Parameters
Jeff = 3.7;
Meff = 1635;
Mmin = 0;
Mmax = 2000;

% Start ramp
w0 = 0;
t0 = 0;
% End ramp
w1 = 12;
t1 = 2;
% End simulation
tEnd = 4;
% Step time
dt = 1e-4;
% PI
kp = 10000;
Ti = 0.5;

%% Setup
t = 0:dt:tEnd;
N = length(t);

% Reference [rad/s]
wref = zeros(1, N);
Mm = zeros(1, N);
% State variables [theta, w, int(ew)] (theta is not actually a state variable, but I
% want it for plotting)
q = zeros(3, N);
% Derivative of state variables [w, d/dt(w), ew]
dq = zeros(3, N);

%% Simulation
for i = 2:N
    if (t(i) < t0)
        wref(i) = w0;
    elseif (t(i) < t1)
        wref(i) = w0 + (t(i) - t0)/(t1 - t0)*(w1 - w0);
    else
        wref(i) = w1;
    end
    
    [dq(:, i), Mm(i)] = sim(q(:, i - 1), wref(i), t(i), Jeff, Meff, Mmin, Mmax, kp, Ti);
    q(:, i) = q(:, i - 1) + dq(:, i)*dt;
end

%% Plot
close all;
figure;
plot(t, [wref; q(2,:)]);

% figure;
% plot(t, Mm);

figure;
plot(t, dq(3, :));

%% Functions
function [dq, Mm] = sim(q, wref, t, Jeff, Meff, Mmin, Mmax, kp, Ti)
    w = q(2);
    ew_int = q(3);

    ew = wref - w;
    
    if Ti == 0
        Mm = kp*ew;
    else
        Mm = kp*ew + kp/Ti*ew_int;
    end
    Mm = lim(Mmin, Mm, Mmax);
    
    dw = (Mm - Meff)/Jeff;
    
    dq = [w, dw, ew];
end

function y = lim(xMin, x, xMax)
    y = max(xMin, min(x, xMax));
end