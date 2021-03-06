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
% State variables [w, int(ew)] 
q = zeros(2, N);
% Derivative of state variables [d/dt(w), ew]
dq = zeros(2, N);

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

yMax = max(q(1, :));
yMin = min(q(1, :));
yOffs = (yMax - yMin)*0.05;

figure;
plot(t, [wref; q(1,:)]);
axis([t(1), t(end), yMin - yOffs, yMax + yOffs]);
xlabel('Time [s]');
ylabel('Velocity [rad/s]');
legend('Reference', 'Response');

yMax = max(dq(2, :));
yMin = min(dq(2, :));
yOffs = (yMax - yMin)*0.05;

figure;
plot(t, dq(2, :));
xlabel('Time [s]');
ylabel('Velocity [rad/s]');
legend('Error');
axis([t(1), t(end), yMin - yOffs, yMax + yOffs]);

%% Functions
function [dq, Mm] = sim(q, wref, t, Jeff, Meff, Mmin, Mmax, kp, Ti)
    w = q(1);
    ew_int = q(2);

    ew = wref - w;
    
    if Ti == 0
        Mm = kp*ew;
    else
        Mm = kp*ew + kp/Ti*ew_int;
    end
    Mm = lim(Mmin, Mm, Mmax);
    
    dw = (Mm - Meff)/Jeff;
    
    dq = [dw; ew];
end

function y = lim(xMin, x, xMax)
    y = max(xMin, min(x, xMax));
end