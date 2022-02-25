clear; clc;

%% Parameters
Jeff = 3.7;
Meff = 1635;
Mmin = 0;
Mmax = 2000;
ttq = 0.15;

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
kp = 1000;
Ti = 0.3;

%% Setup
t = 0:dt:tEnd;
N = length(t);

% Velocity reference [rad/s]
wref = zeros(1, N);
% Torque reference
Mm_ref = zeros(1, N);
% Velocity error
ew = zeros(1, N);
% Proportional gain
Gp = zeros(1, N);
% Integral gain
Gi = zeros(1, N);
% State variables [w, Mm] 
q = zeros(2, N);
% Derivative of state variables [d/dt(w), dMm]
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
    
    [dq(:, i), Mm_ref(i), ew(i), Gp(i), Gi(i)] = sim(q(:, i - 1), wref(i), t(i), dt, Jeff, Meff, Mmin, Mmax, kp, Ti, Gi(i-1), ttq);
    q(:, i) = q(:, i - 1) + dq(:, i)*dt;
    % Add saturation to Mm
    q(2, i) = lim(Mmin, q(2, i), Mmax);
end

w = q(1, :);
Mm = q(2, :);

%% Plot
close all;

figure;
plot(t, [wref; w]);
axis(autoscale(t, [w, wref]));
xlabel('Time [s]');
ylabel('Velocity [rad/s]');
legend('Reference', 'Response');

figure;
plot(t, ew);
xlabel('Time [s]');
ylabel('Velocity [rad/s]');
legend('Error');
axis(autoscale(t, ew));

figure;
plot(t, [Mm_ref; Mm]);
xlabel('Time [s]');
ylabel('Torque [Nm]');
legend('Reference', 'Motor');
axis(autoscale(t, [Mm, Mm_ref]));

figure;
hold on
plot(t, Gp, 'r');
plot(t, Gi, 'b');
plot(t, Mm_ref, 'g');
xlabel('Time [s]');
ylabel('Torque [Nm]');
legend('Proportional', 'Integral', 'Total');
axis(autoscale(t, [Gp, Gi, Mm_ref]));

%% Functions
function [dq, Mm_ref, ew, Gp, Gi] = sim(q, wref, t, dt, Jeff, Meff, Mmin, Mmax, kp, Ti, Gi_0, ttq)
    w = q(1);
    Mm = q(2);

    ew = wref - w;
    
    Gp = kp*ew;
        
    if Ti == 0
        Gi = 0;
    else
        % Anti-windup
        Gi = lim(Mmin  - Gp, Gi_0 + ew*dt*kp/Ti, Mmax - Gp);
    end
    
    Mm_ref = Gp + Gi;
    
    dMm = 1/ttq*(Mm_ref - Mm);
    
    dw = (Mm - Meff)/Jeff;
    
    dq = [dw; dMm];
end

function y = lim(xMin, x, xMax)
    y = max(xMin, min(x, xMax));
end

function axis = autoscale(x, y)
    yMax = max(y);
    yMin = min(y);
    yOffs = (yMax - yMin)*0.05;
    
    axis = [x(1), x(end), yMin - yOffs, yMax + yOffs];  
end
