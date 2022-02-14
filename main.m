clear; clc;

%% Parameters
global m1 m2 m3 k1 k2 k3 b1 b2 b3 a d c L0 t1 t2 t3 F1 F2;

m1 = 1; % [kg]
m2 = 1; % [kg]
m3 = 8; % [kg]

k1 = 4; % [N/m]
k2 = 4; % [N/m]
k3 = 5; % [N/m]

b1 = 2; % [Ns/m]
b2 = 4; % [Ns/m]
b3 = 2; % [Ns/m]

a = 0.4; % [m]
d = 0.2; % [m]
c = 0.6; % [m]
L0 = 0.5; % [m]

t1 = 10; % [s]
t2 = 20; % [s]
t3 = 30; % [s]

F1 = 10; % [N]
F2 = -15; % [N]

%% Initial Conditions
x1_0 = a;
x2_0 = a + d + c;
x3_0 = a;

xDot1_0 = 0;
xDot2_0 = 0;
xDot3_0 = 0;

q = [x1_0; x2_0; x3_0; xDot1_0; xDot2_0; xDot3_0];

%% Time
tEnd = 100; % [s] Simulation end time
dt = 0.01; % [s] Plot timestep
Tspan = 0:dt:tEnd;

%% Simulation
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[T, qT] = ode45(@ODE_Funct, Tspan, q, options);

qDot = zeros(length(T), 6);
for i = 1:length(T)
    qDot(i, :) = ODE_Funct(T(i), qT(i, :)); 
end

%% Plot
close all;
figure;
plot(T, qT(:, 1:3));
legend('x_1', 'x_2', 'x_3')
xlabel('Time [s]');
ylabel('Position [m]');

figure;
plot(T, qT(:, 4:6))
legend('dx_1', 'dx_2', 'dx_3')
xlabel('Time [s]');
ylabel('Velocity [m/s]');

figure;
plot(T, qDot(:, 4:6))
legend('ddx_1', 'ddx_2', 'ddx_3')
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');

