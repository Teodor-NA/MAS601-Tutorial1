clear; clc;

%% Parameters
global m1 m2 m3 k1 k2 k3 b1 b2 b3 a d c L0 t1 t2 t3 F1 F2;

visu = true;

m1 = 1; % [kg]
m2 = 1; % [kg]
m3 = 8; % [kg]

k1 = 4; % [N/m]
k2 = 4; % [N/m]
k3 = 5; % [N/m]

b2 = 4; % [Ns/m]
b1 = 15; % [Ns/m]
b3 = 15; % [Ns/m]

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
N = length(T);

qDot = zeros(N, 6);
fric = zeros(N, 2);
f = zeros(1, N);
for i = 1:N
    qDot(i, :) = ODE_Funct(T(i), qT(i, :));
    
    x1 = qT(i, 1);
    x2 = qT(i, 2);
    x3 = qT(i, 3);
    
    if (x3 - x1 > d)
        fric(i, 1) = -1;
    elseif (x1 > x3 + 8*d)
        fric(i, 1) = 1;
    else
        fric(i, 1) = 0;
    end
    if (x3 - x2 > d) 
        fric(i, 2) = -1;
    elseif (x2 > x3 + 8*d)
        fric(i, 2) = 1;
    else
        fric(i, 2) = 0;
    end
    
    if T(i) < t1
        f(i) = F1*T(i)/t1;
    elseif T(i) < t2
        f(i) = F1 + (F2 - F1)/(t2 - t1)*(T(i) - t1);
    elseif T(i) < t3
        f(i) = F2 + (0 - F2)/(t3 - t2)*(T(i) - t2);
    else
        f(i) = 0;
    end
end

%% Plot
close all;
figure;
plot(T, qT(:, 1:3));
hold on
plot(T, fric)
legend('x_1', 'x_2', 'x_3', 'Contact 1', 'Contact 2')
xlabel('Time [s]');
ylabel('Position [m]');

figure;
plot(T, qT(:, 4:6))
hold on
plot(T, fric)
legend('dx_1', 'dx_2', 'dx_3', 'Contact 1', 'Contact 2')
xlabel('Time [s]');
ylabel('Velocity [m/s]');

figure;
plot(T, qDot(:, 4:6))
hold on
plot(T, fric)
legend('ddx_1', 'ddx_2', 'ddx_3', 'Contact 1', 'Contact 2')
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');

figure;
plot(T, fric)
legend('x_1', 'x_2')
xlabel('Time [s]');
ylabel('Contact');

%% Visualize
tEnd = 50;
fps = 60;
stp = (1/fps)/dt
stp = 5;
spd = 2;
if visu
    close all;
    angle = 0;
    ca = cos(angle);
    sa = sin(angle);
    rot = [ca, -sa; sa, ca];

    figHandle = figure;
    xlabel('Position [m]')
    figHandle.WindowState = 'maximized';
    axHandle = gca;
    axis equal;
    axis([-8, 8, -2.5, 2.5]);
    hold on;
    for i = 1:stp:N   
        tic;
        cla(axHandle);
        title(axHandle, ['t: ', num2str(T(i), '%0.3f'), '[s], F: ', num2str(f(i), '%2.2f'), '[N]']);
        
        if fric(i, 1)
            boxColor1 = 'c';
        else
            boxColor1 = 'b';
        end
        if fric(i, 2)
            boxColor2 = 'm';
        else
            boxColor2 = 'r';
        end
        
        drawBox(axHandle, [qT(i, 1) + d/2; 0], rot, d, d, boxColor1, 4);
        drawBox(axHandle, [qT(i, 2) + d/2; 0], rot, d, d, boxColor2, 4);
        drawBox(axHandle, [qT(i, 3) + 4*d; d+0.05], rot, d, 8*d, 'g', 4);
        if T(i) >= tEnd
            break;
        end
        delta = toc;
        pause(dt*stp/spd - delta);
    end
end