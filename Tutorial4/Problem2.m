clc; clear;

L1 = 1000;
muy = 1.5;
mub = 2.5;
Sy = 230; %e6;
E = 2.1e5; %11;

% x0 = [2240; 191; 9.15; 18.1]*10^-3 % [ yA; yB; rA; rB ] initial guess
x0 = [1000; 1000; 12; 12]; % [ yA; yB; rA; rB ] initial guess
options = optimset('Display', 'off');

theta = 3.6364;

if isempty(theta)
    N = 100;
    theta = linspace(0, 2*pi, N);
    fvals = zeros(1, N);
    singleVal = false;
else
    N = 1;
    singleVal = true;
end

x_last = x0;

for i = 1:N
    i
    
    phi = theta(i) + 3/2*pi;
    F = 10000*[cos(phi); sin(phi)];
    pars = [L1; F; muy; Sy; mub; E];
    
    [x, fval, exitflag] = fmincon(@(x)ObjFunc(x, pars), x_last, [], [], [], [], [], [], ...
        @(x)ConFun(x, pars), options);

    fvals(i) = fval;
    if exitflag == 1
        x_last = x;
    else
        theta(i)*180/pi
    end       
end

%% Plot
close all;
if ~singleVal
    figure;
    plot(theta*180/pi, fvals);
    axis([0, 360, 0, 10]);
    xlabel('\theta [deg]');
    ylabel('Mass [kg]');
end
% % Verify volume
% [AA, AB, LA, LB, uA, uB] = getGeometry(x(1), x(2), L1, x(3), x(4));
% volume = AA*LA + AB*LB;
% % Verify forces
% [FA, FB, sigmaA, sigmaB, FbA, FbB] = getForces(uA, uB, F, x(3), x(4), LA, LB, AA, AB, E);


function f = ObjFunc(vars, pars)
    yA = vars(1);
    yB = vars(2);
    rA = vars(3);
    rB = vars(4);
    
    L1 = pars(1);
    
    [AA, AB, LA, LB, ~, ~] = getGeometry(yA, yB, L1, rA, rB);
    f = (AA*LA + AB*LB)/1e9*7800;
end

function [c, ceq] = ConFun(vars, pars)
    yA = vars(1);
    yB = vars(2);
    rA = vars(3);
    rB = vars(4);
    
    L1 = pars(1);
    F = pars(2:3);
    muy = pars(4);
    Sy = pars(5);
    mub = pars(6);
    E = pars(7);
    
    [AA, AB, LA, LB, uA, uB] = getGeometry(yA, yB, L1, rA, rB);
    [FA, FB, sigmaA, sigmaB, FbA, FbB] = getForces(uA, uB, F, rA, rB, LA, LB, AA, AB, E);

    % A is in compression if it is negative
%     cA = -FA*mub - FbA;
    if FA < 0
      cA = abs(FA)*mub - FbA;
    else
        cA = 0;
    end
    
    % B is in compression if it is positive
%     cB = FB*mub - FbB;
    if FB > 0
        cB = FB*mub - FbB;
    else
        cB = 0;
    end
    
    c = [-rA; -rB; sigmaA*muy - Sy; sigmaB*muy - Sy; cA; cB];
    ceq = [];
end
