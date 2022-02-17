% [Time, qT] = ode45(@odeFunct, tSpan, q);
% qDot = odeFunct(t, q)
function [T, qT] = fwdEuler(callbackFunct, tSpan, qInit)
    dt = tSpan(2) - tSpan(1);
    nSamples = length(tSpan);
    nStates = length(qInit);
    qT = zeros(nSamples, nStates);
    qT(1,:) = qInit';
    for i = 2:nSamples
        dq = callbackFunct(tSpan(i - 1), qT(i - 1, :)');
        qT(i, :) = qT(i - 1, :) + dq'*dt;
    end
    
    T = tSpan;
end