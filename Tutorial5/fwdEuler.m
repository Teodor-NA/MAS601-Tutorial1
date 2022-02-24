function [T, qT, dqT] = fwdEuler(callbackFunct, tSpan, qInit)
    dt = tSpan(2) - tSpan(1);
    nSamples = length(tSpan);
    nStates = length(qInit);
    qT = zeros(nSamples, nStates);
    dqT = zeros(nSamples, nStates);
    qT(1,:) = qInit';
    for i = 2:nSamples
        dqT(i, :) = callbackFunct(tSpan(i - 1), qT(i - 1, :)');
        qT(i, :) = qT(i - 1, :) + dqT(i, :)'*dt;
    end
    
    T = tSpan;
end