function [FA, FB, sigmaA, sigmaB, FbA, FbB] = getForces(uA, uB, F, rA, rB, LA, LB, AA, AB, E)
    M = [uA, uB];
    
    FF = M\(-F);
    
    FA = FF(1);
    FB = FF(2);
    
    sigmaA = abs(FA)/AA;
    sigmaB = abs(FB)/AB;
    
    IA = pi*(rA^4)/4;
    IB = pi*(rB^4)/4;
    
    FbA = (pi^2)*E*IA/(LA^2);
    FbB = (pi^2)*E*IB/(LB^2);
end