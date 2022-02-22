function [AA, AB, LA, LB, uA, uB] = getGeometry(yA, yB, L1, rA, rB)
    lA = [-L1; yA];
    lB = [L1; yB];
    
    LA = sqrt(lA.'*lA);
    LB = sqrt(lB.'*lB);
    
    uA = lA/LA;
    uB = lB/LB;
    
    AA = pi*rA^2;
    AB = pi*rB^2;
end
