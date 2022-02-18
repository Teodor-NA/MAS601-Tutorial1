clc; clear;

x0 = [2; 2];
options = [];
[xa, fvala] = fmincon(@ObjFunc, x0, [], [], [], [], [], [], @(x)ConFun(x, Assignment.a), options)
[xb, fvalb] = fmincon(@ObjFunc, x0, [], [], [], [], [], [], @(x)ConFun(x, Assignment.b), options)
[xc, fvalc] = fmincon(@ObjFunc, x0, [], [], [], [], [], [], @(x)ConFun(x, Assignment.c), options)

function f = ObjFunc(vars, cst)
    x = vars(1);
    y = vars(2);

    f = (x - 5)^2 + (y - 8)^2 + 7;
end

function [c, ceq] = ConFun(vars, cst)
    x = vars(1);
    y = vars(2);
    
    assignment = cst;
    
    if assignment == Assignment.a
        c = [];
        ceq = [];
    elseif assignment == Assignment.b
        c = y - 12 + 2.4*x;
        ceq = [];
    elseif assignment == Assignment.c
        c = [];
        ceq = x - y;
    else
        c = [];
        ceq = [];
    end
end
