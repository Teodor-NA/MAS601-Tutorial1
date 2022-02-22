clc; clear;

syms FA FB uAx uAy uBx uBy Fx Fy

eqs = FA*[uAx; uAy] + FB*[uBx; uBy] + [Fx; Fy] == 0

[FAs, FBs] = solve(eqs, [FA; FB])