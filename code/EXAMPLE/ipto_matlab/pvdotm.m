function f = pvdotm(ri, vi, x)

% magnitude of primer derivative vector

% required by primer.m

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global smu pvi pvdi

% compute state transition matrix

[rf, vf, stm] = stm2(smu, x, ri, vi);

% compute magnitude of primer derivative vector

ppdot = stm * [pvi'; pvdi];

pdot = ppdot(4:6);

f = norm(pdot);
