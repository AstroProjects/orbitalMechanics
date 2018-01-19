function [pvm, pvdm] = pvector (mu, ri, vi, pvi, pvdi, x)

% primer vector and derivative magnitudes

% required by ipto_matlab.m

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute state transition matrix at current time x

[rf, vf, stm] = stm2(mu, x, ri, vi);

% evaluate primer vector fundamental equation

ppdot = stm * [pvi'; pvdi];

% extract primer vector and primer derivative vector

pv = ppdot(1:3);

pvd = ppdot(4:6);

% compute primer vector magnitude

pvm = norm(pv);

% compute primer derivative magnitude

pvdm = dot(pv, pvd) / pvm;
