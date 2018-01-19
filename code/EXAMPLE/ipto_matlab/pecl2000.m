function [r, v] = pecl2000(ntarg, jdate)

% heliocentric planetary state vector

% earth mean ecliptic and equinox of j2000
% coordinate system (jpl ephemeris)

% input

%  jdate = TDB julian date
%  ntarg = "target" body

% output

%  r = position vector (kilometers)
%  v = velocity vector (kilometers/second)

% NOTE: requires equatorial to ecliptic
%       transformation matrix eq2000 via global

% eq2000 = [+1.000000000000d0 -0.000000479966d0  0.000000000000d0]
%          [+0.000000440360d0 +0.917482137087d0 +0.397776982902d0]
%          [-0.000000190919d0 -0.397776982902d0 +0.917482137087d0]
      
% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global eq2000

% call jpl ephemeris

ncent = 11;

result = jplephem (jdate, ntarg, ncent);

% load position and velocity vectors

rtmp = result(1: 3);

vtmp = result(4: 6);

% convert to ecliptic

r = eq2000 * rtmp;

v = eq2000 * vtmp;



