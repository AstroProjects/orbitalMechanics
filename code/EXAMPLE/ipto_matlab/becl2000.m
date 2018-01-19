function [r, v] = becl2000(ntarg, jdate)

% heliocentric planet/asteroid/comet state vector

% earth mean ecliptic and equinox of j2000 coordinate system

% input

%  jdate = TDB julian date
%  ntarg = "target" body

% output

%  r = position vector (kilometers)
%  v = velocity vector (kilometers/second)

% NOTE: requires equatorial to ecliptic
%       transformation matrix eq2000 via global
      
% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global dtr smu au eq2000 boev

if (ntarg == 10)
    
   % asteroid/comet
   
   rper = boev(1);

   ecc = boev(2);

   xinc = dtr * boev(3);

   argper = dtr * boev(4);

   raan = dtr * boev(5);

   jdpp = boev(6);

   % semimajor axis (au)

   sma = rper / (1.0d0 - ecc);

   % time since perhelion passage (seconds)

   tspp = 86400.0d0 * (jdate - jdpp);

   % compute mean anomaly (radians)

   manom = sqrt(smu / abs(au * sma)^3) * tspp;

   % solve Kepler's equation for true anomaly

   [eanom, tanom] = kepler1 (manom, ecc);

   % load orbital elements array

   oev(1) = au * sma;
   oev(2) = ecc;
   oev(3) = xinc;
   oev(4) = argper;
   oev(5) = mod(raan, 2 * pi);
   oev(6) = tanom;

   % determine heliocentric state vector

   [r, v] = orb2eci(smu, oev);

   return
   
end

% call jpl ephemeris

ncent = 11;

result = jplephem (jdate, ntarg, ncent);

% load position and velocity vectors

rtmp = result(1: 3);

vtmp = result(4: 6);

% convert to ecliptic

r = eq2000 * rtmp;

v = eq2000 * vtmp;



