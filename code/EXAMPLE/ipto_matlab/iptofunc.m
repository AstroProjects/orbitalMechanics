function [f, g] = iptofunc (x)
 
% delta-v objective function

% required by ipto_matlab.m

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global smu ip1 ip2 jdate0 eq2000

global otype imcon revmax dv1 dv2

global ri vi rf vf vito vfto sv2

% current julian dates

jdate1 = x(1) + jdate0;

jdate2 = x(2) + jdate0;

% time-of-flight (seconds)

taud = jdate2 - jdate1;

tof = taud * 86400;

% compute initial state vector

[ri, vi] = becl2000(ip1, jdate1);

% compute final state vector

[rf, vf] = becl2000(ip2, jdate2);

% solve Lambert's problem

sv1(1:3) = ri;

sv1(4:6) = vi;
    
sv2(1:3) = rf;

sv2(4:6) = vf;

[vito, vfto] = glambert(smu, sv1, sv2, tof, revmax);

% calculate departure delta-v (kilometers/second)

dv1(1) = vito(1) - vi(1);
dv1(2) = vito(2) - vi(2);
dv1(3) = vito(3) - vi(3);

dvm1 = norm(dv1);

% calculate arrival delta-v (kilometers/second)

dv2(1) = vf(1) - vfto(1);
dv2(2) = vf(2) - vfto(2);
dv2(3) = vf(3) - vfto(3);

dvm2 = norm(dv2);

% load scalar objective function

switch otype
    
case 1
    
   % launch
   
   f(1) = dvm1;
   
case 2
    
   % arrival
   
   f(1) = dvm2;
   
case 3
    
   % launch + arrival
   
   f(1) = dvm1 + dvm2;
   
case 4
    
   f(1) = dvm1 + dvm2;
   
end

if (imcon == 1)
  
   % -----------------------------------
   % compute current mission constraints
   % -----------------------------------
   
   % C3L (km^2/sec^2)
   
   f(2) = norm(dv1) * norm(dv1);

   dveq1 = eq2000' * dv1';

   % DLA (radians)
   
   f(3) = 0.5 * pi - acos(dveq1(3) / norm(dveq1));
   
   % TOF (days)
   
   f(4) = jdate2 - jdate1;
   
   % arrival v-infinity (km/sec)
   
   f(5) = norm(dv2);
   
end

f = f';

% no derivatives

g = [];
