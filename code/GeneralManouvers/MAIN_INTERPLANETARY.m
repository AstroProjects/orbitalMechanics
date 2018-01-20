clc; clear; close all;
%% INTERPLANETARY
%{
Spacecraft trajectory form SOI of planet1 to that of planet 2
%}
%% INPUT DATA
% Define input data files
NAME_INPUT_DATA = 'EarthMars_hyp' ;
% Load data files
eval(NAME_INPUT_DATA);

%% 
global mu
mu = 1.327124e11;
deg = pi/180;
%...Data for planet 1:
planet_id = 3; % (earth)
hour = 0; minute = 0; second = 0;
%...
depart = [planet_id t1(3) t1(2) t1(1) hour minute second];

%...Data for planet 2:
planet_id = 4; % (Mars)
hour = 0; minute = 0; second = 0;
%...
arrive = [planet_id t2(3) t2(2) t2(1) hour minute second];
[planet1, planet2, trajectory] = interplanetary(depart, arrive);

R1 = planet1(1,1:3);
Vp1 = planet1(1,4:6);
jd1 = planet1(1,7);
R2 = planet2(1,1:3);
Vp2 = planet2(1,4:6);
jd2 = planet2(1,7);
V1 = trajectory(1,1:3);
V2 = trajectory(1,4:6);
tof = jd2 - jd1;
%...Use Algorithm 4.1 to find the orbital elements of the
% spacecraft trajectory based on [Rp1, V1]...
coe = coe_from_sv(R1, V1);
% ... and [R2, V2]
coe2 = coe_from_sv(R2, V2);
%...Equations 8.102 and 8.103:
vinf1 = V1 - Vp1;
vinf2 = V2 - Vp2;

%...Echo the input data and output the solution to
% the command window:
fprintf('---------------------------------------------------')
fprintf('\n Example 8.8')
fprintf('\n\n Departure:\n');
[month_name, planet_name] = month_planet_names(depart(3), ...
depart(1));
fprintf('\n Planet: %s', planet_name)
fprintf('\n Year : %g', depart(2))
fprintf('\n Month : %s', month_name)
fprintf('\n Day : %g', depart(4))
fprintf('\n Hour : %g', depart(5))
fprintf('\n Minute: %g', depart(6))
fprintf('\n Second: %g', depart(7))
fprintf('\n\n Julian day: %11.3f\n', jd1)
fprintf('\n Planet position vector (km) = [%g %g %g]', ...
R1(1), R1(2), R1(3))
fprintf('\n Magnitude = %g\n', norm(R1))
fprintf('\n Planet velocity (km/s) = [%g %g %g]', ...
Vp1(1), Vp1(2), Vp1(3))
fprintf('\n Magnitude = %g\n', norm(Vp1))
fprintf('\n Spacecraft velocity (km/s) = [%g %g %g]', ...
V1(1), V1(2), V1(3))
fprintf('\n Magnitude = %g\n', norm(V1))
fprintf('\n v-infinity at departure (km/s) = [%g %g %g]', ...
vinf1(1), vinf1(2), vinf1(3))
fprintf('\n Magnitude = %g\n', norm(vinf1))
fprintf('\n\n Time of flight = %g days\n', tof)
fprintf('\n\n Arrival:\n');
[month_name, planet_name] = month_planet_names(arrive(3), ...
arrive(1));
fprintf('\n Planet: %s', planet_name)
fprintf('\n Year : %g', arrive(2))
fprintf('\n Month : %s', month_name)
fprintf('\n Day : %g', arrive(4))
fprintf('\n Hour : %g', arrive(5))
fprintf('\n Minute: %g', arrive(6))
fprintf('\n Second: %g', arrive(7))
fprintf('\n\n Julian day: %11.3f\n', jd2)
fprintf('\n Planet position vector (km) = [%g %g %g]', ...
R2(1), R2(2), R2(3))
fprintf('\n Magnitude = %g\n', norm(R1))
fprintf('\n Planet velocity (km/s) = [%g %g %g]', ...
Vp2(1), Vp2(2), Vp2(3))
fprintf('\n Magnitude = %g\n', norm(Vp2))
fprintf('\n Spacecraft Velocity (km/s) = [%g %g %g]', ...
V2(1), V2(2), V2(3))
fprintf('\n Magnitude = %g\n', norm(V2))
fprintf('\n v-infinity at arrival (km/s) = [%g %g %g]', ...
vinf2(1), vinf2(2), vinf2(3))
fprintf('\n Magnitude = %g', norm(vinf2))
fprintf('\n\n\n Orbital elements of flight trajectory:\n')
fprintf('\n Angular momentum (kmˆ2/s) = %g', coe(1))
fprintf('\n Eccentricity = %g', coe(2))
fprintf('\n Right ascension of the ascending node')
fprintf(' (deg) = %g', coe(3)/deg)
fprintf('\n Inclination to the ecliptic (deg) = %g', ...
coe(4)/deg)
fprintf('\n Argument of perihelion (deg) = %g', ...
coe(5)/deg)
fprintf('\n True anomaly at departure (deg) = %g', ...
coe(6)/deg)
fprintf('\n True anomaly at arrival (deg) = %g\n', ...
coe2(6)/deg)
fprintf('\n Semimajor axis (km) = %g', coe(7))
if coe(2) < 1
fprintf('\n Period (days) = %g', ...
2*pi/sqrt(mu)*coe(7)^1.5/24/3600)
end
fprintf('\n-----------------------------------------------\n')
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
