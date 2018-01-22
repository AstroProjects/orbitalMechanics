function [C, r] = Main_LaunchWindows(vp, v1, beta)
%% HELIOCENTRIC ORBITAL PROBE ELEMENTS
%{ 
We want to find the azimuts and times corresponding to a launch window
given the Launch Location latitude, max & min inclination for the parking 
orbit, the parking orbit and the velocity of both the spacecraft and the
departure planet

In:
 - vp:      Planet velocity [km/s] -- heliocentric reference
 - v1:      Probe velocity [km/s] -- heliocentric reference
 - obl:     departure planetecliptic obliquity (introduced with a file)
 - SC.r0:   periapsis distance [km] (introduced with a file)
Out:
 - times t1, t2 of the launch window
 - launching azimuth ai, a2 corresponding to the preceding times
Note: 
%}
%% INPUT DATA
% Define input and anonymous functions
% Define input data files
NAME_INPUT_DATA = 'EarthMars1' ;
% Load data files
eval(NAME_INPUT_DATA);

%define the rotation matrices that are going to be used
Rx = @(x) [1 0 0; 0 cos(x) sin(x); 0 -sin(x) cos(x)];
% Ry = @(y) [cos(x) 0 -sin(x); 0 1 0; sin(x) 0 cos(x)]; not necessary
Rz = @(x) [cos(x) sin(x) 0; -sin(x) cos(x) 0; 0 0 1];

%angle between two vectors
angle = @(x,y) atan2(norm(cross(x,y)),dot(x,y));

%radius of the circle of injection points (CIP) at the spacecarft altitude
r=(SC.r0+DP.R)*sin(beta); %[km]

%excess velocity
v=v1-vp;

%% CORE
%1. Compute the ecliptic longitude between X and X0 (see 5b, pg32)
lambda = angle([1 0 0], vp) - pi/2; %[rad]
%2. Compute the velocity in the equatorial frame (see 5b, pg33)
v_Q = Rx(-obl)*Rz(-lambda)*v';
%declination & right ascension angle of the CIP center
dec = asin(v_Q(3)); if dec<0; dec=dec+2*pi; end %[rad]
rasc = atan2(v_Q(2), v_Q(1)); if rasc<0; rasc=rasc+2*pi; end %[rad]

%% RETURN
C=[rasc*24/(2*pi), dec];

