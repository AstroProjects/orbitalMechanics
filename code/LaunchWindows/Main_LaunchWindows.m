function [t, az] = Main_LaunchWindows(L, rp, V, planet)
%% HELIOCENTRIC ORBITAL PROBE ELEMENTS
%{ 
We want to find the azimuts and times corresponding to a launch window
given the Launch Location, the parking orbit (which is assumed circular)
and the excess velocity for the previously determined heliocentric orbit.

In:
 - Launch Location latitude [rad]
 - periapsis distance (rp) [km]
 - excess speed [km/s]
 - departure planet data mu, eccentrecity
Out:
 - times t1, t2 of the launch window
 - launching azimuth ai, a2 corresponding to the preceding times

Note: the r1 vector is also needed, but is introduced using an external
file
%}
global PLANET_DATA
%% INPUT DATA
% Define input data files
NAME_INPUT_DATA = 'EarthMars1' ;
% Load data files
eval(NAME_INPUT_DATA);
eval(PLANET_DATA);
mu=MU(planet); obl=OBL(planet); R=Radius(planet);

%define the rotation matrices that are going to be used
Rx = @(x) [1 0 0; 0 cos(x) sin(x); 0 -sin(x) cos(x)];
% Ry = @(y) [cos(x) 0 -sin(x); 0 1 0; sin(x) 0 cos(x)]; not necessary
Rz = @(x) [cos(x) sin(x) 0; -sin(x) cos(x) 0; 0 0 1];

%define the inclination
inclination = @(x,y) 2*atan(norm(x*norm(y)-norm(x)*y)/...
                     norm(x*norm(y)+norm(x)*y));
incl=[L, pi/2]; %min and max inclination

%define the longitude
lambda=atan(rT_t1(2)/rT_t1(1));


%% PREPROCESSING
%Compute the velocity of the parking orbit
vc=sqrt(mu/R+rp);

%Hyperbolic departure trajectories parameters
%compute the right ascension and the declination of the center of the hiper
Vbis=Rx(-obl)*Rz(-lambda);
C=[atan(Vbis(2)/Vbis(1))+pi, asin(Vbis(3))];
a=mu/norm(V)^2;
e=1+(norm(V)/vc)^2;
beta=acos(1/e);
b=a*sqrt(e^2-1);
r=rp*sin(beta); %radius of the cone at the orbit altitude

%% CORE
n0=[0 0 1]; %normal to the equatorial plane
%!!!TODO, add local restrictions
az=zeros(length(incl)); t=az;
for i=1:length(incl)
%     n=Ry(incl(i)); no fa falta
    %compute the azimuth for the limits
    az(i)=asin(cos(incl(i))/cos(lat));
    I=sper2cart([rp C])+[r*cos(incl(i)) r*sin(incl(i)) 0];
    %pass I to spher
    %pass from (alpha, delta) to (H, delta)
    %pass L from (az, lat) to (H, delta)
    %time = H2-H1
    
    
    
end

