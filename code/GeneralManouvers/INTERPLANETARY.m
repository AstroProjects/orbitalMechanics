clc; clear; close all;
%% INTERPLANETARY
%{
Spacecraft trajectory form SOI of planet1 to that of planet 2
%}
%% INPUT DATA
% Define input data files
NAME_INPUT_DATA = 'EarthMars2' ;
% Load data files
eval(NAME_INPUT_DATA);

%% PLANETARY DEPARTURE
%hyperbolic excess speed
v_inf = sqrt(Sun.mu/O.R_dp)*(sqrt((2*O.R_ap)/(O.R_dp + O.R_ap))-1); %[km/s]

%Space craft speed on circular parking orbit
vc = sqrt(DP.mu/(DP.R + SC.r0)); %[km/s]

%delta-v to step up to the departue hyperbola
delta_v = vc*(sqrt(2+(v_inf/vc)^2)-1); %[km/s]

%perigee of departure hyperbola
rp = DP.R + SC.r0;
beta = acosd(1/(1+(rp*v_inf^2/DP.mu))); %[º]

%Propellant percent of spacecraft mass
perc_mp = 1 - exp(-delta_v/(SC.Isp*SC.g0)); %tan per one