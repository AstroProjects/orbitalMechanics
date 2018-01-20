clc; clear; close all;
%% GENERAL MANOEUVRES
%{
We want to find: Omega, i, omega, a, e theta 1
Methode imposes: t1, t2
%}
%% INPUT DATA
% Define input data files
NAME_INPUT_DATA = 'EarthMoon1' ;
% Load data files
eval(NAME_INPUT_DATA);

%% DEPARTURE GEOCENTRIC CONIC

p = Probe.r0^2*Probe.v0^2/Earth.mu;
e = Probe.r0*Probe.v0^2/Earth.mu - 1;
r1 = sqrt(R^2 + Moon.Re^2 - 2*R*Moon.Re*cosd(Probe.lamb));
c_theta1 = 1/e*((p-r1)/r1);
v1 = sqrt(Probe.v0^2-2*Earth.mu*(1/Probe.r0-1/r1));
c_gamma1 = Probe.v0*Probe.r0/(v1*r1);
s_beta = Moon.Re/r1*sind(Probe.lamb);

% if e < 1
%     t1 = sqrt
% end
% if e > 1
% end