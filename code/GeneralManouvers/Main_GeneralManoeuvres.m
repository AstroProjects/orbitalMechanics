clc; clear; close all;
%% GENERAL MANOEUVRES
%{ 
We want to find: Omega, i, omega, a, e theta 1
Methode imposes: t1, t2
%}
%% INPUT DATA
% Define input data files
NAME_INPUT_DATA = 'EarthMars1' ;
% Load data files
eval(NAME_INPUT_DATA);