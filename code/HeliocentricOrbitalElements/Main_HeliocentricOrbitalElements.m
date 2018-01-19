clc; clear; close all;
%% HELIOCENTRIC ORBITAL PROBE ELEMENTS
%{ 
We want to find: Omega, i, omega, a, e theta1
Methode imposes: t1, t2
%}
%% INPUT DATA
% Define input data files
NAME_INPUT_DATA = 'EarthMars1' ;
% Load data files
eval(NAME_INPUT_DATA);

%Calcul inclinació i (named 'in' because 'i' is usually used when doing iterations)
cosinc=tand(dLamb0)/tand(dThet0); 
inc=acosd(cosinc);%[º]

%Omega calculation (longitud del node ascendent)
Omega=radtodeg(rT_t1(2));
if Omega<0
    Omega=360+Omega;
end

theta1=0; %[º]
error=100;
while (error>10) && (theta1<=360.001)
    theta1=theta1+0.001;
    %Eccentricity calculation. Its value has to be between 0 and 1
    e=(rM_t2(1)-rT_t1(1))/(rT_t1(1)*cosd(theta1)-rM_t2(1)*cosd(theta1+dThet0));
    if (e>=0) && (e<1) %Check that the orbit is eliptic
    % Semieix
    a=(rT_t1(1)*(1+e*cosd(theta1)))/(1-e^2);
    %Temps de viatge
    t2t1=(365.25/(2*pi))*(a^(3/2))*(2*atand(sqrt((1-e)/(1+e))*tand((theta1+dThet0)/2))...
        -((e*sqrt(1-e^2)*sind(theta1+dThet0))/(1+e*cosd(theta1+dThet0)))...
        -2*atand(sqrt((1-e)/(1+e))*tand(theta1/2))+(e*sqrt(1-e^2)*sind(theta1))/(1+e*cosd(theta1)));                                       
    error=abs(t2t1-dt);
    else
        error=1000;
    end
end
if theta1>=360.001
    disp('No theta possible');
end




