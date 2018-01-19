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
%Define type of trajectory (eliptic or hiperbolic)
eliptic=true;

%Calcul inclinació i (named 'in' because 'i' is usually used when doing
%iterations) (equacions pag. 14 tema5b)
dLamb=atan(rM_t2(2)/rM_t2(1))-atan(rT_t1(2)/rT_t1(1));
beta2=asin(rM_t2(3)/norm(rM_t2));
inc=atand(tan(beta2)/sin(dLamb)); %[º]

%Calcul increment the theta (equacions pag. 14 tema5b)
dTheta=acosd(cos(dLamb)*cos(beta2)); %[º]

%Omega calculation (longitud del node ascendent)
Omega=radtodeg(rT_t1(2));

if eliptic == true
    disp('Computing eliptic trajectory...');
    %It has been considered that the equations work in degrees, but it is not sure...
[e, a, theta1]=Computeeliptic(rT_t1(1),rM_t2(1),dt,dTheta); 
    if theta1>=360.001
        disp('No theta possible');
    end
else
    disp('Computing hyperbolic trajectory...');
    %It has been considered that the equations work in degrees, but it is not sure...
  [e, a, theta1]=Computehyperbolic(rT_t1(1),rM_t2(1),dt,dTheta);
end
%argument del periheli(omega)
omega=-theta1;




