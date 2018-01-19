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
lambda1=atan(rT_t1(2)/rT_t1(1));
lambda2=atan(rM_t2(2)/rM_t2(1));
dLamb=lambda2-lambda1;
beta1=asin(rT_t1(3)/norm(rT_t1));
beta2=asin(rM_t2(3)/norm(rM_t2));
inc=atand(tan(beta2)/sin(dLamb)); %[º]

%r1 & r2:
r1=(rT_t1(1)/(cos(beta1)*cos(lambda1)));
r2=(rM_t2(1)/(cos(beta2)*cos(lambda2)));

%Calcul increment the theta (equacions pag. 14 tema5b)
dTheta=acosd(cos(dLamb)*cos(beta2)); %[º]

%Omega calculation (longitud del node ascendent)
Omega=atand(rT_t1(2)/rT_t1(1)); %[º]

if eliptic == true
    disp('Computing eliptic trajectory...');
    %It has been considered that the equations work in degrees, but it is not sure...
    [e, a, theta1]=Computeeliptic(r1,r2,dt,dTheta);
    if theta1>=360.001
        disp('No theta possible');
    end
else
    disp('Computing hyperbolic trajectory...');
    %It has been considered that the equations work in degrees, but it is not sure...
    [e, a, theta1]=Computehyperbolic(r1,r2,dt,dTheta);
    if theta1>=360.001
        disp('No theta possible');
    end
end
%argument del periheli(omega)
omega=-theta1;





