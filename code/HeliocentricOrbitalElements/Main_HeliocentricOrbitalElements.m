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
if inc<0
    inc=360+inc; %Change to a positive angle
end
%r1 & r2:
r1=norm(rT_t1);
r2=norm(rM_t2);

%Calcul increment the theta (equacions pag. 14 tema5b)
dTheta=acosd(cos(dLamb)*cos(beta2)); %[º]

%Omega calculation (longitud del node ascendent)
Omega=radtodeg(lambda1); %[º]
if Omega<0
    Omega=360+Omega;  %Change to a positive angle
end

if eliptic == true
    disp('Computing eliptic trajectory...');
    [e, a, theta1]=Computeeliptic(r1,r2,dt,dTheta);
    if theta1>=360.001
        disp('No theta possible');
    elseif theta1<0
            theta1=360+theta1
        end
else
    disp('Computing hyperbolic trajectory...');
    [e, a, theta1]=Computehyperbolic(r1,r2,dt,dTheta);
    if theta1>=360.001
        disp('No theta possible');
    elseif theta1<0
            theta1=360+theta1
        end
end
%argument del periheli(omega)
omega=-theta1;
if omega<0
    omega=omega+360;
end





