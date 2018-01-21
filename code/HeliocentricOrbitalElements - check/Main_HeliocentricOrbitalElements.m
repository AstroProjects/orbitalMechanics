clc; clear; close all;
%% HELIOCENTRIC ORBITAL PROBE ELEMENTS
%{
We want to find: Omega, i, omega, a, e theta1
Methode imposes: t1, t2
%}
%% INPUT DATA
% Define input data files
NAME_INPUT_DATA = 'EarthMars2' ;
% Load data files
eval(NAME_INPUT_DATA);
%Define type of trajectory (eliptic or hiperbolic)
eliptic=false;

Omega=0; omega=0; %initialize values

%% Initial calculation with r vector (pag30 t5b):
r1=norm(rT_t1);
r2=norm(rM_t2);
beta1=asin(rT_t1(3)/r1); %[rad]
beta2=asin(rM_t2(3)/r2); %[rad]
lambda1=zeroto2pi(atan2(rT_t1(2),rT_t1(1))); %[rad]
lambda2=zeroto2pi(atan2(rM_t2(2),rM_t2(1))); %[rad]
dLamb=zeroto2pi(lambda2-lambda1); %[rad]

% Potser langle era 3r quad ... pero no :(
% lambda1s=rT_t1(2)/(r1*cos(beta1));
% lambda1c=rT_t1(1)/(r1*cos(beta1));
% lambda2s=rM_t2(2)/(r2*cos(beta2));
% lambda2c=rM_t2(1)/(r2*cos(beta2));
% dLamb=dLamb+pi;


%% Inclination, dTheta and Omega calculation (general case, page 18 t5b)
dTheta=acosd(sin(beta1)*sin(beta2)+cos(beta1)*cos(beta2)*cos(dLamb));%[º]
A=asin(cos(beta2)*sin(dLamb)/sind(dTheta)); %[rad]
inc=acosd(sin(A)*cos(beta1))*beta2/abs(beta2); %[º] multiplied by the sign 
                                  %of beta2 because if beta2<0, then i<0
if inc<0
    Omega=pi;
    omega=pi;
    inc=abs(inc); 
elseif inc > 90 
    inc = 180 - inc;
end
L=asin(tan(beta1)/tand(inc));%[rad] 
sigma=atan(tan(beta1)/cos(A));%[rad]
Omega=Omega+lambda1-L; %[rad]
Omega=radtodeg(Omega);
if Omega<0
    Omega=360+Omega;
end



%% Calculation of e, a and theta1 (eliptic or hyperbolic)
if eliptic == true
    disp('Computing eliptic trajectory...');
    [e, a, theta1]=Computeeliptic(r1,r2,dt,dTheta);
    if theta1>=360.001
        disp('No theta possible');
    elseif theta1<0
        theta1=360+theta1;
    end
else
    disp('Computing hyperbolic trajectory...');
    [e, a, theta1]=Computehyperbolic(r1,r2,dt,dTheta);
    if theta1>360
        disp('No theta possible');
    elseif theta1<0
        theta1=360+theta1;
    end
end

omega=omega+2*pi-(degtorad(theta1)-sigma); %[rad]
omega=radtodeg(omega);


