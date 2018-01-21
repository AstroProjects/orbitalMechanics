clc; clear; close all;
%% HELIOCENTRIC ORBITAL PROBE ELEMENTS
%{
We want to find: Omega, i, omega, a, e theta1
Methode imposes: t1, t2
%}
%% INPUT DATA
% Define input data files
addpath('cases');
NAME_INPUT_DATA = 'MarsJupiter1v' ;
% Load data files
eval(NAME_INPUT_DATA);
%Define type of trajectory (eliptic or hiperbolic)
eliptic=true;

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
dTheta=180+acosd(sin(beta1)*sin(beta2)+cos(beta1)*cos(beta2)*cos(dLamb));%[º] &MARSJUPITER
%dTheta=acosd(sin(beta1)*sin(beta2)+cos(beta1)*cos(beta2)*cos(dLamb));
A=asin(cos(beta2)*sin(dLamb)/sind(dTheta)); %[rad]
inc=acosd(sin(A)*cos(beta1))*(beta2-beta1)/abs(beta2-beta1); %[º] multiplied by the sign 
                                               %of incbeta because if incbeta<0, then i<0
if inc<0
    Omega=pi;
    omega=pi;
    inc=abs(inc); 
end
L=asin(tan(beta1)/tand(inc));%[rad] 
sigma=atan(tan(beta1)/cos(A));%[rad]
Omega=Omega+lambda1-L; %[rad]
Omega=radtodeg(Omega);
if Omega<0
    Omega=360+Omega;
end

if inc > 90
    inc = 180 - inc;
end


%% Calculation of e, a and theta1 (eliptic or hyperbolic)
if eliptic == true
    disp('Computing eliptic trajectory...');
    [e, a, theta1]=Computeeliptic(r1,r2,dt,dTheta);
    theta1=radtodeg(theta1);
    if theta1>360
        disp('No theta possible');
    elseif theta1<0
        theta1=360+theta1;
    end
else
    disp('Computing hyperbolic trajectory...');
    [e, a, theta1]=Computehyperbolic(r1,r2,dt,dTheta);
    theta1=radtodeg(theta1);
    if theta1>360
        disp('No theta possible');
    elseif theta1<0
        theta1=360+theta1;
    end
end

omega=omega+2*pi-(degtorad(theta1)-sigma); %[rad]
omega=radtodeg(omega);


%% PRINT
fprintf('\n\na = %g AU;\te = %g;\ttheta_0 = %g deg;\nw = %g deg;\ti = %g deg;\tW = %g deg\n',...
    a,e,theta1,omega,inc,Omega);

%% Heliocentric velocities
method=1;
if method==2
    %using the 2nd method described in tema2 (pg59)
    cO=cosd(Omega); sO=sind(Omega);
    co=cosd(omega); so=sind(omega);
    ci=cosd(inc); si=sind(inc);
    P=[cO*co-sO*ci*so, sO*co+cO*ci*so, si*so];
    Q=[-cO*so-sO*ci*co, -sO*so+cO*ci*co si*co];
    p=a*(e^2-1);
    if eliptic; p=-p; end

    v1=1/1e3*sqrt(muS/(p*1.496e11))*(-sind(theta1)*P+(e+cosd(theta1))*Q);
    v2=1/1e3*sqrt(muS/(p*1.496e11))*(-sind(dTheta+theta1)*P+(e+cosd(dTheta+theta1))*Q);
    
    %prova de verificació resultats
    r1bis=p/(1+e*cosd(theta1))*(cosd(theta1)*P+sind(theta1)*Q);
    r2bis=p/(1+e*cosd(theta1+dTheta))*(cosd(theta1+dTheta)*P+sind(theta1+dTheta)*Q);
elseif method==1
    %compute M: see T2 pg18 & pg29 NOT REALLY necessary..... :(
    Meli= @(e,th) 2*atan((tand(th/2)*sqrt((e+1)/(e-1))))-...
          (e*sqrt(1-e^2)*sind(th))/(1+e*cosd(th));
    Mhyp = @(e,th) sqrt(e^2-1)*(e*sind(th)/(1+e*cosd(th))-1/sqrt(e^2-1)*...
          log((tand(th/2)+sqrt((e+1)/(e-1)))/(tand(th/2)-sqrt((e+1)/(e-1)))));
    if eliptic
        M1=Meli(e,theta1);
        M2=Meli(e,theta1+dTheta);
    else
        M1=Mhyp(e,theta1);
        M2=Mhyp(e,theta1+dTheta);
    end
    %prev vars
    p=a*(e^2-1);
    if eliptic; p=-p; end
    r = @(th) p*1.496e11/(1+cos(th));   u = @(th) omega+th;
    cO=cosd(Omega); sO=sind(Omega);
    ci=cosd(inc); si=sind(inc);
    
    r1bis=r(theta1)/1.496e11*[cO*cosd(u(theta1))-sO*ci*sind(u(theta1)) ...
        sO*cosd(u(theta1))+cO*ci*sind(u(theta1)) si*sind(u(theta1))];
    r2bis=r(theta1+dTheta)/1.496e11*[cO*cosd(u(theta1+dTheta))-sO*ci*sind(u(theta1+dTheta)) ...
        sO*cosd(u(theta1+dTheta))+cO*ci*sind(u(theta1+dTheta)) si*sind(u(theta1+dTheta))];
    
    v1=r1bis*e/r(theta1)*sqrt(muS/(p*1.496e11))*sind(theta1)+r(theta1)*...
        [-cO*sind(u(theta1))-sO*ci*cosd(u(theta1))...
        -sO*sind(u(theta1))+cO*ci*cosd(u(theta1))...
        si*cosd(u(theta1))]*sqrt(muS*(p*1.496e11))/r(theta1)^2;
    v2=r2bis*e/r(theta1+dTheta)*sqrt(muS/(p*1.496e11))*sind(theta1+dTheta)+r(theta1+dTheta)...
        *[-cO*sind(u(theta1+dTheta))-sO*ci*cosd(u(theta1+dTheta))...
        -sO*sind(u(theta1+dTheta))+cO*ci*cosd(u(theta1+dTheta))...
        si*cosd(u(theta1+dTheta))]*sqrt(muS*(p*1.496e11))/r(theta1+dTheta)^2;
    v1=v1/1e3;
    v2=v2/1e3;
end
    

fprintf('\n\nVs(t1) = [%g, %g, %g] km/s\nVs(t2) = [%g, %g, %g] km/s\n',...
    v1(1),v1(2),v1(3),v2(1),v2(2),v2(3));

fprintf('\n\nCOMPROVACIÓ\nr1(t1) = [%g, %g, %g] AU\nr2(t2) = [%g, %g, %g] AU\n',...
    r1bis(1),r1bis(2),r1bis(3),r2bis(1),r2bis(2),r2bis(3));
%% Planet velocities