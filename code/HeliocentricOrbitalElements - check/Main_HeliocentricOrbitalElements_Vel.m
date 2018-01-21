clc; clear; close all;
%% HELIOCENTRIC ORBITAL PROBE ELEMENTS
%{
We want to find: Omega, i, omega, a, e theta1
Methode imposes: t1, t2
%}
%% INPUT DATA
% Define input data files
NAME_INPUT_DATA = 'MarsJupiter1' ;
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
dTheta=180+acosd(sin(beta1)*sin(beta2)+cos(beta1)*cos(beta2)*cos(dLamb));%[º]
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

%% Planet velocities
global mu
mu = 1.327124e11;
deg = pi/180;
%...Algorithm 8.1:
[coe, r, v, jd] = planet_elements_and_sv ...
    (planet_id(1), t1(3), t1(3), t1(3), 12, 0, 0);
%...Convert the planet_id and month numbers into names for output:
[month_name, planet_name] = month_planet_names(t1(2), ...
    planet_id(1));
%...Echo the input data and output the solution to
% the command window:
fprintf('---------------------------------------------------')
fprintf('\n Departure planet')
fprintf('\n Planet: %s', planet_name)
fprintf('\t Date: %g/%g/%g', t1(1),t1(2),t1(3))
fprintf('\t\t Time(UT): %gh %gmin %gs', 12,0,0)
fprintf('\t\t Julian day: %11.3f', jd)
fprintf('\n Position vector (km) = [%g %g %g]', r(1), r(2), r(3))
fprintf('\n Velocity (km/s) = [%g %g %g]',v(1), v(2), v(3))

%...Algorithm 8.1:
[coe, r, v, jd] = planet_elements_and_sv ...
    (planet_id(2), t2(3), t2(3), t2(3), 12, 0, 0);
%...Convert the planet_id and month numbers into names for output:
[month_name, planet_name] = month_planet_names(t2(2), ...
    planet_id(2));

fprintf('\n\n Arrival planet')
fprintf('\n Planet: %s', planet_name)
fprintf('\t Date: %g/%g/%g', t2(1),t2(2),t2(3))
fprintf('\t\t Time(UT): %gh %gmin %gs', 12,0,0)
fprintf('\t\t Julian day: %11.3f', jd)
fprintf('\n Position vector (km) = [%g %g %g]', r(1), r(2), r(3))
fprintf('\n Velocity (km/s) = [%g %g %g]',v(1), v(2), v(3))
fprintf('\n---------------------------------------------------\n')

%% Probe velocity
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
