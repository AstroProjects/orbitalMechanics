%% INTERPLANETARY DATA
%Departure, Arrival and delta t
t1 = [19, 7, 2020]; %day, month, year
t2 = [25, 1, 2021]; %day, month, year
dt = 190; %days

dLamb0 = 29.837; %[º]
dThet0 = 141.683;%[º]

rT_t1 = [0.4537, -0.9094, 0.0000]; %[AU]
rM_t2 = [0.3148, 1.5078, 0.0239]; %[AU]

%%
R = 384400; %[km]
G = 6.67428e-11; %[m^3·kg^-1·s^-2]

%SpaceCraft (SC)
SC.r0 = 300; %[km]
SC.t0 = 0;
SC.v0 = 10.84; %[km/s]
SC.lamb = 37.5; %[º]
SC.Isp = 300; %[s]
SC.g0 = 9.81e-3; %[km/s^2]


%Earth
DP.R = 6371; %[km]
DP.M = 5.972e24; %[kg]
DP.mu = DP.M*G/1000^3; %[km^3/s^2]
%Mars
AP.Re = 66183; %[km]
AP.R = AP.Re/38; %[km]
%Sun
Sun.R = 6.960e8; %[m]
Sun.M = 1.989e30; %[kg]
Sun.mu = Sun.M*G/1000^3; %[km^3/s^2]
%Orbital parameters
O.R_dp = 149.6e6; %[km] Orbital radii of the departure planet
O.R_ap = 227.9e6; %[km] Orbital radii of the arrival planet