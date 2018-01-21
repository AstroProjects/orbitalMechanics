%% INTERPLANETARY DATA
%Departure, Arrival and delta t
t1 = [05, 6, 2026]; %day, month, year
t2 = [25, 4, 2029]; %day, month, year
dt = 1055; %days

dLamb0 = 101.387; %[º]
dThet0 = 182.859;%[º]


rT_t1 = [1.3277, 0.4901, 0.0223]; %[AU]
rM_t2 = [5.0135, 2.1380, 0.0505]; %[AU]

%ID
planet_id = [ 4, 5];

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