%% INTERPLANETARY DATA
%Earth to Mars example with hyperbolic trajectory
%Departure, Arrival and delta t
t1 = [06, 3, 2020]; %day, month, year
t2 = [09, 6, 2020]; %day, month, year
dt = 95; %days

dLamb0 = 81.006; %[º]
dThet0 = 135.670;%[º]


rT_t1 = [-0.9609, 0.2466, 0.0000]; %[AU]
rM_t2 = [0.7285, -1.1980, -0.0430]; %[AU]

%rT_t1 = [-1.43752e+08 3.69097e+07 -247.055]; %[km]
%rM_t2 = [1.09014e+08 -1.79191e+08 -6.43273e+06]; %[km]

muS=1.32712440018e20;