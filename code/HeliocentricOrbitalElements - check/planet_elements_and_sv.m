function [ coe, r, v, jd ] = planet_elements_and_sv( planet_id, year, month, day, hour, minute, second)
%calculation of the state vector of a planet at a given epoch
%    This function calculates the orbital elements and the state 
%    vector of a planet from the date (year, month, day)
%    and universal time (hour, minute, second).
global mu
deg = pi/180;

%...Equation 5.48:
j0 = J0(year, month, day);
ut = (hour + minute/60 + second/3600)/24;
%...Equation 5.47
jd = j0 + ut;
%...Obtain the data for the selected planet from Table 8.1:
[J2000_coe, rates] = planetary_elements(planet_id);
%...Equation 8.104a:
t0 = (jd - 2451545)/36525;
%...Equation 8.104b:
elements = J2000_coe + rates*t0;
a = elements(1);
e = elements(2);
%...Equation 2.61:
h = sqrt(mu*a*(1 - e^2));
%...Reduce the angular elements to within the range 0 - 360 degrees:
incl = elements(3);
RA = zero_to_360(elements(4));
w_hat = zero_to_360(elements(5));
L = zero_to_360(elements(6));
w = zero_to_360(w_hat - RA);
M = zero_to_360((L - w_hat));
%...Algorithm 3.1 (for which M must be in radians)
E = kepler_E(e, M*deg);
%...Equation 3.10 (converting the result to degrees):
TA = zero_to_360...
(2*atan(sqrt((1 + e)/(1 - e))*tan(E/2))/deg);
coe = [h e RA incl w TA a w_hat L M E/deg];
%...Algorithm 4.2 (for which all angles must be in radians):
[r, v] = sv_from_coe([h e RA*deg incl*deg w*deg TA*deg]);
end

