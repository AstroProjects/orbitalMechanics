function [ y ] = spher2cart( x )
%PASS from spherical to cartesian given the vector in the form
% x=[r, lon, lat]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y=x(1)*[cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
end

