function [ j0 ] = J0( year, month, day )
%JULIAN day number J0 in days
%   Julian day number at 0 hr UT 
% 1901 ? year ? 2099
% 1 ? month ? 12
% 1 ? day ? 31

j0 = 367*year - fix(7*(year+fix((month+9)/12))/4)+fix(275*month/9)+day+1721013.5;

end