function tmatrix = mme2000 (jdate)

% eme2000-to-Mars-mean-equator and 
% IAU node of epoch transformation

% input

%  xjdate = julian date

% output

%  tmatrix = transformation matrix

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dtr = pi / 180.0d0;

t = (jdate - 2451545.0d0) / 36525.0d0;

% iau 2000 pole orientation

rasc_pole = 317.68143d0 - 0.1061d0 * t;

decl_pole = 52.88650d0 - 0.0609d0 * t;

phat_mars(1) = cos(rasc_pole * dtr) * cos(decl_pole * dtr);

phat_mars(2) = sin(rasc_pole * dtr) * cos(decl_pole * dtr);

phat_mars(3) = sin(decl_pole * dtr);

% unit pole vector in the j2000 system

phat_j2000(1) = 0.0d0;

phat_j2000(2) = 0.0d0;

phat_j2000(3) = 1.0d0;

% iau-defined x direction

x_iau = cross(phat_j2000, phat_mars);

x_iau = x_iau / norm(x_iau);

% y-direction

yhat = cross(phat_mars, x_iau);

yhat = yhat / norm(yhat);

% load elements of transformation matrix

tmatrix(1, 1) = x_iau(1);
tmatrix(1, 2) = x_iau(2);
tmatrix(1, 3) = x_iau(3);

tmatrix(2, 1) = yhat(1);
tmatrix(2, 2) = yhat(2);
tmatrix(2, 3) = yhat(3);

tmatrix(3, 1) = phat_mars(1);
tmatrix(3, 2) = phat_mars(2);
tmatrix(3, 3) = phat_mars(3);


