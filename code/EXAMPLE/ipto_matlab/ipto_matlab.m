% ipto_matlab.m        March 31, 2015

% two impulse ballistic interplanetary trajectory optimization

% JPL DE424 ephemeris and 64 bit SNOPT algorithm (Sep 15, 2014 version)

% DE424 valid from 12/24/1999 to 2/1/2200

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global iephem km ephname eq2000

global smu ip1 ip2 jdate0 jdate1

global dtr otype imcon revmax au boev

global dv1 dv2 dvm1 dvm2 ri vi rf vf vito

% initialize jpl ephemeris

ephname = 'de424.bin';

iephem = 1;

km = 1;

% J2000 equatorial-to-ecliptic transformation matrix

eq2000 = [[1.000000000000000 0  0]; ...
    [0   0.917482062069182   0.397777155931914]; ...
    [0  -0.397777155931914   0.917482062069182]];

% angular conversion factors

rtd = 180.0d0 / pi;

dtr = pi / 180.0d0;

atr = dtr / 3600.0d0;

% gravitational constant of the sun (km^3/sec^2)

smu = 132712440040.944;

% define "reference" julian date (1/1/2000)

jdate0 = 2451544.5;

% define celestial body name vector

pname = ['Mercury       '; 'Venus         '; 'Earth         '; 'Mars          '; ...
    'Jupiter       '; 'Saturn        '; 'Uranus        '; 'Neptune       '; ...
    'Pluto         '; 'asteroid/comet'];

% begin simulation

clc; home;

fprintf('\n            program ipto_matlab\n');

fprintf('\n< interplanetary trajectory optimization >\n\n');

% request departure calendar date

fprintf('\ndeparture conditions - start date\n');

[month, day, year] = getdate;

jdate1 = julian(month, day, year);

while(1)
    
    fprintf('\nplease input the departure date search boundary in days\n');
    
    ddays1 = input('? ');
    
    if (ddays1 >= 0)
        break;
    end
    
end

% request arrival calendar date

fprintf('\n\narrival conditions - start date\n');

[month, day, year] = getdate;

jdate2 = julian(month, day, year);

% request search information

while(1)
    
    fprintf('\nplease input the arrival date search boundary in days\n');
    
    ddays2 = input('? ');
    
    if (ddays2 >= 0)
        break;
    end
    
end

% request departure and arrival pecl2000s

for i = 1:1:2
    
    fprintf('\n celestial body menu\n');
    
    fprintf('\n  <1>  Mercury');
    fprintf('\n  <2>  Venus');
    fprintf('\n  <3>  Earth');
    fprintf('\n  <4>  Mars');
    fprintf('\n  <5>  Jupiter');
    fprintf('\n  <6>  Saturn');
    fprintf('\n  <7>  Uranus');
    fprintf('\n  <8>  Neptune');
    fprintf('\n  <9>  Pluto');
    fprintf('\n  <10> asteroid/comet');
    
    if (i == 1)
        
        while(1)
            fprintf('\n\nplease select the departure celestial body\n');
            
            ip1 = input('? ');
            
            if (ip1 >= 1 && ip1 <= 10)
                break;
            end
        end
    end
    
    if (i == 2)
        
        while(1)
            
            fprintf('\n\nplease select the arrival celestial body\n');
            
            ip2 = input('? ');
            
            if (ip2 >= 1 && ip2 <= 10)
                break;
            end
        end
    end
end

if (ip2 == 10)
    
    [datafile, pathname] = uigetfile('*.dat', 'please input the name of the asteroid/comet data file');
    
    fid = fopen(datafile, 'r');
    
    % read 25 lines of information
    
    for i = 1:1:25
        
        cline1 = fgetl(fid);
        
        switch i
            
            case 10
                
                tl = size(cline1);
                
                ci = strfind(cline1, ',');
                
                % extract month, day and year
                
                month = str2double(cline1(1:ci(1)-1));
                
                day = str2double(cline1(ci(1)+1:ci(2)-1));
                
                year = str2double(cline1(ci(2)+1:tl(2)));
                
            case 13
                
                boev(1) = str2double(cline1);
                
            case 16
                
                boev(2) = str2double(cline1);
                
            case 19
                
                boev(3) = str2double(cline1);
                
            case 22
                
                boev(4) = str2double(cline1);
                
            case 25
                
                boev(5) = str2double(cline1);
                
        end
        
    end
    
    % julian date of perihelion passage
    
    boev(6) = julian(month, day, year);
    
    fclose(fid);
end

% less than one rev transfer

revmax = 0;

% mission constraints option

while(1)
    
    fprintf('\n\nwould you like to enforce mission constraints (y = yes, n = no)\n');
    
    slct = input('? ', 's');
    
    if (slct == 'y' || slct == 'n')
        break;
    end
    
end

imcon = 0;

if (slct == 'y')
    
    imcon = 1;
    
    while(1)
        
        fprintf('\nplease input the lower bound for departure C3 (kilometers^2/second^2)\n');
        
        flow_c3 = input('? ');
        
        if (flow_c3 > 0)
            break;
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the upper bound for departure C3 (kilometers^2/second^2)\n');
        
        fupp_c3 = input('? ');
        
        if (fupp_c3 > flow_c3)
            break;
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the lower bound for departure DLA (degrees)\n');
        
        flow_dla = input('? ');
        
        if (flow_dla >= -90.0)
            break;
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the upper bound for departure DLA (degrees)\n');
        
        fupp_dla = input('? ');
        
        if (fupp_dla <= +90.0)
            break;
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the lower bound for time-of-flight (days)\n');
        
        flow_tof = input('? ');
        
        if (flow_tof > 0)
            break;
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the upper bound for time-of-flight (days)\n');
        
        fupp_tof = input('? ');
        
        if (fupp_tof > flow_tof)
            break;
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the lower bound for arrival v-infinity (kilometers/second)\n');
        
        flow_vinf = input('? ');
        
        if (flow_vinf > 0)
            break;
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the upper bound for arrival v-infinity (kilometers/second)\n');
        
        fupp_vinf = input('? ');
        
        if (fupp_vinf > flow_vinf)
            break;
        end
        
    end
    
end

% request type of optimization

while(1)
    
    fprintf('\n      optimization menu\n');
    
    fprintf('\n <1> minimize departure delta-v\n');
    
    fprintf('\n <2> minimize arrival delta-v\n');
    
    fprintf('\n <3> minimize total delta-v\n');
    
    fprintf('\n <4> no optimization\n');
    
    fprintf('\n selection (1, 2, 3 or 4)\n');
    
    otype = input('? ');
    
    if (otype == 1 || otype == 2 || otype == 3 || otype == 4)
        break;
    end
    
end

if (otype < 4)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % find optimal solution %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    xg(1) = jdate1 - jdate0;
    
    xg(2) = jdate2 - jdate0;
    
    xg = xg';
    
    % bounds on control variables
    
    xlwr(1) = xg(1) - ddays1;
    xupr(1) = xg(1) + ddays1;
    
    xlwr(2) = xg(2) - ddays2;
    xupr(2) = xg(2) + ddays2;
    
    xlwr = xlwr';
    xupr = xupr';
    
    % bounds on objective function
    
    flow(1) = 0.0d0;
    fupp(1) = +Inf;
    
    if (imcon == 1)
        
        % bounds on nonlinear constraints
        
        flow(2) = flow_c3;
        fupp(2) = fupp_c3;
        
        flow(3) = dtr * flow_dla;
        fupp(3) = dtr * fupp_dla;
        
        flow(4) = flow_tof;
        fupp(4) = fupp_tof;
        
        flow(5) = flow_vinf;
        fupp(5) = fupp_vinf;
        
        fmul = zeros(5, 1);

        fstate = zeros(5, 1);
    
    else
        
        fmul = zeros(1, 1);

        fstate = zeros(1, 1);   
        
    end
    
    flow = flow';
    fupp = fupp';
    
    xmul = zeros(2, 1);

    xstate = zeros(2, 1);
    
    % find optimum
    
    snscreen on;
    
    [x, ~, inform, xmul, fmul] = snopt(xg, xlwr, xupr, xmul, xstate, ...
        flow, fupp, fmul, fstate, 'iptofunc');
        
    if (imcon == 1 && inform ~= 1)
        
        fprintf('\n\ncheck solution!!\n');
        
        fprintf('\n\nall mission constraints may not be satisfied\n');
        
    end
    
    % solution julian dates
    
    jdate1 = x(1) + jdate0;
    
    jdate2 = x(2) + jdate0;
    
    taud = jdate2 - jdate1;
    
    % evaluate current solution
    
    [f, g] = iptofunc (x);
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % no optimization - solve TPBVP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    taud = jdate2 - jdate1;
    
    tof = taud * 86400.0;
    
    % compute initial state vector
    
    [ri, vi] = becl2000(ip1, jdate1);
    
    % compute final state vector
    
    [rf, vf] = becl2000(ip2, jdate2);
    
    % solve Lambert's problem
    
    sv1(1:3) = ri;
    
    sv1(4:6) = vi;
    
    sv2(1:3) = rf;
    
    sv2(4:6) = vf;
    
    [vito, vfto] = glambert(smu, sv1, sv2, tof, revmax);
    
    % calculate departure delta-v
    
    dv1(1) = vito(1) - vi(1);
    dv1(2) = vito(2) - vi(2);
    dv1(3) = vito(3) - vi(3);
    
    dvm1 = norm(dv1);
    
    % calculate arrival delta-v
    
    dv2(1) = vf(1) - vfto(1);
    dv2(2) = vf(2) - vfto(2);
    dv2(3) = vf(3) - vfto(3);
    
    dvm2 = norm(dv2);
    
end

%%%%%%%%%%%%%%%%%
% print results %
%%%%%%%%%%%%%%%%%

% convert solution julian dates to calendar dates and ut

[cdstr1, utstr1] = jd2str(jdate1);

[cdstr2, utstr2] = jd2str(jdate2);

fprintf('\n\n       program ipto_matlab\n');

switch otype
    
    case 1
        
        fprintf('\n   minimize departure delta-v\n');
        
    case 2
        
        fprintf('\n    minimize arrival delta-v\n');
        
    case 3
        
        fprintf('\n     minimize total delta-v\n');
        
    case 4
        
        fprintf('\n        no optimization\n');
        
end

fprintf('\ndeparture celestial body     ');

disp(pname(ip1, 1:14));

fprintf('\ndeparture calendar date      ');

disp(cdstr1);

fprintf('departure TDB time           ');

disp(utstr1);

fprintf('\ndeparture julian date        %10.4f', jdate1);

fprintf('\n\n\narrival celestial body       ');

disp(pname(ip2, 1:14));

fprintf('\narrival calendar date        ');

disp(cdstr2);

fprintf('arrival TDB time             ');

disp(utstr2);

fprintf('\narrival julian date          %10.4f', jdate2');

fprintf('\n\ntransfer time                %10.4f  days \n ', taud);

% print orbital elements of the initial orbit

fprintf('\nheliocentric orbital conditions prior to the first maneuver');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev = eci2orb1 (smu, ri, vi);

oeprint2(smu, oev);

svprint(ri, vi);

% print orbital elements of the transfer orbit

fprintf('\nheliocentric orbital conditions after the first maneuver');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev = eci2orb1 (smu, ri, vito');

oeprint2(smu, oev);

svprint(ri, vito);

rito = ri;

fprintf('\nheliocentric orbital conditions prior to the second maneuver');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

tau = 86400 * (jdate2 - jdate1);

[rft, vft] = twobody2(smu, tau, ri, vito');

oev = eci2orb1 (smu, rft, vft);

oeprint2(smu, oev);

svprint(rft, vft);

% print orbital elements of the final orbit

fprintf('\nheliocentric orbital conditions after the second maneuver');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev = eci2orb1 (smu, rft, vft + dv2);

oeprint2(smu, oev);

svprint(rft, vft + dv2);

% orbital elements and state vector of the target body

[rfb, vfb] = becl2000(ip2, jdate2);

fprintf('\nheliocentric orbital conditions of arrival body');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev = eci2orb1 (smu, rfb, vfb);

oeprint2(smu, oev);

svprint(rfb, vfb);

% compute departure velocity vector in equatorial frame

dveq1 = eq2000' * dv1';

dveqm1 = norm(dveq1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute orientation of the departure hyperbola
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

decl1 = 90.0d0 - rtd * acos(dveq1(3) / dveqm1);

rasc1 = rtd * atan3(dveq1(2), dveq1(1));

% compute arrival v-infinity velocity vector in equatorial frame

dveq2 = eq2000' * (-dv2');

dveqm2 = norm(dveq2);

dveq2prt = dveq2;

if (ip2 == 4)
    
    % Mars mean equator and IAU node of epoch
    
    tmatrix = mme2000 (jdate2);
    
    dvwrk = dveq2;
    
    dveq2 = tmatrix * dvwrk;
    
end

% compute orientation of the arrival hyperbola

decl2 = 90.0d0 - rtd * acos(dveq2(3) / dveqm2);

rasc2 = rtd * atan3(dveq2(2), dveq2(1));

fprintf('\ndeparture delta-v and energy requirements');
fprintf('\n(mean equator and equinox of J2000)');
fprintf('\n-----------------------------------\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * dveq1(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * dveq1(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * dveq1(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second', 1000.0 * dveqm1);

fprintf('\n\nenergy                      %12.6f  kilometers^2/second^2\n', dveqm1 * dveqm1);

fprintf('\nasymptote right ascension   %12.6f  degrees', rasc1);

fprintf('\nasymptote declination       %12.6f  degrees\n\n', decl1);

fprintf('\narrival delta-v and energy requirements');
fprintf('\n(mean equator and equinox of J2000)');
fprintf('\n-----------------------------------\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * dveq2prt(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * dveq2prt(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * dveq2prt(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second', 1000.0 * dveqm2);

fprintf('\n\nenergy                      %12.6f  kilometers^2/second^2\n', dveqm2 * dveqm2);

if (ip2 == 4)
    
    fprintf('\nMars-mean-equator and IAU node of epoch\n');
    
end

fprintf('\nasymptote right ascension   %12.6f  degrees', rasc2);

fprintf('\nasymptote declination       %12.6f  degrees\n', decl2);

fprintf('\n\ntotal delta-v               %12.6f  meters/second', 1000.0 * (dveqm1 + dveqm2));

fprintf('\n\ntotal energy                %12.6f  kilometers^2/second^2\n', ...
    dveqm1 * dveqm1 + dveqm2 * dveqm2);

%%%%%%%%%%%%
% graphics %
%%%%%%%%%%%%

while(1)
    
    fprintf('\n\nwould you like to plot this trajectory (y = yes, n = no)\n');
    
    slct = input('? ', 's');
    
    if (slct == 'y' || slct == 'n')
        break;
    end
    
end

if (slct == 'y')
    
    while(1)
        
        fprintf('\nplease input the plot step size (days)\n');
        
        deltat = input('? ');
        
        if (deltat > 0)
            break;
        end
        
    end
    
    % compute orbital periods
    
    [r1, v1] = becl2000(ip1, jdate1);
    
    oev1 = eci2orb1(smu, r1, v1);
    
    period1 = 2 * pi * oev1(1) * sqrt(oev1(1) / smu) / 86400;
    
    [r2, v2] = becl2000(ip2, jdate1);
    
    oev2 = eci2orb1(smu, r2, v2);
    
    period2 = 2 * pi * oev2(1) * sqrt(oev2(1) / smu) / 86400;
    
    xve = oev1(1) / au;
    
    if (oev2(1) > oev1(1))
        
        xve = oev2(1) / au;
        
    end
    
    % determine number of data points to plot
    
    npts1 = fix(period1 / deltat);
    
    npts2 = fix(period2 / deltat);
    
    npts3 = fix((jdate2 - jdate1) / deltat);
    
    [rti, vti] = orb2eci(smu, oev);
    
    % create departure orbit data points
    
    for i = 0:1:npts1
        
        jdate = jdate1 + i * deltat;
        
        [r1, v1] = becl2000(ip1, jdate);
        
        x1(i + 1) = r1(1) / au;
        
        y1(i + 1) = r1(2) / au;
        
    end
    
    % compute last data point
    
    [r1, v1] = becl2000(ip1, jdate1 + period1);
    
    x1(npts1 + 1) = r1(1) / au;
    
    y1(npts1 + 1) = r1(2) / au;
    
    % create arrival orbit data points
    
    for i = 0:1:npts2
        
        jdate = jdate1 + i * deltat;
        
        [r2, v2] = becl2000(ip2, jdate);
        
        x2(i + 1) = r2(1) / au;
        
        y2(i + 1) = r2(2) / au;
        
    end
    
    % compute last data point
    
    [r2, v2] = becl2000(ip2, jdate1 + period2);
    
    x2(npts2 + 1) = r2(1) / au;
    
    y2(npts2 + 1) = r2(2) / au;
    
    % create transfer orbit data points
    
    for i = 0:1:npts3
        
        tau = 86400 * i * deltat;
        
        [rft, vft] = twobody2(smu, tau, ri, vito');
        
        x3(i + 1) = rft(1) / au;
        
        y3(i + 1) = rft(2) / au;
        
    end
    
    % compute last data point
    
    tau = 86400 * (jdate2 - jdate1);
    
    [rft, vft] = twobody2(smu, tau, ri, vito');
    
    x3(npts3 + 1) = rft(1) / au;
    
    y3(npts3 + 1) = rft(2) / au;
    
    % plot orbits and transfer trajectory
    
    figure(1);
    
    hold on;
    
    plot(x1, y1, '.b');
    
    plot(x1, y1, '-b');
    
    plot(x2, y2, '.g');
    
    plot(x2, y2, '-g');
    
    plot(x3, y3, '.r');
    
    plot(x3, y3, '-r');
    
    % plot and label vernal equinox direction
    
    line ([0.05, 1.1 * xve], [0, 0], 'Color', 'black');
    
    text(1.15 * xve, 0, '\Upsilon');
    
    % label launch and arrival locations
    
    [r2, v2] = becl2000(ip2, jdate1);
    
    plot(r2(1) / au, r2(2) / au, '*r');
    
    text(r2(1) / au + 0.05, r2(2) / au + 0.05, 'L');
    
    [r2, v2] = becl2000(ip1, jdate2);
    
    plot(r2(1) / au, r2(2) / au, '*r');
    
    text(r2(1) / au + 0.05, r2(2) / au + 0.05, 'A');
    
    plot(x3(1), y3(1), '*r');
    
    text(x3(1) + 0.05, y3(1) + 0.05, 'L');
    
    plot(x3(npts3 + 1), y3(npts3 + 1), '*r');
    
    text(x3(npts3 + 1) + 0.05, y3(npts3 + 1) + 0.05, 'A');
    
    % label launch and arrival dates
    
    text(0.85 * xve, -xve + 0.8, 'Launch ', 'FontSize', 8);
    
    text(0.875 * xve, -xve + 0.7, pname(ip1, 1:14), 'FontSize', 8);
    
    text(0.875 * xve, -xve + 0.6, cdstr1, 'FontSize', 8);
    
    text(0.85 * xve, -xve + 0.4, 'Arrival', 'FontSize', 8);
    
    text(0.875 * xve, -xve + 0.3, pname(ip2, 1:14), 'FontSize', 8);
    
    text(0.875 * xve, -xve + 0.2, cdstr2, 'FontSize', 8);
    
    % label plot and axes
    
    xlabel('X coordinate (AU)', 'FontSize', 12);
    
    ylabel('Y coordinate (AU)', 'FontSize', 12);
    
    title('Interplanetary Trajectory Optimization', 'FontSize', 16);
    
    % plot aphelion and perihelion of departure planet
    
    oev1(6) = 0.0;
    
    [r, v] = orb2eci(smu, oev1);
    
    plot(r(1) / au, r(2) / au, 'sb');
    
    oev1(6) = pi;
    
    [r, v] = orb2eci(smu, oev1);
    
    plot(r(1) / au, r(2) / au, 'ob');
    
    % plot aphelion and perihelion of arrival planet
    
    oev2(6) = 0.0;
    
    [r, v] = orb2eci(smu, oev2);
    
    plot(r(1) / au, r(2) / au, 'sg');
    
    oev2(6) = pi;
    
    [r, v] = orb2eci(smu, oev2);
    
    plot(r(1) / au, r(2) / au, 'og');
    
    if (ip1 == 3)
        
        % plot line of nodes (Earth is departure planet)
        
        oev2(6) = -oev2(4);
        
        [r, v] = orb2eci(smu, oev2);
        
        x4(1) = r(1) / au;
        
        y4(1) = r(2) / au;
        
        oev2(6) = oev2(6) + pi;
        
        [r, v] = orb2eci(smu, oev2);
        
        x4(2) = r(1) / au;
        
        y4(2) = r(2) / au;
        
        plot(x4, y4, ':g');
        
    end
    
    % plot sun
    
    plot(0, 0, 'hy', 'MarkerSize', 10);
    
    axis equal;
    
    zoom on;
    
    % the next line creates a color eps graphics file with tiff preview
    
    print -depsc -tiff -r300 ipto_matlab.eps
    
end

%%%%%%%%%%%%%%%%%%%%%%%%
% primer vector analysis
%%%%%%%%%%%%%%%%%%%%%%%%

% perform primer vector initialization

tof = 86400.0 * (jdate2 - jdate1);

[pvi, pvdi] = pviniz(smu, tof, rito, vito', dv1, dv2);

% number of graphic data points

npts = 300;

% plot behavior of primer vector magnitude

dt = tof / npts;

for i = 1:1:npts + 1
    
    t = (i - 1) * dt;
    
    if (t == 0)
        
        % initial value of primer magnitude and derivative
        
        pvm = norm(pvi);
        
        pvdm = dot(pvi, pvdi) / pvm;
        
    else
        
        % primer vector and derivative magnitudes at time t
        
        [pvm, pvdm] = pvector(smu, rito, vito', pvi, pvdi, t);
        
    end
    
    % load data array
    
    xp1(i) = t / 86400.0;
    
    yp1(i) = pvm;
    
    yp2(i) = pvdm;
    
end

figure(2);

hold on;

plot(xp1, yp1, '-r', 'LineWidth', 1.5);

plot(xp1, yp1, '.r');

title('Primer Vector Analysis', 'FontSize', 16);

xlabel('simulation time (days)', 'FontSize', 12);

ylabel('primer vector magnitude', 'FontSize', 12);

grid;

print -depsc -tiff -r300 primer.eps;

% plot behavior of magnitude of primer derivative

figure(3);

hold on;

plot(xp1, yp2, '-r');

plot(xp1, yp2, '.r');

title('Primer Vector Analysis', 'FontSize', 16);

xlabel('simulation time (days)', 'FontSize', 12);

ylabel('primer derivative magnitude', 'FontSize', 12);

grid;

print -depsc -tiff -r300 primer_der.eps;
