function [e, a, theta1]=Computeeliptic(r1,r2,dt,dThet)
%This function iterates the value of theta1 in order to solve the system of
%three equations. Theory can be seen in page 23-25 of Tema5b.
theta1=0; %[º]
error=100;
while (error>0.001) && (theta1<=360.001)
        theta1=theta1+0.001;
        %Eccentricity calculation. Its value has to be between 0 and 1
        e=(r2-r1)/(r1*cosd(theta1)-r2*cosd(theta1+dThet));
        if (e>=0) && (e<1) %Check that the orbit is eliptic
            % Semieix
            a=(r1*(1+e*cosd(theta1)))/(1-e^2);
            %Temps de viatge
            t2t1=(365.25/(2*pi))*(a^(3/2))*(2*atan(sqrt((1-e)/(1+e))*tand((theta1+dThet)/2))...
                -((e*sqrt(1-e^2)*sind(theta1+dThet))/(1+e*cosd(theta1+dThet)))...
                -2*atan(sqrt((1-e)/(1+e))*tand(theta1/2))+(e*sqrt(1-e^2)*sind(theta1))/(1+e*cosd(theta1)));
            error=abs(t2t1-dt);
        else
            error=1000;
        end
    end