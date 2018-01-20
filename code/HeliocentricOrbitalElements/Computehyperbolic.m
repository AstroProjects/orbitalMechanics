function [e, a, theta1]=Computehyperbolic(r1,r2,dt,dThet)
%This function iterates the value of theta1 in order to solve the system of
%three equations. Theory can be seen in page 27-28 of Tema5b.
theta1=0; %[º]
error=100;
while (error>0.1) && (theta1<=360.001)
    theta1=theta1+0.001;
    %Eccentricity calculation. Its value has to higher than 1
    e=(r2-r1)/(r1*cosd(theta1)-r2*cosd(theta1+dThet));
    if e>1 %Check that the orbit is hyperbolic
        % Semieix
        a=(r1*(1+e*cosd(theta1)))/(e^2-1);
        %Temps de viatge
        t2t1=(365.25/(2*pi))*a^(3/2)*( ( (e*sqrt(e^2-1)*sind(theta1+dThet)) / (1+e*cosd(theta1+dThet)) )...
            -log(abs(  ( tand((theta1+dThet)/2)+sqrt( (e+1)/(e-1) ) ) / ( tand((theta1+dThet)/2)-sqrt( (e+1)/(e-1) )   ) ))...
            -  ( ( e*sqrt(e^2-1)*sind(theta1)  ) / ( 1+e*cosd(theta1) ) )...
            + log(abs(  ( tand((theta1)/2)+sqrt( (e+1)/(e-1) ) ) / ( tand((theta1)/2)-sqrt( (e+1)/(e-1) )   )  )) );
        error=abs(t2t1-dt);
    else
        error=1000;
    end
end