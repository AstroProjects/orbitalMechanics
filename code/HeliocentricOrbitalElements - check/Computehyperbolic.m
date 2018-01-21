function [e, a, theta1]=Computehyperbolic(r1,r2,dt,dThet)
%This function iterates the value of theta1 in order to solve the system of
%three equations. Theory can be seen in page 27-28 of Tema5b.
%% PREPROCESS
dTheta=degtorad(dThet); %[rad]
theta1 = degtorad(0); %[rad]

error=100;
delta=1e-5;
dtheta=1e-4;

%avoid looping between 2 values
errorArray=zeros(2,1);
loopVal=0;
maxLoopVal=3;

%% LOOP
while (error>delta) && (theta1<=2*pi)
    %1. Eccentricity calculation. Its value has to be higher than 1
    e=(r2-r1)/(r1*cos(theta1)-r2*cos(theta1+dTheta));
%     e = abs(e);

    if e>1 %Check that the orbit is hyperbolic
        %2. Semieix
        a=(r1*(1+e*cos(theta1)))/(e^2-1);
        
        %3. Temps de viatge
        t2t1=(365.25/(2*pi))*a^(3/2)*(((e*sqrt(e^2-1)*sin(theta1+dTheta))...
            /(1+e*cos(theta1+dTheta)))-log((tan((theta1+dTheta)/2)...
            +sqrt((e+1)/(e-1)))/(tan((theta1+dTheta)/2)-sqrt((e+1)/(e-1))))...
            - ((e*sqrt(e^2-1)*sin(theta1))/(1+e*cos(theta1)))...
            + log((tan(theta1/2)+sqrt((e+1)/(e-1)))/...
            (tan(theta1/2)-sqrt((e+1)/(e-1)))));
        
        %4. Compute the error
        if imag(t2t1)==0
            error=abs(t2t1-dt);
            
            %avoid looping between 2 values
            if error==errorArray(1) && loopVal==maxLoopVal
                fprintf('\nlooping between 2 values of theta.\nSTOP LOOP');
                break;
            end
            if error==errorArray(1)
                loopVal=loopVal+1;
                dtheta=dtheta/1e1;
            end
            errorArray=[errorArray(2) error];
            
            %include the sign of the difference to perform a kind of
            %intelligent convergence
            theta1=theta1+dtheta*-((t2t1-dt)/abs(t2t1-dt));
            if theta1<0; theta1=theta1+2*pi; end
        else
            theta1=theta1+dtheta*10;
        end
    else
        error=1000;
        theta1=theta1+dtheta*10;
    end
end