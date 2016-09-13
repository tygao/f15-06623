function nonlin_shooting()
guess = 1;  
phi2 = 0.01;
t =fzero(@residue,guess);
sprintf('The initial value y(0) is %f',t)
x0 =[t;0]; %initial value
tspan =[0,1];
[x,y] = ode45(@prob,tspan,x0);
plot(x,y(:,1));
xlabel('x');
ylabel('y');
title('Problem')
    % define function of first order differential equation system
    function dydx = prob(x,y)
        dydx = zeros(2,1);
        dydx(1) = y(2);
        if x==0
            dydx(2) = phi2*y(1);
        else
            dydx(2) = -2/x*y(2)+phi2*y(1);
        end
    end
    % define function of residue 
    function r = residue(x)
        x0 = [x;0]; % initial condition with guess
        tspan=[0,1];
        [t,y] = ode45(@prob,tspan,x0);
        cal = y(end,1);  % calculate value for bondary condition
        r = cal-1;     % boundary residue
    end
        
end