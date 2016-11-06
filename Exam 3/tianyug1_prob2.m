function tianyug1_prob2
% This is a BVP problem
% First scale the equation with r = R *r';
% get y''+1/x*y' (h*A*R^2)/(1-e)/k * y= 0 
R = 1.3;
h = 0.001;
e = 0.36;
A = 15;
k = 0.0034;

alpha = (h*A*R^2)/(1-e)/k;

% Choose shooting method to solve this BVP
guess = 4; %guess for b.C
t =fsolve(@residue,guess);
%sprintf('The initial value y''(1) is %f',t)
x0 =[t;0]; %initial value
tspan =[0 1];
[x,y] = ode45(@prob,tspan,x0);

plot(x,y(:,1));
xlabel('Scaled radius = r/R');
ylabel('Scaled temperature');
title('Problem 2')
    % define function of first order differential equation systems
    function dydx = prob(x,y)
        dydx = zeros(2,1);
        dydx(1) = y(2);
        if x ==0;
            dydx(2) = 0;
        else
            dydx(2) = alpha*y(1)-1./x*y(2);
        end
        
    end
    % define function of residue 
    function r = residue(x)
        x0 = [x;0]; % initial condition with guess
        tspan=[0 1];
        [t,y] = ode45(@prob,tspan,x0);
        cal = y(end,1);  % calculate value for bondary condition
        r = cal-1;     % boundary condition residue
    end
        

end
