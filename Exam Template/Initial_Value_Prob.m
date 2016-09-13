function IVP()
Cf = 3;
Q = 60e-6;
H = -2e8;
rho = 1000;
k0 = 4.5e6;
Tf = 298;
V = 18e-3;
Cp= 4000;
E = 7550;
k = k0*exp(-E/Tf);
tau = V/Q;
alpha = k*tau;
beta = H*Cf/rho/Cp/Tf;

y0 = [0,1.5]; % initial value
tspan = [0 5]; % integrate interval
[t,y]=ode45(@prob,tspan,y0);
plot(t,y);
xlabel('x');
ylabel('y');
title('Problem');
    function dydx = prob(~,y) % differential equations
        dydx = zeros(2,1);
        k = k0*exp(-E/Tf)*exp(1/y(2));
        alpha = k*tau;
        dydx(1) = 1-(alpha+1)*y(1);
        dydx(2) = 1-y(2)-alpha*beta*y(1);
    end

end