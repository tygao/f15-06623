function Assignment_5_3
%% Problem 3a: Spline Interpolation
M=[0.5, 1.0, 1.5, 2.0, 2.5];
mu=[0.135, 0.165, 0.195, 0.207, 0.223];
xint = 1.25;
ya = spline(M,mu,xint)
tspan = linspace(min(M),max(M));
yaplot = spline(M,mu,tspan);

%% Problem 3b: Linear Regression
A = [ones(5,1),(1./M)'];
b = (1./mu)';
x = inv(A'*A)*A'*b;

mu_0 = 1/x(1)
alpha = x(2)*mu_0

yb = mu_0*xint/(xint+alpha)
ybplot = mu_0.*tspan./(tspan+alpha);

%% Problem 3c: Nonlinear Fit
guess =[0.2,0.3]; % initial guess for optimization
par = fminunc(@func3b,guess)

yc = func3c(xint)


%% Problem 3d: Plot
figure
plot(M,mu,'o',tspan,yaplot,tspan,ybplot,tspan,func3c(tspan),'--');
legend('Data','Spline Interpolation','Linear Fit','Nonlinear Fit','location','northwest');
xlabel('M');
ylabel('\mu');
title('Problem 3');

    function y = func3a(x)  % evalutaion function
        y = x(1).*M./(M+x(2));
    end
    function d = func3b(x)  % function of squared difference
        yfit = func3a(x);
        e = yfit-mu;
        d = 1/2*dot(e,e);
    end
    function y = func3c(x)
        y = par(1).*x./(x+par(2));
    end
end