function opt_nonlinfit()
xdata = [100 150 338 1140 2563 3844 8650];
ydata = [5.23 5.09 3.75 2.18 0.64 0.23 0.01];
% data provieded

guess =[2,2e-3,1e-3]; % initial guess for optimization
par = fminunc(@func2,guess)
plot(xdata,ydata,'o',xdata,func1(par));
legend('Data','Fit');
xlabel('x');
ylabel('y');
title('Problem Optimization: Non-linear Fit of Data');

    function y = func1(x)  % evalutaion function
        y = x(1)./(x(2)+exp(x(3).*xdata));
    end
    function d = func2(x)  % function of squared difference
        yfit = func1(x);
        e = yfit-ydata;
        d = 1/2*dot(e,e);
    end
end
