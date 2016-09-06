function Assignment_5_4
%% a Linear Interpolation
T = 125:50:325;
P = [1.802, 9.25, 27.51, 63.45, 130.2];
Tint = [150, 300];
Tspan = linspace(125,325);
Pa = interp1(T,P,Tint,'linear')
Paplot = interp1(T,P,Tspan,'linear');

%% b Cubic Spline Interpolation
Pb = spline(T,P,Tint)
Pbplot = spline(T,P,Tspan);

%% c Nonlinear Fit
guess =[1,1,1]; % initial guess for optimization
para = fminunc(@func4b,guess);
func4c(Tint)
figure
plot(T,P,'o',Tspan,func4c(Tspan),'--',Tspan,Paplot,Tspan,Pbplot);
legend('Data','Nonlinear Fit','Linear Interpolation','Cubic Spline Interpolation','location','northwest');
xlabel('x');
ylabel('y');
title('Problem 4');

    function y = func4a(x)  % evalutaion function
        y = exp(x(1)+x(2)./(T+x(3)));
    end
    function d = func4b(x)  % function of squared difference
        yfit = func4a(x);
        e = yfit-P;
        d = 1/2*dot(e,e);
    end
    function y = func4c(x)  % evalutaion function
        y = exp(para(1)+para(2)./(x+para(3)));
    end
end
