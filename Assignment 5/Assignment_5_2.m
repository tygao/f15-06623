t = linspace(0,pi,10);
ydata = 2*sin(3*t)+2;

xfit = linspace(0,pi,1000);
y = 2*sin(3*xfit)+2;
yfit = spline(t,ydata,xfit);
figure
plot(t,ydata,'o',xfit,y,'b--',xfit,yfit)
title('Problem 2');
xlabel('x');
ylabel('y');
legend('Data Point','y=2*sin(3*t)+2','Spline Interpolation');