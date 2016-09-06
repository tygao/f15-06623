
%{
k = k0*exp(-Ea/(RT))
lnk0 - 1/(RT)*Ea = lnk
hence, 
lnk0 - 1/(RT1) * Ea = lnk1;
lnk0 - 1/(RT2) * Ea = lnk2;
...
lnk0 - 1/(RT5) * Ea = lnk5;

Assume, x1 = lnk0, x2= Ea; b=[lnk1; lnk2;...;lnk5];
A=[1, -1/(RT1); 1, -1/(RT2);... ; 1, -1/(RT5);

Ax=b, linear system.
%}
clc;
clear all;
close;
T=[430, 450, 470, 480, 490]';
k=[0.0026, 0.0118, 0.0460, 0.0873, 0.1800]';
A=ones(5,1);
A(:,2)=-1/8.314./T;
b=log(k);
x=inv(A'*A)*A'*b; % least-squares solution to the normal equations: 
%from the book "Modeling and Analysis Principles for Chemical & Biological
%Engineers".
k0=exp(x(1)); %calculate pre-exponential factor.
t=linspace(400,500); %prepare for ploting
y=k0.*exp(-x(2)/8.312./t); % y is least-squares fitting value.
figure
hold on;
plot(t,y);
plot(T,k,'o');
legend('k-T fit', 'k-T measured')
title('Rate constant versus temperature')
xlabel('Temperature, K');
ylabel('Rate constant, s^-1')
sprintf('Pre-exponential factor k0= %g', k0)
sprintf('Activation Energy Ea= %g J/mol', x(2))




