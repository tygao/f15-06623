function tianyug1_prob1
% First linear regresion the data and then use fsolve to solve the question.
clc;
clear all;
close;
T = [2200, 2350, 2450, 2800, 3050, 3250, 3320, 3410]';
K1 = [0.043, 0.058, 0.060, 0.085, 0.106, 0.116, 0.135, 0.136]';
K2 = [0.003, 0.004, 0.005, 0.010, 0.013, 0.018, 0.020, 0.021]';
A1=ones(8,1);
A1(:,2)=-1/8.314./T;
b1=log(K1);
x1=inv(A1'*A1)*A1'*b1; % least-squares solution t
k01=exp(x1(1)); %calculate pre-exponential factor.
t=linspace(2200,3410); %prepare for ploting
y1=k01.*exp(-x1(2)/8.312./t); % y is least-squares fitting value.

A2=ones(8,1);
A2(:,2)=-1/8.314./T;
b2=log(K2);
x2=inv(A2'*A2)*A2'*b2; % least-squares solution t
k02=exp(x2(1)); %calculate pre-exponential factor.
y2=k02.*exp(-x2(2)/8.312./t); % y is least-squares fitting value.

figure
hold on;
plot(t,y1,t,y2);
plot(T,K1,'o',T,K2,'o');
legend('K1-T fit','K2-T fit', 'K1-T measured', 'K2-T measured')
title('Equilibrium constant versus temperature')
xlabel('Temperature, K');
ylabel('Equilibrimu constant')

% From the least-square Fit we can caculate K at T = [2500:200: 3500]
Tcal = 2500:200:3500;
K1cal = k01.*exp(-x1(2)/8.312./Tcal);
K2cal = k02.*exp(-x2(2)/8.312./Tcal);

guess = [0.1,0.1];
Ecal = ones(6,2);
for i = 1:6;

    F=@(e)[(4*e(1)^2)*((1/3)+e(1)-e(2))/((1/3)-2*e(1))^2-K1cal(i);
            (4*e(2)^2)/((1/3+e(1)-e(2))*(1/3-e(2)))-K2cal(i);];

    
    [ecal]=fsolve(F,guess);
    Ecal(i,:)=ecal;
end
sprintf('The reaction extent is ')
Ecal
end
