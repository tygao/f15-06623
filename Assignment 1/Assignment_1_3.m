%{
Written by: Tianyu Gao
Born on: Sept 9, 2015
%}
clear all
clc
close
%% Problem a)
k1=3/60;
k2=0.1;
c0=[3 0 0]';
figure
hold on
for t= 1:1:200
    A=[1+k1*t 0 0; -k1*t 1+k2*t 0; 0 -k2*t 1]; % A is the coefficients matrix
    c=A\c0;

    m(t,:)=t; % m to record t
    n(:,t)=c; % n to record c
    
end
plot(m,n(3,:),'-');
title('Concentration of c vesus the residence time')
xlabel('Residence time t min ');
ylabel('Concentration of c mol/L');
%% Problem b)
figure;
hold on
t=0.5;
k1=0.5;
i=1;
for k2= 0:0.1:50;
    A=[1+k1*t 0 0; -k1*t 1+k2*t 0; 0 -k2*t 1];
    c=A\c0;
    p(i)=k2;
    q(i)=c(3,:);
    i=i+1;
end
plot(p,q);
title('Concentration of c versus k2');
xlabel('Rate constant k2 min^-1');
ylabel('Concentration of c mol/L');