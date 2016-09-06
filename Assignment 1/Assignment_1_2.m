%{
Written by: Tianyu Gao
Born on: Sept 9, 2015

" A linear system of equations Ax=b and its martix A whose condition number is small are well-conditioned.
A large condition number indicates ill-conditiong"
%}
clear all
clc
close
figure
j=1;
hold on
for i=-100:0.5:100
    A=[-3 -2 1; 2 i 1; 4 1 -2];
    x(j)=i;
    y(j)=cond(A);
    j=j+1;
end
plot(x,y)
title('condition number versus alpha (step length=0.5)')
xlabel('alpha');
ylabel('condition number');
[max,pos]=max(y);
x_max=x(pos);
sprintf('When alpha equas %g, this system is ill-conditioned', x_max)
% due to the extreme large maxmium, the other points are hard to see, 
% so change the domain and loop step size.
j=1;
for i=-50:1:50
    A=[-3 -2 1; 2 i 1; 4 1 -2];
    x1(j)=i;
    y1(j)=cond(A);
    j=j+1;
end
figure
plot(x1,y1)
title('condition number versus alpha (step length=1)');
xlabel('alpha');
ylabel('condition number');

