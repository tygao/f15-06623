%{
Written by Tianyu Gao
Born on Sept, 21, 2015
Power method for eigenvalues
%}
clc
clear all

A=[2,-1,0;-1,2,-1;0,-1,3];
x0 = [1, 0, 0]';
N = 1000; % maximum iteration times.
eps = 10^-6; % error
i = 1;
x = A * x0;
x = x/norm(x,Inf);
s = x0;

while norm(s-x)>=eps && i<N % if p-q less than error, we find the eigenvalue
    %and if i more than maximun times, stop iteration.
    s = x;
    y= A*x; 
    x=y/norm(y,Inf);% normalize vector x to converge;
    i = i + 1;
end
if i == N
    sprintf('Fail: reaching maximum iteration times.')
else
    value = norm(y,Inf)
    sprintf('The max eigenvalue is %f',(value))
    disp('The eigenvector is')
    x/norm(x)
end