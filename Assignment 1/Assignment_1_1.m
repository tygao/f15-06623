%{
Written by: Tianyu Gao
Born on: Sept 6, 2015
%}
clear all
clc
close
x=(1:2:20); % assign x as a vector
Z=zeros(10,3); % assign Z as a 10-by-3 matrix
z1=50*rand(1,10)'; %z1 as first column of Z
z2=z1.^2; %z2 as second column of Z
z3=z1.^(1/2);  %z3 as third column of Z
Z=[z1,z2,z3];  %assign Z
plot(x,z2,'o');
grid;
xlabel('x');
ylabel('z2');
title('x versus the second column of Z');