%{
Written by Tianyu Gao
Born on Sept 16, 2015
%}

%%
% $\frac{dVC_A}{dt}=vC_{A,0}-vC_{A}-k_1C_AV$
%
% $\frac{dVC_B}{dt}=vC_{B,0}-vC_{B}+k_1C_AV-k_2C_BV+k_3C_DV-k_4C_BV$
% 
% $\frac{dVC_C}{dt}=vC_{C,0}-vC_{C}+k_4C_BV$
%
% $\frac{dVC_D}{dt}=vC_{D,0}-vC_{D}+k_2C_BV-k_3C_DV$
%
% Steady State and $\tau=\frac{V}{v}$
%
% $C_{A,0}=C_{A}+k_1C_A\tau$
%
% $C_{B,0}=C_{B}-k_1C_A\tau+k_2C_B\tau-k_3C_D\tau+k_4C_B\tau$
%
% $C_{C,0}=C_{C}-k_4C_B\tau$
%
% $C_{D,0}=C_{D}-k_2C_B\tau+k_3C_D\tau$
% 
% 

clc
clear all
close all

k = [0.1, 0.2, 0.1, 0.8]; % /sec
t = 10/1; % sec
Ci = [5; 0; 0; 1]; %feed 
A = [1+k(1)*t, 0, 0, 0; -k(1)*t, 1+k(2)*t+k(4)*t, 0,-k(3)*t; 0, -k(4)*t, 1, 0; 0, -k(2)*t, 0, 1+k(3)*t,];

[L, U, P]=lu(A); % P is permutation matrix, which L * U = P * A
disp('Solution is')
Css = U\(L\(P*Ci))
disp('Triangular matrix mutiple the feed vector is')
U*Css 