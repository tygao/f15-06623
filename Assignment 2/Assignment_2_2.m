%{
Written by Tianyu Gao
Born on Sept, 20
%}
clc
clear all
%% Problem a
A = [3 2 2 1; 2 3 1 2; -1 1 2 0; 2 4 3 5];
[V,D] = eig(A);
V
D
sprintf('Norms of these  four eignevectors are %.2f, %.2f, %.2f and %.2f ', ...
 (norm(V(:,1))),(norm(V(:,2))),(norm(V(:,3))),(norm(V(:,4))))
% So we can see that Matlab uses 2-norm to normalize eigenvectors.

%% Problem b
A = sym(A);
[Ve,De]=eig(A);
Ve
De
% When use sym function. Matlab will output the exact result.