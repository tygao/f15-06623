%{
room temperature 25 centigrade, H2O  coefficient of viscosity u= 0.8949e-3 Pas
b =150 microns = 1.5e-4 m
delta P/L = 1atm / 3 inch = 101325 Pa/ (3*0.0254) m = 1.3297e+6 Pa/m;
%}
%% a), 4 nodes
y=linspace(0,1.5e-4,6)';
deltay= 1.5e-4/5;
G=-1.3297e+6/0.8949e-3*deltay^2.*ones(4,1)
A=[2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2] % coefficients A is a sparse and banded matrix
v=A\G;
v=[0;v;0]; % add v(0)=0, v(B)=0 as the boundry condition.
figure
plot(v,y); % negetive velocity indicate the direction 
title('Velocity profile (4 nodes)')

%% b, 15 nodes;
n=15;
y=linspace(0,1.5e-4,n+2);
deltay=1.5e-4/(n+1);
G=-1.3297e+6/0.8949e-3*deltay^2.*ones(n,1);
A= zeros(n,n);
for i=2:1:n
    A(i,i)=2;
    A(i-1,i)=-1;
    A(i,i-1)=-1;
end
A(1,1)=2;
A(n,n)=2;
v=A\G;
v=[0;v;0];
figure
plot(v,y);
title('Velocity profile (15 nodes)')









