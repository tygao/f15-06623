function Assignment3
%Author: Tianyu Gao
%% Problem 1
% reaction extent x = -([S]-[S]0) = ([P] - [P]0)
% [S] = [S]0 - [P] + [P]0 = [S]0 - [P]
%g% $\frac{d[P]}{dt}=k_{cat}[E]_0\frac{[S]}{K_m+[S]}=k_{cat}[E]_0\frac{[S]_0-[P]}{K_m+[S]_0-[P]}$
S0 = 1; % mM
E0 = 0.1 * S0;

% a: Pepsin Catalyzed
kcat= 0.5; % s-1
Km = 0.3; % mM
P0 = 0; % initial value
tspan = [0 60]; % integration reange to show trends.
[T,p] = ode45(@rate,tspan ,P0);

s = S0 - p;
es = E0*s ./(Km+s);
e = E0*ones(numel(T),1);
figure
plot(T,p,T,s,T,es,T,e);
legend('[P]','[S]','[ES]','[E]');
xlabel('time[S]');
ylabel('Concentration [M]')
title('Pepsin catalyzed')

% b: Fumarase catalyzed 
% This becomes a stiff system.
kcat= 0.08; % s-1
Km = 5e-3; % mM
P0 = 0; % initial value
tspan = [0 180];
[T,p] = ode15s(@rate,tspan ,P0);

s = S0 - p;
es = E0*s ./(Km+s);
e = E0*ones(numel(T),1);
figure 
plot(T,p,T,s,T,es,T,e);
legend('[P]','[S]','[ES]','[E]');
xlabel('time[S]');
ylabel('Concentration [M]')
title('Fumarase catalyzed')

% Assume that the reaction between E, S and ES
% reach equalibrium instantly, due to the lack of rate constant.
function r=rate(t,P)
    S = S0 - P;
    r = kcat*E0*S/(Km+S);
end



%% Problem 2
% $\frac{dT}{dt}=\frac{2\sigma}{\rho d}\frac{T_F^4-T^4}{c_p(T)+T\frac{dc_p}{dT}}$
%
% $c_p(T) = 355.2 + 0.1004T$
%
% $\frac{dc_p}{dT} = 0.1004$
% 
% $\theta = T/T_F  \ \ \  T = \theta T_F$
%
% $\tau = t/t^*  \ \ \  t = \tau t^*$
%
% $\frac{d\theta}{d\tau}=\frac{t^*}{T_F}\frac{2\sigma}{\rho d}\frac{T_F^4(1-\theta)}{355.2+0.2008\theta T_F}$ 
%
rho = 8933; % kg/m^3
d = 0.002; % m
TF = 1200; % K
sigma = 5.676e-8; % W/m^2K^4
T0 = 300; % K
theta0 = T0/TF; % initial value
% Choose 300 s as t* for scaling
tstar = 300;  
tspan = [0,1]; % intergration interval 

myfun = @(x,y)tstar/TF*2*sigma/rho/d*TF^4*(1-y)/(355.2+0.2008*y*TF);
[t,y]= ode45(myfun,tspan,theta0);

plot(t,y);
xlabel('$\tau = t/t^*$','interpreter','latex')
ylabel('$\theta = T/T_F$','interpreter','latex')
title('$\theta (\tau)$ Plot','interpreter','latex')

%% Problem 3
% $\frac{dx_1}{dt} = -[x-cos(\theta_0)]+[x_2 - sin(\theta_0)]$
%
% $0 = x_1^2+x_2^2-1$
%
% $\left(\begin{array}{cc} 1 & 0\\ 0 & 0 \end{array}\right) \left(\begin{array}{c} \frac{dx_1}{dt}\\ \frac{dx_2}{dt} \end{array}\right)$
% $=\left(\begin{array}{c} -[x-cos(\theta_0)]+[x_2 - sin(\theta_0)]\\ x_1^2+x_2^2-1 \end{array}\right)$
%
M = [1 0;0 0];
theta0 = 0;
tspan = [0 1];
options = odeset('mass',M);
[t,x] = ode15s(@problem3,tspan,[0 1],options);
figure
subplot(1,2,1);
plot(x(:,1),x(:,2));
xlabel('x1');
ylabel('x2');
title('x2 verse x1 (initial value = [0,1]');
[t,x] = ode15s(@problem3,tspan,[0 0.8],options);
subplot(1,2,2);
plot(x(:,1),x(:,2));
xlabel('x1');
ylabel('x2');
title('x2 verse x1 (initial value = [0,0.8]');
% From the plot we can see matlab gives the same answer to different initial
% value, this shows the algebraic equation override the initial value.
function dxdt=problem3(t,x)
    dxdt = zeros(2,1);
    dxdt(1) = -(x(1)-cos(theta0)) + x(2)-sin(theta0);
    dxdt(2) = x(1)^2 + x(2)^2 -1;
end


%% Problem 4
% $\frac{dC_A}{dt}=-k_1C_A$
%
% $\frac{dC_B}{dt}=k_1C_A+k_3C_C-k_2C_B-k_4C_B$
%
% $\frac{dC_C}{dt}=k_2C_B-k_3C_C$
%
% $\frac{dC_D}{dt}=k_4C_B$
%
%
C0 = [10 0 0 0]'; % M initial value
k = [1 1 0.5 1]'; % 1/min rate constant
tspan = [0 20];
[t1,y]=ode45(@batch,tspan,C0);
figure
subplot(4,1,1);
plot(t1,y(:,1));
title('Concentration vs Time');
xlabel('Time [min]');

subplot(4,1,2);
plot(t1,y(:,2));
xlabel('Time [min]');

subplot(4,1,3);
plot(t1,y(:,3));
ylabel('C concentration [M]');
subplot(4,1,4);
plot(t1,y(:,4));
xlabel('Time [min]');
ylabel('D concentration [M]');
Cabatch = y(:,1);
function dxdt = batch(t,x)
    dxdt=zeros(4,1); % return a column vector
    dxdt(1) = -k(1)*x(1);
    dxdt(2) = k(1)*x(1)-k(2)*x(2)+k(3)*x(3)-k(4)*x(2);
    dxdt(3) = k(2)*x(2)- k(3)*x(3);
    dxdt(4) = k(4)*x(2);
end

%% Problem 5
% $\frac{dC_A}{dt}=QC_{A0}-QC_A-k_1C_AV$
%
% $\frac{dC_B}{dt}=QC_{B0}-QC_Bk_1C_AV+k_3C_CV-k_2C_BV-k_4C_BV$
%
% $\frac{dC_C}{dt}=QC_{C0}-QC_Ck_2C_BV-k_3C_CV$
%
% $\frac{dC_D}{dt}=QC_{D0}-QC_Dk_4C_BV$
Cfeed = [10 0 0 0]'; % M initial value
C0 = [10 0 0 0]'; % Initial concentration in CSTR
k = [1 1 0.5 1]'; % 1/min rate constant
V = 1; % L reactor volume
Q = 1; % L/min residence time 1 minute
[t2,y]=ode45(@CSTR,tspan,C0);
Cacstr = y(:,1);
subplot(4,1,1);
plot(t2,y(:,1));
xlabel('Time [min]');
ylabel('A [M]');
title('Concentration vs Time');
subplot(4,1,2);
plot(t2,y(:,2));
xlabel('Time [min]');
ylabel('B [M]');
subplot(4,1,3);
plot(t2,y(:,3));
xlabel('Time [min]');
ylabel('C [M]');
subplot(4,1,4);
plot(t2,y(:,4));
xlabel('Time [min]');
ylabel('D [M]');

figure
plot(t1,Cabatch,t2,Cacstr);
xlim([0,5]);
xlabel('Time [min]');
ylabel('Concentration [M]');
legend('Cabatch','Cacstr');
title('Batch vs CSTR');

function dxdt = CSTR(t,x)
    dxdt=zeros(4,1); % return a column vector
    dxdt(1) = Q*Cfeed(1)-Q*x(1)-k(1)*x(1)*V;
    dxdt(2) = Q*Cfeed(2)-Q*x(2)+k(1)*x(1)*V-k(2)*x(2)*V+k(3)*x(3)*V-k(4)*x(2)*V;
    dxdt(3) = Q*Cfeed(3)-Q*x(3)+k(2)*x(2)*V- k(3)*x(3)*V;
    dxdt(4) = Q*Cfeed(4)-Q*x(4)+k(4)*x(2)*V;
end


%% Problem 6
% $\frac{dL}{dt}=\frac{A_ih}{\rho_1c_{p,1}V_1}(C-L)$
%
% $\frac{dC}{dt}=\frac{A_oh}{\rho_2c_{p,2}V_2}(32-C)+\frac{A_ih}{\rho_2c_{p,2}V_2}(L-C)$

L0 = 150 ; % F initial condition;
C0 = 150; % F
%Parameters
rho1 = 62;
rho2 = 139;
Cp1 = 1.00;
Cp2 = 0.2;
V1 = 0.03;
V2 = 0.003;
Ai = 0.4;
Ao = 0.5;
h = 8.8;
tspan = [0 8];
[t,y]=ode45(@problem6,tspan,[L0,C0]);

plot(t,y(:,1),t,y(:,2));
xlabel('Time [hr]');
ylabel('Temperature [F]');
legend('Liquid','Container');
title('Temperature vs Time');
function dxdt =problem6(t,x)
    dxdt = zeros(2,1);
    dxdt(1) = Ai*h/rho1/Cp1/V1*(x(2)-x(1));
    dxdt(2) = Ao*h/rho2/Cp2/V2*(32-x(2))+Ai*h/rho2/Cp2/V2*(x(1)-x(2));
end


%% Problem 7
% $y''' - y''-2y = 2x^2+2x$
%
% define:
% $y'=u\ \ u'=v$
%
% $y' =u$
%
% $u'=v$
%
% $v'=2x^2+2x+2y-v$

% Problem a:
y0=[-1;0;-4];
tspan = [0 1];
[x,y]=ode45(@problem7, tspan,y0);
figure
plot(x,y);
legend('y','u','v');
title('Initial Value Problem');
% Problem b:
% Shooting methods for BVP
guess = 1;    % make a guess then solve until boundary value been satisfied.
i = fzero(@odefzero7,guess); % initial value of y(0)
sprintf('The initial value y(0) is %f',i)
yi=[i,0,-4];
[x,y]=ode45(@problem7,tspan,yi);
figure
plot(x,y);
legend('y','u','v');
title('Boundary Value Problem');

function r=odefzero7(i)
    yi=[i,0,-4];
    [x,y]=ode45(@problem7, tspan,yi);
    r=y(numel(x),1)+1; % boundry condition
end

function dydx=problem7(x,y)
    dydx=zeros(3,1);
    dydx(1) = y(2);
    dydx(2) = y(3);
    dydx(3) = 2*x^2+2*x+2*y(1)-y(3);
end

%% Problem 8
% $y''-y'+y=3e^(2x)-2sinx$
%
% define
% $y'=u$
%
% $y'=u$
%
% $u'=3e^{2x}-2sinx+u-y$
tspan = [1,2];
y1 = 6.308447;
y2 = 55.430436;
guess = 1; % guess of initial value
i = fzero(@odefzero8,guess);
yi = [y1;i];
[t,y] = ode45(@problem8,tspan,yi);
figure
plot(t,y(:,1));
xlabel('x');
ylabel('y');
title('y vesus x');
function r=odefzero8(i)
    yi = [y1;i];
    [t,y] = ode45(@problem8,tspan,yi);
    r = y(numel(t),1)-y2; 
end

function dydx=problem8(x,y)
    dydx=zeros(2,1);
    dydx(1)=y(2);
    dydx(2)=3*exp(2*x)-2*sin(x)+y(2)-y(1);
end


%% Problem 9
% $-u'' + \pi ^2 u = 2\pi ^2 sin(\pi x)
%
% Finite difference: Compared to $y''=p(x)y' + q(x)y + r(x)$
%
% we get $p(x)=0  \ q(x) = \pi ^2 \ r(x) = -2\pi ^2sin(\pi x)$
%
% define N = 10, then h = 1/10
%
% $d_i(1<i<N) = 2+\frac{\pi}{10}^2$
%
% $u_i = l_i =-1$
%
% $b =[\alpha, -h^2r_1, -h^2r_2, ... , -h^2r_{N-1}, \beta]^T$
%
 
N = 100;
h = 1/N;
alpha = 0;
beta = 0;
b = zeros(N-1,1);
xspan = linspace(0,1,N+1);
for i = 1:N-1
    b(i) = -h^2*(-2*pi^2)*sin(i*pi/N);
end
b = [alpha;b;beta];
A = diag((2+(pi/N)^2)*ones(1,N+1))+diag(-ones(1,N),1)+diag(-ones(1,N),-1);
A(1,1)=1;
A(N+1,N+1)=1;
A(1,2)=0;
A(N+1,N)=0;
x=A\b;
plot(xspan,x);
xlabel('x');
ylabel('u');
title('Problem 9');

% Shooting Methods:
u01 = 0;
uder01 = 0;
u02 = 0;
uder02 = 1;
u10 = [u01; uder01];
u20 = [u02; uder02];
tspan = [0; 0.25; 0.5; 0.75; 1];
[x1,u1]= ode45(@IVP1, tspan, [0 0]);
[x2,u2]= ode45(@IVP2, tspan, [0 1]);
u1(:,1)
u2(:,2)
% Neither u1(1) nor u2(1) equals to 0; we need another parameter c
c = (0-u1(5,1))/u2(5,1)
exact = sin(pi);
estim = u1(5,1)+c*u2(5,1);
error = estim - exact

% Initial Value Problem 1:
    function dudx = IVP1(x,u)
    dudx=zeros(2,1);
    dudx(1) = u(2);
    dudx(2) = pi^2*u(1)-2*pi^2*sin(pi*x);
    end
% Initial Value Problem 2:
    function dudx = IVP2(x,u)
    dudx=zeros(2,1);
    dudx(1) = u(2);
    dudx(2) = pi^2*u(1);  
    end

end