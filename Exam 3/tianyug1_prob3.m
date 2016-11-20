function tianyug1_prob3

% This is a IVP problem;
D = 0.02;
A1 = 2;
A2 = 1;
A3 = 0.5;
g = 9.8;
Q0 = 1/60; %m3/s
tspan = [0 3*10^5];
h0 = [0,0,0];
[t,y] = ode45(@prob,tspan,h0);
plot(t,y);
legend('h1','h2','h3');
xlabel('Time (s)')
ylabel('Hight (m)');
title('Problem 3');

    function dhdt = prob(t,h)
        dhdt=zeros(3,1);
        Q1=(2*g*h(1))^0.5*pi*D^2/4;
        Q2=(2*g*h(2))^0.5*pi*D^2/4;
        Q3=(2*g*h(3))^0.5*pi*D^2/4;
        dhdt(1) = 1/A1*(Q0-Q1);
        dhdt(2) = 1/A2*(Q1-Q2);
        dhdt(3) = 1/A3*(Q2-Q3);
    end
% Time is large, Might be something wrong with unit scale 
end
