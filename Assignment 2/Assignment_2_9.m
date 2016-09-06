function Assignment_2_9
clear all
x0=[0.1 0.1 0.1 0.1 0.1 0.1 0.1]'; %initial guess
x =fsolve(@myfun,x0);
vpa(x,4)
end

function F=myfun(x)
l=[100,100,200,75,100,75,50];
F=[x(1)-x(2)-x(6); % q1 = q2+q6
x(1)-x(7);         % q1 = q7
x(2)-x(3)-x(4);    % q2 = q3 + q4
x(2)-x(5);         % q5 = q3 + q4 = q2
l(3)*x(3)^2-l(4)*x(4)^2; 
l(2)*x(2)^2+l(4)*x(4)^2+l(5)*x(5)^2-l(6)*x(6)^2;
l(1)*x(1)^2+l(6)*x(6)^2+l(7)*x(7)^2-5.2*10^5*pi^2*0.2^5/8/0.02/998];
end