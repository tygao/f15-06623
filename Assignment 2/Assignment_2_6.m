function Assignment_2_6()
clear all
x=[0;0];
i=0;
N = 10000;
eps = 10^-6;
J = Jac(x);
F = Fun(x);
while norm(F,Inf)>eps && i<N
J = Jac(x);
F = Fun(x);
detla = J\(-1*F);
x = x + detla;
i = i+1;
end

if i >= N
    disp('Fail: reaching maximum iteration times')
else
    vpa(x,8)
end
figure ;
hold on;
f1=ezplot('y-(x-1).^2',[-10 10]);
set(f1,'color','r');
f2=ezplot('(y+4).^2 - tan(x)',[-10 10]);
set(f2,'color','b');
plot(x(1),x(2),'*y');
legend('y-(x-1)^2=0','(y+4)^2-tanx=0','Solution Point')
title('y versus x');
hold off;
end

function J = Jac(x) 
J = [-2*x(1)+2, 1;
    -sec(x(1)).^2, 2*x(2)+8];
end
function F = Fun(x)
F = [x(2)-(x(1)-1).^2;
    (x(2)+4).^2 - tan(x(1))];
end

