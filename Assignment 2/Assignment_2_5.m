function Assignment_2_5
m=1;  % initial guess 1
n=2;  % initial guess 2
N = 1000;  % maxiumn iteration times
eps = 10^-6; % tolerance

for i = 1:1:N
    p = myfun(m); 
    q = myfun(n);
    delta = -q*(n-m)/(q-p); % I use Secant methods.
    m = n;
    n = n + delta;
    if abs(n-m)<eps  % stop iteration when |x(k)-x(k-1)|<eps
        break;
    end
end
if i == N
    disp('Fail: reaching maxiumn iteration times');
else
sprintf('When h = %g, velocity requirements is satisfied.',n)
end
end
    


function F = myfun(x)
V = 5; % m/s
t = 3; % s
L = 5; % m
g = 9.81; %m/s^2
F = V/(2*g*x)^0.5-tanh(t*(2*g*x)^0.5/2/L);
end
