function nonlinfin()
% non-linear finit differences
N = 100; % total nodes
h = 1/N; % step size
D = 10;
L = 1000;
v = 0.1;
k = 5e-5;
guess = 0.1*ones(N,1); % guess for fsolve
y = fsolve(@prob,guess); % fsolve for nonlinear system
x = linspace(0,1,N);  % x for plot
plot(x,y);
xlabel('x');
ylabel('y');
title('Problem');


    function y = prob(x)
        % function of nonlinear finite difference system
        y = zeros(N,1);
        for i = 2:N-1;
            y(i) =(x(i+1)+x(i-1)-2*x(i))/h^2+v*L/D*(x(i+1)-x(i-1))/2/h-L^2*k*x(i)^2/D;
        end
        y(1) = 1-x(1);  % Dirichle BC
        y(N) = 0.1-x(N);  % Dirichle BC
        
    end
            
end