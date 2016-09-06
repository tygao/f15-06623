function Assignment_5_1
%% Problem 1a: Vandermonde matrix 
%
% $\left(\begin{array}{ccccc} 1 & x_1 & x_1^2 & x_1^3 & x_1^4\\ 1 & x_2 & x_2^2 & x_2^3 & x_2^4 \\ 1 & x_3 & x_3^2 & x_3^3 & x_3^4 \\1 & x_4 & x_4^2 & x_4^3 & x_4^4\\1 & x_5 & x_5^2 & x_5^3 & x_5^4\end{array}\right)$

x = [0.1, 0.5, 1, 1.5, 2]';
y = log(x);
A = [x.^0, x, x.^2, x.^3,x.^4];
b = A\y;
%sprintf('The coefficients of polynomial is')
%b
tspan=linspace(0,2);
sola=[];
for i=1:numel(tspan);    
    p = [tspan(i).^0, tspan(i), tspan(i).^2, tspan(i).^3, tspan(i).^4];
    q = p*b;
    sola = [sola,q]; 
end

%% Problem 1b: Lagrange 
%
% $f_5(x)=\sum\limits_{i=1}^{5}L_i(x)f(x_i)$
%
% $L_i(x) = \prod\limits_{j=1,j\neq i}^{5}\frac{x-x_j}{x_i-x_j}$
A = zeros(1,5);
for i = 1:5;
    A(i) = y(i);
    for j = 1:5;
        if j~=i;
            A(i) = A(i)/(x(i)-x(j));
        end
    end
end


%{    
A(1) = y(1)/((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))*(x(1)-x(5)));
A(2) = y(2)/((x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))*(x(2)-x(5)));
A(3) = y(3)/((x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))*(x(3)-x(5)));
A(4) = y(4)/((x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3))*(x(4)-x(5)));
A(5) = y(5)/((x(5)-x(1))*(x(5)-x(2))*(x(5)-x(3))*(x(5)-x(4)));
%}
sprintf('The coefficients of Lagrange polynomial is')
A

function y = lagrange(t)
    B = ones(5,1);
    for i = 1:5;
        for j = 1:5;
            if j~=i;
                B(i)=B(i)*(t-x(j));
            end
        end
    end
    y = A*B;
end
solb=[];
for i=1:100;    
    solb = [solb,lagrange(tspan(i))]; 
end


%% Problem 1c: Newton
% $f[x_i,x_j] = \frac{f(x_i)-f(x_j)}{x_i-x_j}$
%
% $f[x_i,x_j,x_k] = \frac{f[x_i,x_j]-f[x_j,x_k]}{x_i-x_k}$
%
% $...$
A = zeros(1,4);
for i=1:4
    A(i) = (y(i)-y(i+1))/(x(i)-x(i+1));
end

B = zeros(1,3);
for i = 1:3
    B(i) = (A(i)-A(i+1))/(x(i)-x(i+2));
end

C = zeros(1,2);
for i = 1:2;
    C(i) = -(B(i)-B(i+1))/(x(i)-x(i+3));
end
D = (C(1)-C(2))/(x(1)-x(5));
sprintf('The coefficients of  Newton polynomial is')
E =[y(1),A(1),B(1),C(1),D]


solc=[];
for i=1:100;    
    solc = [solc,newton(tspan(i))]; 
end

function y = newton(t)
    y = E(1)+E(2)*(t-x(1))+E(3)*(t-x(1))*(t-x(2))-E(4)*(t-x(1))*(t-x(2))*(t-x(3))-E(5)*(t-x(1))*(t-x(2))*(t-x(3))*(t-x(4));
end

%% Problem 1d: Plot 
figure;
plot(x,y,'bo',tspan,sola,'r',tspan,solb,'-.',tspan,solc,'--',tspan,log(tspan),'b')
legend('Data Points','Prob a','Prob b','Prob c','y = ln(x)');
title('Problem 1');
% There no obvious difference between three polynomial interpolation
% methods.

end