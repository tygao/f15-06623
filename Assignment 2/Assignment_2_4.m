%{
Written by Tianyu Gao
%}
function Assignmeng_2_4()
%% Problem a:
% $f(x) = x^3 - 5x^2 + 7x - 3$
%
% $f(x) = (x - 1)^2(x-3)$
%
% The roots of $f(x)$ is 1,1,3
%% Problem b,c:
% $f(x) = x^3 - 5x^2 + 7x - 3$
%
% $f^{(1)}(x) = 3x^2 - 10x$
%
% $f^{(2)}(x) = 6x$
clc
clear all
x1 = 0; % initial guess for standard update function
x2 = 0; % initial guess for modified update function
X1 = zeros(5,1);
X2 = zeros(5,1);

for i = 1:5
    y1 = NR(x1);
    x1 = x1 + y1;
    X1(i) = x1;
    
end
for i = 1:5
    y2 = mNR(x2);
    x2 = x2 + y2;
    X2(i) = x2;
end
disp('Problem b: Standard function:')
vpa(X1,8)
disp('Problem C: Modified function:')
vpa(X2,8)



%% Proble d:
x1 = 4; % initial guess for standard update function
x2 = 4; % initial guess for modified update function
X1 = zeros(5,1);
X2 = zeros(5,1);

for i = 1:5
    y1 = NR(x1);
    x1 = x1 + y1;
    X1(i) = x1;
    
end
for i = 1:5
    y2 = mNR(x2);
    x2 = x2 + y2;
    X2(i) = x2;
end
disp('Problem d: Standard function:')
vpa(X1,8)
disp('Problem d: Modified function:')
vpa(X2,8)
end

%% Function
 function [y] = fx(x)    % f(x)
 y = x^3 - 5*x^2 + 7*x -3;
 end
 function [y] = f1x(x)   % first derivative function 
 y = 3*x^2 - 10 * x + 7;
 end
 function [y] = f2x(x)  % second deritivative function 
 y = 6*x - 10;
 end
 function [y] = NR(x)   % standard update function
 y = -fx(x)/f1x(x);
 end
 function [y]= mNR(x)   % modified update function
 y = -fx(x)*f1x(x)/(f1x(x)^2 - fx(x)*f2x(x));
 end
 % The result is different from that of b,c, which means different initial
 % guesses may yeild different roots.