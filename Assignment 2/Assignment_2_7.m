function Assignment_2_7()
%clear all
P = 760; %mmHg
T1 = fzero(@(x) exp(-3848.09/(x+273.15)+17.5318)-P,10); % boiling point of alpha
T2 = fzero(@(x) exp(-4328.12/(x+273.15)+17.913)-P,10);  % boiling point of beta
T = linspace(T1,T2);
x = zeros(100,1); % assign vector x;
y = zeros(100,1); % assign vector y
for i = 1:1:100   % iteration times equalls 100, decided by 'linspace'.
    P1 = P1(T(i)); % function P1 to calculate pressure of alpha at temperature T(i)
    P2 = P2(T(i)); % function P2 to calculate pressure of beta at temperature T(i)
    p = fzero(@(x)P-P1*x-P2*(1-x),0.5); % solve for liquid composition
    q = fzero(@(y)1/P-y/P1-(1-y)/P2,0.5); % solve for gas composition
    x(i)=p; 
    y(i)=q;
end
figure;
hold on;
plot(x,T,y,T);
axis([0,1,T1-5,T2+5]);
xlabel('Mole fraction of A');
ylabel('Temperture [$\circ$C]','Interpreter','LaTex');
Tb = fzero(@(x) 0.5*exp(-3848.09/(x+273.15)+17.5318)+0.5*exp(-4328.12/(x+273.15)+17.913)-P,90);
% solve for bubble point for equimolar mixture. 
Td = fzero(@(x) 0.5/exp(-3848.09/(x+273.15)+17.5318)+0.5/exp(-4328.12/(x+273.15)+17.913)-1/P,100);
% solve for dew point for equimolar mixture. 
plot(0.5,Tb,'o')
plot(0.5,Td,'o')
legend('bubble','dew','buble point solution','dew point solution')
plot(0,Tb,'*')
plot(0,Td,'*')
text(0.02,Tb,'bubble temperture')
text(0.02,Td,'dew temperture')
end

function P=P1(T)
A = -3848.09;
B = 17.5318;
P = exp(A/(T+273.15)+B);
end
function P=P2(T)
A = -4328.12;
B = 17.913;
P = exp(A/(T+273.15)+B);
end
