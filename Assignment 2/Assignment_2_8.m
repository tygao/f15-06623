function Assignment2_8
x0=[0.05,0.02]';
x =fsolve(@extenta,x0);
disp('Equilibrium Extent on Condition A is')
vpa(x,4)
x =fsolve(@extentb,x0);
disp('Equilibrium Extent on Condition B is')
vpa(x,4)
end

%% problem a:
function E=extenta(x) % vector x is the reaction extents of two reactions.
K =[0.1071;0.01493];  %Equilibrium constant
Ni=[1/3, 0, 1/3, 1/3, 0]';  % Feed mole vector
N(1) = Ni(1)-2*x(1);        % Equilibrium mole vector
N(2) = Ni(2)+2*x(1);
N(3) = Ni(3)+x(1)-x(2);
N(4) = Ni(4)-x(2);
N(5) = Ni(5)+2*x(2);
Nt = N(1)+N(2)+N(3)+N(4)+N(5);   % Totle mole
y(1) = N(1)/Nt;                  % mole fraction
y(2) = N(2)/Nt;
y(3) = N(3)/Nt;
y(4) = N(4)/Nt;
y(5) = N(5)/Nt;
E = [y(2)^2*y(3)/y(1)^2-K(1);
    y(5)^2/y(3)/y(4)-K(2)];
end

%% Problem b:
function E=extentb(x) % vector x is the reaction extents of two reactions.
K =[0.1071;0.01493];  %Equilibrium constant
Ni=[2, 0, 1/3, 1/3, 0]';  % Feed mole vector
N(1) = Ni(1)-2*x(1);        % Equilibrium mole vector
N(2) = Ni(2)+2*x(1);
N(3) = Ni(3)+x(1)-x(2);
N(4) = Ni(4)-x(2);
N(5) = Ni(5)+2*x(2);
Nt = N(1)+N(2)+N(3)+N(4)+N(5);   % Totle mole
y(1) = N(1)/Nt;                  % mole fraction
y(2) = N(2)/Nt;
y(3) = N(3)/Nt;
y(4) = N(4)/Nt;
y(5) = N(5)/Nt;
E = [y(2)^2*y(3)/y(1)^2-K(1);
    y(5)^2/y(3)/y(4)-K(2)];
end