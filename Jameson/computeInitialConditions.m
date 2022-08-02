function u = computeInitialConditions(IC,x)
if (IC == 1)
    u(x <= 0) = 1;
    u(x > 0 ) = -1;
elseif (IC == 2)
    u = sin(pi*x);
else
    u(x <= 0) = -x( x<= 0);
    u(x> 0) = 0;
end
figure(1)
plot(x,u)
axis([-2 2 -1 1])
title('Initial Condition at time t = 0')
xlabel('x')
ylabel('u(x,0)')
end