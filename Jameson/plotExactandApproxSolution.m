function plotExactandApproxSolution(u,x,method, limiter, t_n, IC)
uexact = zeros(size(u));
figure(3)
str = sprintf('Comparison of Approximate %s Solution and Exact Solution to Burgers Equation at time T = %.2f',...
    method,t_n);
if (IC == 1)
    uexact(x <= 0)=1;
    uexact(x > 0)= -1 ;
    plot(x,uexact,'-',x,u,'o','MarkerSize',4)
    axis([-1 1 -2 2])
    legend('Exact Solution', limiter)
    title(str)
elseif (IC == 3)
    if t_n < 1
        uexact( x <= 0) = x(x <=0) / (t_n -1);
        uexact(x > 0) = 0;
    else
        uexact(x <= 0 ) = 1;
        uexact(x > 1) = 0;
    end
    plot(x,uexact,'-',x,u,'-ro','MarkerSize',4)
    hold on;
    plot(x,u,'-r')
    axis([-1 1 -0.5 1.5])
    legend('Exact Solution', limiter)
    title(str)
else
    %uexact periodic 0.7314x/t
    plot(x,u,'-o','MarkerSize',4)
    hold on;
    plot(x,u,'-b')
    legend(limiter)
    str = sprintf('Approximate %s Solution to Burgers Equation at time T = %.2f', method, t_n);
    title(str)
    axis([-2 2 -1 1])
end
xlabel('x')
ylabel('u(x,T)')
end