function comparisonSlopePlot(dx,t_n)
x = -1:dx:1;
uj = JST(dx,0.5,3);
uexact = zeros(size(uj));
uexact( x <= 0) = x(x <=0) / (t_n -1);
uexact(x > 0) = 0;
figure(1);
plot(x,uexact,x,uj);
legend('uexact','JST')
xlabel('x')
ylabel('u(x,t)')
title('Comparison Plot of Exact Solution versus JST Solution at time t < 1')

figure(2)
u_min = SLIP(dx,1,3,'minmod');
u_van = SLIP(dx,1,3,'vanleer');
u_super = SLIP(dx,1,3,'superbee');
plot(x,uexact,x,u_min,x,u_van,x,u_super);
legend('uexact','minmod', 'vanleer', 'superbee')
xlabel('x')
ylabel('u(x,t)')
title('Comparison Plot of Exact Solution versus SLIP with various Limiters Solutions at time t < 1')

figure(3)
u_min = USLIP(dx,1,3,'minmod');
u_van = USLIP(dx,1,3,'vanleer');
u_super = USLIP(dx,1,3,'superbee');
plot(x,uexact,x,u_min,x,u_van,x,u_super);
legend('uexact','minmod', 'vanleer', 'superbee')
xlabel('x')
ylabel('u(x,t)')
title('Comparison Plot of Exact Solution versus USLIP with various Limiters Solutions at time t < 1')
end
