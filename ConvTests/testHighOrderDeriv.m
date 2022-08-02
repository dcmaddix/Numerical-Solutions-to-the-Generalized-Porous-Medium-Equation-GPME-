close all
N = 128;
dx = 1 / N;
lb = 0;
ub = 1.0;
x = lb:dx:ub;
f = @(x) (x.^7 + x.^6 +x.^5);
p = f(x)';
i = 2:N;
%First Deriv
dp_dx_exact = (7*x.^6 + 6*x.^5 + 5*x.^4)';
%chnage so take all of ghostnodes
ghostNodes(1) = f(lb-dx);
ghostNodes(2) = f(ub+dx);
%odd for left bdry, even for right bdry needed for 4th order 3rd and 4th
%deriv
ghostNodes(3) = f(lb-2*dx);
ghostNodes(4) = f(ub+2*dx);
%Error goes down by 16 = 4^2 when decrease dx by 2 shows 4th order accurate
dp_dx = central_4thOrder(p, ghostNodes, dx);
%dp_dx = (p(i+1) - p(i-1)) / (2*dx);
figure(1)
plot(x(i),dp_dx_exact(i), x(i), dp_dx)
dx^(1/2)*norm(dp_dx_exact(i) - dp_dx)
legend('exact', 'numerical')
title('First Derivative 4th Order')

%Second Deriv
d2p_dx2_exact = (42*x.^5 + 30*x.^4 + 20*x.^3)';
%verified 4th order accurate error decreases by 16 when dx decreased by 2
d2p_dx2 = central_2ndDeriv_4thOrder(p, ghostNodes, dx);
%d2p_dx2 = (p(i+1)-2*p(i) + p(i-1)) / (dx^2);
figure(2)
plot(x(i),d2p_dx2_exact(i),x(i), d2p_dx2)
%max(abs((d2p_dx2_exact(i) - d2p_dx2))) / max(abs(d2p_dx2))
dx^(0.5) * norm(d2p_dx2_exact(i) -d2p_dx2)
legend('exact', 'numerical')
title('Second Derivative 4th Order')

%Third Derivative
d3p_dx3_exact = (210*x.^4 + 120*x.^3 + 60*x.^2)';
%verified 2nd order accurate error decreases by 4 when dx decreased by 2
%d3p_dx3 = central_3rdDeriv_2ndOrder(p, ghostNodes, dx);%2nd order works correct order of accuracy
%Error goes down by 16 = 4^2 when decrease dx by 2 shows 4th order accurate
d3p_dx3 = central_3rdDeriv_4thOrder(p, ghostNodes, dx);
figure(3)
plot(x(i),d3p_dx3_exact(i), x(i), d3p_dx3)
dx^(0.5) * norm(d3p_dx3_exact(i) - d3p_dx3)
%max(abs((d3p_dx3_exact(i) - d3p_dx3)))
legend('exact', 'numerical')
title('Third Derivative 4th Order')

 f = @(x) (1/pi^4)*(sin(pi*x));%(840*x.^3 +260*x.^2+120*x);% (1/pi^4)*(sin(pi*x)); %cos(x) -sin(x) -cos(x) sin(x)
 ghostNodes(1) = f(lb-dx);
 ghostNodes(2) = f(ub+dx);
% %odd for left bdry, even for right bdry needed for 4th order 3rd and 4th
% %deriv
 ghostNodes(3) = f(lb-2*dx);
 ghostNodes(4) = f(ub+2*dx);
 p = f(x)';
%Fourth Derivative
%d4p_dx4_exact = (840*x.^3 +360*x.^2+120*x)';
%d4p_dx4_exact = (120*x)';
d4p_dx4_exact = sin(pi*x)';
%verified 2nd order accurate eror decreases by 4 as dx decreased by 2
%4th order to 512 err 1e-6
d4p_dx4 = central_4thDeriv_2ndOrder(p, ghostNodes, dx);
%d4p_dx4 = (p(5:end) - 4*p(4:end-1) + 6*p(3:end-2) -4*p(2:end-3) +p(1:end-4))*dx^(-4);
%4th order on course grid then blows up dx^4 too big 4th order to 128
d4p_dx4 = central_4thDeriv_4thOrder(p, ghostNodes, dx); 
figure(4)
plot(x(i), d4p_dx4_exact(i), x(i), d4p_dx4)
dx^(0.5)*norm(d4p_dx4_exact(i) - d4p_dx4)
legend('exact', 'numerical')
title('Fourth Derivative 4th Order')
