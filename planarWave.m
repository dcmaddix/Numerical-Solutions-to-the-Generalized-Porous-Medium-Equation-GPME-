%planar wave solution
N = 25;
kpow = 2;
dx = 1/ N;
dy = dx;
eps = 1e-9; %approach 0 to get true analytical solution
c = 1; %mimetic paper 0.4
%define domain bounds
lb_x = 0.0;
ub_x = 1.0; %3 in mimetic case
x = lb_x:dx:ub_x;
lb_y = 0.0;
ub_y = 1.0;
y = lb_y:dy:ub_y;
[x y] = meshgrid(x,y);
%search for indices less than ct
Nx = size(x,1);
Ny = size(x,2);
p = eps^(1/kpow)*ones(Nx,Ny); %initial condition is (eps)^(1/kpow) 1e-3
%compute at x_cent = 0.12
[coordx, coordy] = find(x == 0.12 & y == 0.08); %(.8,.64)
nt = 0.3; %t = 3;
dt = dx^2/32; %32 to compare with other cases
t = 0:dt:nt;
p_t = zeros(1,length(t)-1);
saveon = false;
nw = ones(2,1);
nw = nw / norm(nw);
k = p;
%can track same position for time plot: ADD
for n = 1:length(t);
    p_t(n) = p(coordx,coordy);
    k = kpow*c*(pos(c*t(n) - (nw(1)*x + nw(2)*y)));
    p = k.^(1/kpow);
    p(p == 0) = eps;
end
f = surf(x,y,p);
view([14.5 16]);
colorbar;
xlabel('x')
ylabel('y')
zlabel('p')
title('2D Planar Wave solution on x-axis')
if (saveon)
    saveas(f, sprintf('~/Documents/Stanford/Research/2Dres/planar/exactsolN%i_t%i.jpg', N,nt));
end
f = figure(2);
plot(t,p_t)
ylabel('p_t')
xlabel('t')
title('Pressure at (x,y) vs time')
if (saveon)
    saveas(f, sprintf('~/Documents/Stanford/Research/2Dres/planar/exactptN%i_t%i.jpg', N,nt));
end
save(sprintf('data/2D/m%i/planar/exact%i.mat',kpow,N),'x','y','t', 'p_t','p')