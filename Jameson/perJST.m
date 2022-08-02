function perJST(dx,t_n, IC) %this function solves the 1D Burger's equation 
%u_t + uu_x = 0 on [-1,1] using the JST scheme
%u-initial condition at time t = 0
%dx-spatial timestep use with CFL condition to compute timestep
%t_n is the final timestep that we evolve to
eps = 10^-3;
K = 1/2;
q = 3;
%x = -1:dx:1;
x=-2:dx:2; %spatial vector
dt = dx^2;
lambda = dt/dx;
t = dt:dt:t_n;
%initial burgers condition
u = computeInitialConditions(IC,x);
n = length(x); %we need 2 ghost cells second order method on each side
f=@(v)1/2*v.^2;%make flux function for Burgers
r = 3/2;
R = @(u,v) abs( (u-v) / (max( (abs(u) + abs(v)), eps*dx^r))) ^ q;
%if u and v have opposite signs like at shocks or extremum R will have value 1
unew = u;
k = @(j) (mod(j-1,n) + 1);
for niter = 1:length(t) %time integration
    %for j = 3:n-2 %stay in bounds make 2 loop ghost cells 1,2, n-1 and n
    for j = 1:n
        %periodic boundayr conditions 
        del_u32R = u(k(j+2)) - u(k(j+1));
        del_u12R = u(k(j+1)) - u(k(j));
        del_u12L = u(k(j)) - u(k(j-1));
        del_u32L = u(k(j-1)) - u(k(j-2));
        
        %compute numerical speed of shock
        aR = 0.5 * (u(k(j+1)) + u(k(j)));
        aL = 0.5 * (u(k(j)) + u(k(j-1)));
       
        qR = R(del_u32R, del_u12L);
        qL = R(del_u12R, del_u32L);
        
        eps_2R = 0.5 * abs(aR) * qR;
        eps_2L = 0.5 * abs(aL) * qL;
        
        eps_4R = K*abs(aR) * (1-qR); %turn off at extremum and shocks
        eps_4L = K*abs(aL) * (1-qL);
        
        dR = eps_2R*del_u12R - eps_4R*(del_u32R - 2*del_u12R + del_u12L);
        dL = eps_2L*del_u12L - eps_4L*(del_u12R - 2*del_u12L + del_u32L);
        
        hR = 0.5 * (f(u(k(j+1))) + f(u(k(j)))) - dR;
        hL = 0.5 * (f(u(k(j))) + f(u(k(j-1)))) - dL;
        
        unew(k(j)) = u(k(j)) + lambda * (hL- hR);
    end
    %BC 0 Neuman
    % periodic boundary conditions:
%     if (IC == 2)
%         unew(1) = unew(end - 3);   % copy value from rightmost unknown to ghost cell on left
%         unew(2) = unew(end-2);
%         unew(end) = unew(4);
%         unew(end - 1) = unew(3);
%     end
   u = unew;
   plot(x,u,'o','MarkerSize',4)
   axis([-2 2 -1 2])
   pause(1e-4)
end 
plotExactandApproxSolution(u,x, 'JST', 'JST', t_n, IC)
end