function u = MUSCL(dx,t_n, IC, limiter, u, f, x) %this function solves the 1D Burger's equation
%u_t + uu_x = 0 on [-1,1] using the MUSCL (monotone upstream-centered
%schems for conservation laws)
%x = computeDomain(dx,IC);
dt = dx^2;
lambda = dt/ dx;
t = dt:dt:t_n;
%u = computeInitialConditions(IC,x);
n = length(x); %we need 2 ghost cells second order method on each side
L = computeLimiter(limiter);
unew = u;
%initialize temporary solution
%time tieration explicit forward euler method
for niter = 1:length(t)
    for j = 3:n-2
        del_u32R = u(j+2) - u(j+1);
        del_u12R = u(j+1) - u(j);
        del_u12L = u(j) - u(j-1);
        del_u32L = u(j-1) - u(j-2);
        
        uL_plus = u(j) +  0.5 * L(del_u12R, del_u12L);
        uR_plus = u(j+1) - 0.5 * L(del_u32R, del_u12R);
        
        uL_min = u(j-1) + 0.5 * L(del_u12L, del_u32L);
        uR_min = u(j) - 0.5 * L(del_u12R, del_u12L);

        %compute f+ and f- using EO splitting
        [f_plus, f_min] = computeEOsplitting(uL_plus, uR_plus, f, j);

        hR = f_plus + f_min;
        [f_plus, f_min] = computeEOsplitting(uL_min, uR_min, f, j);
        hL = f_plus + f_min;
        
        unew(j) = u(j) + lambda * (hL - hR); 
    end
    unew = computeBC(unew,IC);
    u = unew;
%     
%     plot(x,u,'o','MarkerSize',4)
%    axis([-2 2 -1 2])
%   pause(1e-4)
end
plotExactandApproxSolution(u,x, 'MUSCL', limiter, t_n, IC)
end