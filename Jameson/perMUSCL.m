function perMUSCL(dx,t_n, IC, limiter) %this function solves the 1D Burger's equation
%u_t + uu_x = 0 on [-1,1] using the MUSCL (monotone upstream-centered
%schems for conservation laws)
x=-1:dx:1; %spatial vector
dt = dx / 4;
lambda = dt/ dx;
t = dt:dt:t_n;
u = computeInitialConditions(IC,x);
n = length(x); %we need 2 ghost cells second order method on each side
L = computeLimiter(limiter);
unew = u;
%initialize temporary solution
%time tieration explicit forward euler method
k = @(j) (mod(j-1,n) + 1);
for niter = 1:length(t)
    %for j = 3:n-2
    for j = 1:n
        del_u32R = u(k(j+2)) - u(k(j+1));
        del_u12R = u(k(j+1)) - u(k(j));
        del_u12L = u(k(j)) - u(k(j-1));
        del_u32L = u(k(j-1)) - u(k(j-2));
        
        uL_plus = u(k(j)) +  0.5 * L(del_u12R, del_u12L);
        uR_plus = u(k(j+1)) - 0.5 * L(del_u32R, del_u12R);
        
        uL_min = u(k(j-1)) + 0.5 * L(del_u12L, del_u32L);
        uR_min = u(k(j)) - 0.5 * L(del_u12R, del_u12L);
        
        %compute f+ and f- using EO splitting
        [f_plus, f_min] = computeEOsplitting(uL_plus, uR_plus);

        hR = f_plus + f_min;
        [f_plus, f_min] = computeEOsplitting(uL_min, uR_min);
        hL = f_plus + f_min;
        
        unew(k(j)) = u(k(j)) + lambda * (hL - hR); 
    end
%     if (IC == 2) %periodic
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
plotExactandApproxSolution(u,x, 'MUSCL', limiter, t_n, IC)
end