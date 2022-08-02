function [u u_t t] = JST(dx,nt, dt, x, u, x_coord, upwindFlux, intK, lb, ub, t0) %this function solves the 1D Burger's equation 
%output u_t as well
%u_t + uu_x = 0 on [-1,1] using the JST scheme
%u-initial condition at time t = 0
%dx-spatial timestep use with CFL condition to compute timestep
%t_n is the final timestep that we evolve to
eps = 10^-3;
K = 1/2;
q = 3;
lambda = dt / dx;
t = t0 : dt : t0 + nt;
%initial burgers condition
%u = computeInitialConditions(IC,x);
n = length(x); %we need 2 ghost cells second order method on each side
r = 3/2;
R = @(u,v) abs( (u-v) ./ (max( (abs(u) + abs(v)), eps*dx^r))) .^ q;
%if u and v have opposite signs like at shocks or extremum R will have
%value 1
j = 2:n-1;
cell = find(x < x_coord+1e-6 & x > x_coord - 1e-6, 1, 'first');
if (isempty(cell)) %make sure xcoord exists
     fprintf('Ind warning\n')
     return
end
u_t = zeros(1, length(t));
for niter = 1:length(t) - 1 %time integration
    u_t(niter) = u(cell);
    ghost_uL = interp1(x,u,lb - dx, 'linear', 'extrap');
    ghost_uR = interp1(x,u,ub + dx, 'linear', 'extrap');
    %Onlt need to compute integral once
    f = intK([ghost_uL; u ;ghost_uR]); %compute integral for each point only once and shift
    if (upwindFlux) %flux is derivative of integral
        f = -(f(2:end-1) - f(1:end-2)) / dx;
    else
        f = -(f(3:end) - f(2:end-1)) / dx;
    end
    
    del_u12R = u(j+1) - u(j);
    del_u12L = u(j) - u(j-1);
    del_u32R = [u(4:n); ghost_uR] - u(j+1); %3:n
    del_u32L = u(j-1) - [ghost_uL; u(1:n-3)];

    %compute numerical speed of shock: 
    aL = (f(j) - f(j-1)) ./ max(abs(u(j) - u(j - 1)), 1e-8);
    %When they are equal numerator is 0 so divide by small positive
    %tolerance: otherwise stay 0
    aR = (f(j+1) -  f(j)) ./ max(abs(u(j+1) - u(j)), 1e-8);

    qR = R(del_u32R, del_u12L);
    qL = R(del_u12R, del_u32L);

    eps_2R = 0.5 * abs(aR) .* qR;
    eps_2L = 0.5 * abs(aL) .* qL;

    eps_4R = K * abs(aR) .* (1-qR); %turn off at extremum and shocks
    eps_4L = K * abs(aL) .* (1-qL);

    dR = eps_2R.*del_u12R - eps_4R.*(del_u32R - 2*del_u12R + del_u12L);
    dL = eps_2L.*del_u12L - eps_4L.*(del_u12R - 2*del_u12L + del_u32L);

    %average flux    
    hL = 0.5 * (f(j) + f(j-1)) - dL;
    hR = 0.5 * (f(j+1) + f(j)) - dR;

    u(j) = u(j) + lambda * (hL - hR);
end 
u_t(end) = u(cell);
%boundaries are Dircihlet and fixed in time so do not change-Only update
%interior
end