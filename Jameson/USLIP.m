function [u u_t t] = USLIP(dx, nt, dt, x, u, x_coord, limiter, upwindFlux, intK, lb, ub, t0) %this function solves the 1D Burger's equation
%u_t + uu_x = 0 on [-1,1] using the USLIP (upstream limited positive) in
%contrast to the SLIP (symmetric limited positive) schemes
lambda = dt/ dx;
t = t0 : dt : t0 + nt;
n = length(x); %we need 2 ghost cells second order method on each side
L = computeLimiter(limiter);
j = 2:n-1;
dR = zeros(length(j), 1);
dL = zeros(length(j), 1);
u_t = zeros(1, length(t));
cell = find(x < x_coord+1e-6 & x > x_coord - 1e-6, 1, 'first');
if (isempty(cell)) %make sure xcoord exists
     fprintf('Ind warning\n')
     return
end
%time tieration explicit forward euler method
for niter = 1:length(t) - 1
    u_t(niter) = u(cell);
    ghost_uL = interp1(x,u,lb - dx, 'linear', 'extrap');
    ghost_uR = interp1(x,u,ub + dx, 'linear', 'extrap');
    f = intK([ghost_uL; u ;ghost_uR]);
    if (upwindFlux) %flux is derivative of integral
        f = -(f(2:end-1) - f(1:end-2)) / dx;
    else
        f = -(f(3:end) - f(2:end-1)) / dx;
    end
    del_u12R = u(j+1) - u(j);
    del_u12L = u(j) - u(j-1);
    del_u32R = [u(4:n); ghost_uR] - u(j+1); %3:n
    del_u32L = u(j-1) - [ghost_uL; u(1:n-3)];

    %compute nuemrical speed of shock
    aL = (f(j) - f(j-1)) ./ max(abs(u(j) - u(j - 1)), 1e-8);
    %if 0 will multiply by 0 but ok since f returns 0 at these points
    % if u(j-1) = u(j) in flat region and so u(j) = u(j+1) since
    % monotonically nonincreasing
    aL = aL .* sign(u(j) - u(j - 1));
    %When they are equal numerator is 0 so divide by small positive
    %tolerance: otherwise stay 0
    aR = (f(j+1) -  f(j)) ./ max(abs(u(j+1) - u(j)), 1e-8);
    aR = aR .* sign(u(j + 1) - u(j));

    alphaR = 0.5*abs(aR);
    alphaL = 0.5*abs(aL);
    ind_pos =  aR > 0;
    dR(ind_pos) =  alphaR(ind_pos) .*(del_u12R(ind_pos)...
                                - L(del_u12R(ind_pos), del_u12L(ind_pos)));
    dR(~ind_pos) = alphaR(~ind_pos) .* (del_u12R(~ind_pos) ...
                                - L(del_u12R(~ind_pos), del_u32R(~ind_pos)));
    ind_pos = aL > 0;
    dL(ind_pos) = alphaL(ind_pos) .*(del_u12L(ind_pos) ...
                                - L(del_u12L(ind_pos), del_u32L(ind_pos)));

    dL(~ind_pos) = alphaL(~ind_pos) .* (del_u12L(~ind_pos)...
                                - L(del_u12L(~ind_pos), del_u12R(~ind_pos)));
    
    %average flux    
    hL = 0.5 * (f(j) + f(j-1)) - dL;
    hR = 0.5 * (f(j+1) + f(j)) - dR;

    u(j) = u(j) + lambda*(hL - hR);
end
    u_t(end) = u(cell);
end