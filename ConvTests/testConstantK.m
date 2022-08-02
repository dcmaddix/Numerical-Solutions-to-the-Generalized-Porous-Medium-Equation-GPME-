function testConstantK(N)
    %this function where N is the number of gridpoints
    %tests the case where k = 1 is constant
    %in 1D solving d/dx(-kdp/dx) = 1 k = 1 => -d^2p/dx^2 = 1
    %subject to homogenous Dirichlet BC p(0) = p(1) = 0
    %define region omega by upper and lower bounds
    h = 1 / N; %uniform grid spacing
    lb = 0;
    ub = 1;
    k = ones(N+1,1)';
    F = ones(N-1, 1); %constant source term
    %we have N interior grid points k is the cell centered permeability
    %constants. So there are N-1 on the interior and to extend to the boundary
    %we get N+1 total points

    %call scheme I with homogenous Dirichlet BC at both endpoints
    bdry = zeros(2,1);
    [p_harm p_arith p_FV, ~, u_harm u_arith u_FV] = testSchemeI(N,k, bdry, F);
    
    %define x at gridpoints
    x = linspace(h/2,1-h/2, N); %to plot against velocity
    %define x at cell-centers
    x_cent = [lb 0.5 * (x(1:end-1) + x(2:end)) ub]; %to plot against pressures
    %compare to exact solution which is a simple quadratic
    p_exact = 0.5*(-x_cent.^2 + x_cent)'; %only in constant coeff case
    
    figure(1);
    plot(x_cent,p_harm, x_cent, p_arith, x_cent, p_FV, x_cent, p_exact)
    title('Plot of pressure from Mimetic Scheme I with various k_i versus the exact solution')
    legend('Harmonic', 'Arithmetic', 'Finite Volume with Harmonic', 'Exact')
    xlabel('x')
    ylabel('pressure p')
    %1st order in this metric as decrease error by 1/2 h decreases by 1/2
    %compute relative errors as measured by Euclidean 2-norm
    %Note that in paper achieve 2nd order convergence in P usign their
    %special Ph norm
    press_error_FV = norm(p_FV - p_exact)/ norm(p_FV)
    press_error_harm = norm(p_harm - p_exact) / norm(p_exact)
    press_error_arith = norm(p_arith - p_exact) / norm(p_arith)
    
    %plot u = -grad(p) as defined at nodes rather than cell centers
    u_exact = (x - 0.5)';
    figure(2);
    plot(x, u_harm, x, u_arith, x, u_FV, x, u_exact)
    title('Plot of velocity Mimetic Scheme I with various k_i versus the exact solution')
    legend('Harmonic', 'Arithmetic', 'Finite Volume with Harmonic','Exact')
    xlabel('x')
    ylabel('velocity u')
    %first order as measured by relative two norm
    vel_error_FV = norm(u_FV - u_exact) / norm(u_FV) %first order
    vel_error_harm = norm(abs(u_harm - u_exact)) / norm(u_harm)
    vel_error_arith = norm(abs(u_arith - u_exact)) / norm(u_arith)
end