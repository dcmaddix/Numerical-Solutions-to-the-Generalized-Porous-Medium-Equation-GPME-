function testLinearK(N,k_lb) %input the lower bound on the k spectrum here k is linear
    %tests case where k is smoothly varying in this case linearly, so we
    %are solving d/dx(k(x)dp/dx) = 1 with variable coefficients subject to
    %homogenous Dirichlet BC p(0) = p(1) =1 on the domain [0,1]
    h = 1 / N; %uniform grid spacing
    lb = 0;
    ub = 1;
    k = linspace(k_lb,1,N+1);
    F = ones(N-1,1); %constant forcing term
    %call scheme I with homogenous Dirichlet BC at both endpoints
    bdry = zeros(2,1);
    [p_harm p_arith p_FV, ~, u_harm u_arith u_FV] = testSchemeI(N,k, bdry, F);
    %for plotting
    %define x at gridpoints
    x = linspace(h/2,1-h/2, N); %to plot against velocity
    %define x at cell-centers
    x_cent = [lb 0.5 * (x(1:end-1) + x(2:end)) ub]; %to plot against pressures
    %plot exact pressure solution
    p_exact = -(1 / (log(k_lb) * (1-k_lb))) * log((1-k_lb)*x_cent + k_lb) -...
         x_cent / (1-k_lb)  + 1 / (1-k_lb);
    %plot exact velocity solution
    u_exact = -(-x -k_lb/(1-k_lb) - 1/log(k_lb)) ./ ((1-k_lb)*x+k_lb);
     
    figure(1);
    plot(x_cent,p_harm, x_cent, p_arith, x_cent, p_FV, x_cent, p_exact)
    title('Plot of pressure from Mimetic Scheme I with various k_i')
    legend('Harmonic', 'Arithmetic', 'Finite Volume with Harmonic', 'Exact Solution')
    xlabel('x')
    ylabel('pressure p')
    
    figure(2);
    plot(x, u_harm, x, u_arith, x, u_FV,x, u_exact)
    title('Plot of velocity Mimetic Scheme I with various k_i')
    legend('Harmonic', 'Arithmetic', 'Finite Volume with Harmonic', 'Exact Solution')
    xlabel('x')
    ylabel('velocity u')
    
    figure(3);
    plot(x_cent,k)
    xlabel('x')
    ylabel('permeability k')
    title('Plot of permeabilities versus position')
    
end
