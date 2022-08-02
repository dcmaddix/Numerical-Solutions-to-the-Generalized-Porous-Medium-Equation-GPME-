%Mimetic 2016 paper implementing Scheme I for various ki functions for the
%steady transport equation, i.e. a Poisson equation solver
%div(k*grad(p)) = F with Dirichlet boundary conditions
function [p u] = mimeticSteadyTransportCoeff(k_i, N, k, bdry, F)

    %INPUTS
    %k_i - function to approx cell-centered diffusion coeff k_tilde at node
    %point usually depends on the arithmetic and harmonic cell-centered
    %averages
    %N - number of interior grid points
    %k - array of constant cell-centered permeability coefficients
    %(k_tilde^(i-1/2) for cell i-1/2 and k_tilde^(1+1/2) for cell i + 1/2 as
    %defined in the paper
    %bdry - 2x1 array giving the prescribed Dirichlet pressure values at the
    %endpoints of the interval
    %F - forcing function
    
    h = 1 / N; %uniform grid spacing in 1D
    
    %solving for pressures-initialize array
    p = zeros(N+1,1);
    
    %specify Dir BC at endpoints
    p(1) = bdry(1);
    p(end) = bdry(2);
   
    %define Ti which are the tridiagonal entries of the matrix. Formula
    %simplifies to below for unitorm spacing k
    k_plus = 2*(k_i(k(2:end-1), k(3:end)).^2 ./ (k(2:end-1) + k(3:end)))';
    k_minus = 2*( k_i(k(1:end-2), k(2:end-1)).^2 ./ (k(2:end-1) + k(1:end-2)))';
    
    %reorder for matrix off-diagonal for spdiag function
    k_minus_new = k_minus;
    k_minus_new(1:end-1) = k_minus(2:end);
    k_plus_new = k_plus;
    k_plus_new(2:end) = k_plus(1:end-1);
    
    %tridiagonal as specified by the second argument
    D_k = spdiags([-k_minus_new, k_minus + k_plus, -k_plus_new], [-1 0 1], N-1, N-1);
    
    %multiply RHS by h^2 before adding boundary contribution
    F = h^2 * F;
    
    %Add Dirichlet BC contribution to RHS
    F(1) = F(1) + k_minus(1)*p(1);
    F(end) = F(end) + k_plus(end)*p(end);
    
    %solve (N-1)x(N-1) linear system for interior pressures
    %Poisson solver
    p(2:end-1) = D_k \ F;
    
    %use p to construct ui = -GRAD(ph)i at node points-slide 12
    u = (2*k_i(k(1:end-1), k(2:end))' .* (p(1:end-1) - p(2:end))) ./ (h*(k(1:end-1) + k(2:end))');
end