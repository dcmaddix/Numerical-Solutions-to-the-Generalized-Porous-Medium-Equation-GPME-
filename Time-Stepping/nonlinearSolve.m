%returns function G of the nonlinear solve for implicit schemes
function G = nonlinearSolve(p, k_p, k_i, p_old,alpha,p1,pend, F0, timeScheme, N,BC_type, delta_t)
    %define spatial and temporal spacing
    h = 1 / N; %uniform grid spacing in 1D
    k = k_p([p1; p; pend]);
    k_plus = 2*(k_i(k(2:end-1), k(3:end)).^2 ./ (k(2:end-1) + k(3:end)))';
    k_minus = 2*( k_i(k(1:end-2), k(2:end-1)).^2 ./ (k(2:end-1) + k(1:end-2)))';
    %reorder for matrix off-diagonal
    k_minus_new = k_minus;
    k_minus_new(1:end-1) = k_minus(2:end);
    k_plus_new = k_plus;
    k_plus_new(2:end) = k_plus(1:end-1); %to get spdiag to work
    
    %tridiagonal as specified by the second argument changed sign here
    D_k = spdiags([-k_minus_new, k_minus + k_plus, -k_plus_new], [-1 0 1], N-1, N-1);
    %3 term recursion (T_i + T_(i-1)) p_i - T_i p_i+1 - T_(i-1)_p_(i-1)
    %homogenous Neumann sub into last equation p(n+1) = p(n)
    %T_i cancels and get simplified euqaiton T_(i-1)p_i -
    %T_(i-1)p_(i-1)
    if (strcmp(BC_type, 'Neu')) 
        D_k(N-1,N-1) = -D_k(N-1,N-2); 
    end
    F0(1) = F0(1) + k_minus(1)*p1/h^2; %need to reset to 0 otherwise keep adding on
    if (strcmp(BC_type, 'Dir'))
       F0(end) = F0(end) + k_plus(end)*pend / h^2; %Dirichlet RHS
    end
    %fullt implicit timesteps
    if (strcmp(timeScheme, 'BackwardEuler')) %Backward Euler
        G = (alpha / delta_t)*(p - p_old) + (1/h^2)*D_k*p - F0;
    elseif (strcmp(timeScheme, 'CrankNicolson')) %Crank Nicolson
        G = (alpha / delta_t)*(p - p_old) + (0.5/h^2)*D_k*p - F0 + (0.5/h^2)*D_k*p_old;
    end
end