%Mimetic 2016 paper implementing Scheme I for various ki functions
function p = mimeticTransportCoeffNonlinear(k_i, N, k_p, bdry_left, bdry_right,...
                                        F0, alpha, BC_type, timeScheme, nt, p0, delta_t)
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
    %BC_type Dirichlet or Neumann
    %assuming Dirichlet on left bdry
    %on right could be Neumann or Dirichlet
     
    %initial condition
    p = p0;
    
    for t = delta_t: delta_t : nt
        %write simple finite volume solver and separate BC
        p(1) = bdry_left(t); %compute boundaries at NEW timestep to used in matrix
        %p(1) = 10;
        if (strcmp(BC_type, 'Dir'))
            %set Dirichlet RHS
            p(end) = bdry_right(t);
        end
        G = @(p_new) nonlinearSolve(p_new, k_p, k_i, p(2:end-1),alpha,p(1),p(end), F0, timeScheme, N, BC_type, delta_t);
        eps = 0;
        p(2) = p(2) + eps;
        p(3) = p(3) + eps;
        p(2:end-1) = fsolve(G,p(2:end-1));
        %enforce 0 Neumann BC
        if (strcmp(BC_type,'Neu'))
            p(end) = p(end-1);
        end
    end
end
