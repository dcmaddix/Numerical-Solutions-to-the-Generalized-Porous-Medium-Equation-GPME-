%Mimetic 2016 paper implementing Scheme I for various ki functions using
%predictor corrector for timestepping to prevent locking
function p = mimeticTransportCoeffPredCorrPress(N, k_pow, bdry_left, bdry_right,...
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
    F = F0; %initial source term
    h = 1/ N;
    lb = 0.0;
    ub = 1.0;
    x = linspace(h/2,1-h/2, N); %to plot against velocity
    %define x at cell-centers
    x_cent = [lb 0.5 * (x(1:end-1) + x(2:end)) ub]; %to plot against pressures
    k_exact = zeros(N+1,1);
    %define spatial and temporal spacing
    %define domain bounds
    %define averages-take arithmetic in first half timestep as a predictor and then
    %harmonic in second
    k_H = @(k1,k2) (2*k1.*k2) ./ (k1+k2);
    k_A = @(k1,k2) (k1+k2) / 2;
    k_FV = @(k1,k2) (k_H(k1,k2).*k_A(k1,k2)).^(0.5);
    for t = 0: delta_t : nt
        %write simple finite volume solver and separate BC
        p(1) = bdry_left(t);
        if (strcmp(BC_type, 'Dir'))
            %set Dirichlet RHS
            p(end) = bdry_right;
        end
        pold = p;
        for t_inner = 2:-1:1 %first take half a timestep and use that k for the final
            k = p.^(k_pow);
            %Solve eqtn at 1/2 a timestep for k
            %%%%%%%%%%%End k solve%%%%%%%%%%%%%%%
            if (t_inner == 2)
                k_i = k_A;
            else
                k_i = k_FV;
            end
            k_plus = 2*(k_i(k(2:end-1), k(3:end)).^2 ./ (k(2:end-1) + k(3:end)))';
            k_minus = 2*( k_i(k(1:end-2), k(2:end-1)).^2 ./ (k(2:end-1) + k(1:end-2)))';
            %reorder for matrix off-diagonal
            k_minus_new = k_minus;
            k_minus_new(1:end-1) = k_minus(2:end);
            k_plus_new = k_plus;
            k_plus_new(2:end) = k_plus(1:end-1); %to get spdiag to work
            D_k = spdiags([-k_minus_new', (k_minus + k_plus)', -k_plus_new'], [-1 0 1], N-1, N-1);
            %Now use new k to solve pressure eqtn

            if (strcmp(timeScheme, 'BackwardEuler')) %Baskward Euler
                M = alpha/(delta_t / t_inner)*speye(N-1,N-1) + (1/h^2)*D_k;
            elseif (strcmp(timeScheme, 'CrankNicolson')) %Crank Nicolson
                M = alpha/(delta_t / t_inner)*speye(N-1,N-1) + (0.5/h^2)*D_k;
            end

            %Add Dirichlet BC contribution to RHS-fix boundary by dividing by
            %h^2 as done with matrix term D_k
            F(1) = F(1) + k_minus(1)*p(1) / h^2;

            if (strcmp(BC_type, 'Dir'))
                F(end) = F(end) + k_plus(end)*p(end) / h^2; %Dirichlet RHS
            end

            %solve (N-1)x(N-1) linear system for interior pressures
            if (strcmp(timeScheme, 'ForwardEuler'))
                p(2:end-1) = pold(2:end-1) + ((delta_t / t_inner) / alpha)*(F - (1/h^2)*D_k*pold(2:end-1));

            elseif(strcmp(timeScheme, 'BackwardEuler'))
                p(2:end-1) = M \ ((alpha / (delta_t / t_inner))*pold(2:end-1) + F); %check RHS correct
            elseif(strcmp(timeScheme,'CrankNicolson'))
                p(2:end-1) = M \ ((alpha/(delta_t / t_inner))*pold(2:end-1) + F - (0.5/h^2)*D_k*pold(2:end-1));
            end

            %enforce 0 Neumann BC
            if (strcmp(BC_type,'Neu'))
                p(end) = p(end-1);
            end
            F = F0; %set back at time step otherwise accumulating BC terms
        end
        %compute exact k
        k_exact(x_cent < t) = 3*(t-x_cent(x_cent < t));
        k_exact(x_cent >= t) = 1e-9;
        p_exact = k_exact.^(1/k_pow);
    end
    h*norm(p_exact-p)
    figure(10)
    hold on;
    plot(x_cent,p,'*',x_cent,p_exact,'*')
    legend('p','p_{exact}')
    title(sprintf('Plot of pressure FV with harmonic averaging'))
end
