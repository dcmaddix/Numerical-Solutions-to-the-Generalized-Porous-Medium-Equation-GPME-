%Mimetic 2016 paper implementing Scheme I for various ki functions using
%predictor corrector for timestepping to prevent locking. This function uses a
%predictor correction scheme where the fintie volume divergence form of the
%equation for k is solved at the half timestep, namely
%dk/dt - k^(1-1/n) d/dx(k^(1/n) dkdx) = 0.

function [p error] = mimeticTransportCoeffPredCorrFVk(k_i_keqtn, k_i_peqtn, N, k_pow, bdry_left, bdry_right,...
                                        F0, alpha, BC_type, timeScheme, nt, p0, delta_t, x_coord, errorType, ...
                                        upwind, weight)
%INPUTS
    %k_i_keqtn - function to approx cell-centered diffusion coeff k_tilde at node
           %point usually depends on the arithmetic and harmonic cell-centered
           %averages for the k estimate eqtn
    
    %k_i_peqtn - function to approx cell-centered diffusion coeff k_tilde at node
           %point usually depends on the arithmetic and harmonic cell-centered
           %averages for the pressure eqtn
    
    %N - number of interior grid points
    
    %k_pow - Assuming k = p^(k_pow)
         
    %bdry_left - function handle for left bdry function which grows in time
    
    %bdry_right - prescribed right boundary value for Dirichlet BC
            
    %F0 - forcing function
    
    %alpha - compressibility constant
    
    %BC_type - Dirichlet or Neumann for right bdry assuming Dirichlet on left bdry
    
    %timeScheme - ForwardEuler, BackwardEuler, CrankNicolson
    
    %nt - Final time
    
    %p0 - Initial Condition
    
    %delta_t - Timestep
    
    %nt_inner - Number of inner predictor iterations to do
    
    %x_coord - x position to store pressure at every time for pressure
                %versus time plot to check for oscillations
    
    %upwind - boolean if turned on upwind the (dk/dx)^2 term leading to
              %smoother and improved results. If turned off, do second order central
              %difference for this term
     
    %errorType - l2norm, maxnorm
              
%OUTPUTS:
    %p - pressure for plotting
    %error - returns 2x1 vector with pressure error and permability error,
    %respectively

%ALGORITHM:
    p = p0; %set pressure to initial solution
    F = F0; %initial source term
    h = 1/ N; %uniform grid spacing in 1D
    lb = 0.0; %domain lower bound
    ub = 1.0; %domain upper bound
    x = linspace(h/2,1-h/2, N); %to plot against velocity vector which is defined at cell faces
    %define x at cell-centers
    x_cent = [lb 0.5 * (x(1:end-1) + x(2:end)) ub]; %to plot against pressures defined at cell centers
    %store xcoord to plot at various times
    ind = find(x_cent < x_coord+1e-6 & x_cent > x_coord - 1e-6);
    if (length(ind) ~= 1)
        ind = ind(1);
    end
    k_exact = zeros(N+1,1);
    %time vector
    t = 0:delta_t:nt;
    %plotting p at specific x coord versus time
    p_t = zeros(1,length(t));
    p_t_exact = zeros(1,length(t)); %store exact solution vs time at specific xcoord
    %define arithmetic average function
    k_A = @(k1,k2) 0.5*(k1+k2);
    for i = 1:length(t)
        %write simple finite volume solver and separate BC
        p(1) = bdry_left(t(i));
        if (strcmp(BC_type, 'Dir'))
            %set Dirichlet RHS
            p(end) = bdry_right;
        end
        k = p.^(k_pow);
        p_t(i) = p(ind);
        %add in predicted pressure for right k
        %take half a timestep with k and then use that k for a full
        %timestep with p
        %important hoåw update the nonlinear k in time
        %analytical k
%             
        %define Ti which are the tridiagonal entries of the matrix. Formula
        %simplifies to below for uniform spacing k
        k_plus = (k_i_keqtn(k(2:end-1).^(1/k_pow), k(3:end).^(1/k_pow)).^2 ...
                ./ k_A(k(2:end-1).^(1/k_pow), k(3:end).^(1/k_pow)))';
            
        k_minus = (k_i_keqtn(k(1:end-2).^(1/k_pow), k(2:end-1).^(1/k_pow)).^2 ...
                ./ k_A(k(2:end-1).^(1/k_pow), k(1:end-2).^(1/k_pow)))';
        %reorder for matrix off-diagonal
        k_minus_new = k_minus;
        k_minus_new(1:end-1) = k_minus(2:end);
        k_plus_new = k_plus;
        k_plus_new(2:end) = k_plus(1:end-1); %to get spdiag to work
       
        %tridiagonal as specified by the second argument changed sign here
        D_k = spdiags([-(k_minus_new)', (k_minus + k_plus)'...
                        , (-k_plus_new)'], [-1 0 1], N-1, N-1);
        %3 term recursion (T_i + T_(i-1)) k_i - T_i k_i+1 - T_(i-1)_k_(i-1)
        %homogenous Neumann sub into last equation k(n+1) = k(n)
        %T_i cancels and get simplified euqaiton T_(i-1)k_i -
        %T_(i-1)k_(i-1)
        if (strcmp(BC_type, 'Neu')) 
            D_k(N-1,N-1) = -D_k(N-1,N-2); %k = p^1/3 p(end) = p(end-1) transfers over for j too
        elseif (strcmp(BC_type, 'Dir')) %end-1 for upwinding
            if (upwind) %k(1:end-2)
                k_coeff = k(end-2);
            else %central  k(2:end-1)
                k_coeff = k(end-1);
            end
            F(end) = F(end) + k_coeff ^((k_pow - 1) / k_pow)*k_plus(end)*k(end) / h^2; %Dirichlet RHS
        end
        
        if (upwind)
            k_coeff = k(1);
        else %central scheme
            k_coeff = k(2);
        end
        
        F(1) = F(1) + k_coeff^((k_pow-1) / k_pow) * k_minus(1)*k(1) / h^2;

        if (upwind) %changing k coeff discretization
            k(2:end-1) = k(2:end -1) + ((weight*delta_t) / alpha)*(F - k(1:end-2).^((k_pow-1)/k_pow).*...
                       (1/h^2).*(D_k*k(2:end-1))); %explicit update for k estimate
        else
            k(2:end-1) = k(2:end -1) + ((weight*delta_t) / alpha)*(F - k(2:end-1).^((k_pow-1)/k_pow).*...
                       (1/h^2).*(D_k*k(2:end-1))); %explicit update for k estimate
        end

        F = F0; %reset to inital source term values

        %%%Solve for exact k%%%
        k_exact(x_cent < t(i)) = 3*(t(i)-x_cent(x_cent < t(i)));
        k_exact(x_cent >= t(i)) = 1e-9;
        p_exact = k_exact.^(1/k_pow);
        p_t_exact(i) = p_exact(ind);
        
        %%%%%%%%%%%End k solve%%%%%%%%%%%%%%%
        k_plus = (k_i_peqtn(k(2:end-1), k(3:end)).^2 ./ k_A(k(2:end-1), k(3:end)))';
        k_minus = ( k_i_peqtn(k(1:end-2), k(2:end-1)).^2 ./ k_A(k(2:end-1), k(1:end-2)))';
        %reorder for matrix off-diagonal
        k_minus_new = k_minus;
        k_minus_new(1:end-1) = k_minus(2:end);
        k_plus_new = k_plus;
        k_plus_new(2:end) = k_plus(1:end-1); %to get spdiag to work
        D_k = spdiags([-k_minus_new', (k_minus + k_plus)', -k_plus_new'], [-1 0 1], N-1, N-1);
        %Now use new k to solve pressure eqtn

        if (strcmp(timeScheme, 'BackwardEuler')) %Baskward Euler
            M = alpha/delta_t*speye(N-1,N-1) + (1/h^2)*D_k;
        elseif (strcmp(timeScheme, 'CrankNicolson')) %Crank Nicolson
            M = alpha/delta_t*speye(N-1,N-1) + (0.5/h^2)*D_k;
        end

        %Add Dirichlet BC contribution to RHS-fix boundary by dividing by
        %h^2 as done with matrix term D_k
        F(1) = F(1) + k_minus(1)*p(1) / h^2;

        if (strcmp(BC_type, 'Dir'))
            F(end) = F(end) + k_plus(end)*p(end) / h^2; %Dirichlet RHS
        end

        %solve (N-1)x(N-1) linear system for interior pressures
        if (strcmp(timeScheme, 'ForwardEuler'))
            p(2:end-1) = p(2:end-1) + (delta_t / alpha)*(F - (1/h^2)*D_k*p(2:end-1));

        elseif(strcmp(timeScheme, 'BackwardEuler'))
            p(2:end-1) = M \ ((alpha / delta_t)*p(2:end-1) + F); %check RHS correct
        elseif(strcmp(timeScheme,'CrankNicolson'))
            p(2:end-1) = M \ ((alpha/delta_t)*p(2:end-1) + F - (0.5/h^2)*D_k*p(2:end-1));
        end

        %enforce 0 Neumann BC
        if (strcmp(BC_type,'Neu'))
            p(end) = p(end-1);
        end
        F = F0; %set back at time step otherwise accumulating BC terms
    end
    if (strcmp(errorType, 'l2norm'))
        error(1) = h^0.5*norm(p_exact - p);
        error(2) = h^0.5*norm(k_exact - k);
    elseif (strcmp(errorType, 'maxnorm'))
        [max_p ind_p] = max(abs(p_exact - p));
        error(1) = max_p / max(abs(p_exact));
        [max_k ind_k] = max(abs(k_exact - k));
        error(2) = max_k / max(abs(k_exact));
    else
        fprintf('Undefined Error Measure Type\n')
        return
    end
    
    figure(1)
    hold on;
    plot(x_cent,p,'*',x_cent,p_exact,'*')
    legend('p','p_{exact}')
    title(sprintf('Plot of numerical pressure versus exact pressure'))
    xlabel('position x')
    ylabel('pressure p')
    
    figure(2);
    plot(t,p_t,t, p_t_exact)
    xlabel('time t')
    ylabel(sprintf('pressure at position x = %.1f', x_coord))
    title(sprintf('Plot of pressure versus time at x = %.1f', x_coord))
    legend('numerical pressure','exact pressure')
end
