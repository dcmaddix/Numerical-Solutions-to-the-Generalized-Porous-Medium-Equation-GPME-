%Mimetic 2016 paper implementing Scheme I for various ki functions using
%predictor corrector for timestepping to prevent locking for cell-centered
%permeabilites
%Solving alpha*dp/dt -  d/dx(kdp/dx) = F
%p = k^(1/n) leads to predictor equation for k namely
%alpha(dk/dt) - (1/n)(dk/dx)^2 -kd^2k/dx^2 = F
function [p error] = mimeticTransportCoeffPredCorr(k_i, N, k_pow, ...
                                                   bdry_left, bdry_right,...
                                                   F0, alpha, BC_type, ...
                                                   timeScheme, nt, p0, ...
                                                   delta_t, h, nt_inner, x_coord, ...
                                                   upwind, errorType, use_exact)
%INPUTS:
    %k_i - function to approx cell-centered diffusion coeff k_tilde at node
           %point usually depends on the arithmetic and harmonic cell-centered
           %averages
    
    %N - number of interior grid points
    
    %k_pow - Assuming k = p^(k_pow)
         
    %bdry_left - function handle for left bdry function which grows in time
    
    %bdry_right - function handle for right boundary value for Dirichlet BC
            
    %F0 - forcing function
    
    %alpha - compressibility constant
    
    %BC_type - Dir or Neu for right bdry assuming Dirichlet on left bdry
    
    %timeScheme - ForwardEuler, BackwardEuler, CrankNicolson
    
    %nt - Final time
    
    %p0 - Initial Condition
    
    %delta_t - Timestep
    
    %h - Spatial step
    
    %nt_inner - Number of inner predictor iterations to do
    
    %x_coord - x position to store pressure at every time for pressure
                %versus time plot to check for oscillations
    
    %upwind - boolean if turned on upwind the (dk/dx)^2 term leading to
              %smoother and improved results. If turned off, do second order central
              %difference for this term
     
    %errorType - l2norm, maxnorm
    
    %use_exact - Boolean Flag as whether or not to use exact k for
                %unlocking harmonic
              
%OUTPUTS:
    %p - pressure for plotting
    
    %error - returns 2x1 vector with pressure error and permability error
   
%Algorithm:
    p = p0; %set pressure to initial solution
    F = F0; %initial source term
    
    lb = 0.0; %domain lower bound
    ub = 1.0; %domain upper bound
    
    %define spatial vector at cell faces to plot against velocity
    x = linspace(lb + h/2,ub - h/2, N);
    
    %define x at cell-centers at cell centers to plot against pressure
    x_cent = [lb 0.5 * (x(1:end-1) + x(2:end)) ub];
    
    %Plot pressure at this x-coordinate for each time to get pressure vs
    %time plot
    ind = find(x_cent < x_coord+1e-6 & x_cent > x_coord - 1e-6);
    if (length(ind) > 1)
        ind = ind(1);
    elseif (isempty(ind)) %make sure xcoord exists
        ind = 1;
        fprintf('Ind warning\n')
    end
    
    %define for exact solution
    eps = 1e-9;
    %define k_A as arithmetic average
    k_A = @(k1,k2) 0.5 * (k1 + k2);
    k_exact = zeros(N+1,1)'; %compare against exact k
    p_exact = zeros(N+1,1);
    k = zeros(N+1,1)';
    
    t = 0: delta_t: nt; %time vector where initial condition is given at t=0
    p_t = zeros(1,length(t)); %plotting p at specific x coord versus time
    p_t_exact = zeros(1,length(t));
    %compute constant lambda
    lambda = delta_t / h^2
    
    for i = 1:length(t) - 1
        %write simple finite volume solver and separate BC
        pold = p; %prediction corrector store old p value
        kold = pold.^k_pow;
        %only depends on i so does not need to be in inner loop
        p_exact(x_cent < t(i)) = (k_pow*(t(i)-x_cent(x_cent < t(i))) ).^(1/k_pow);
        p_exact(x_cent >= t(i)) = eps^(1/k_pow); %see if fixes harmonic problem
        p_t_exact(i) = p_exact(ind); %store at every timestep
        p_t(i) = p(ind); %store at initial timestep need to add at end
        k = kold';
        %Inner t-predictor loop only take 1/2 a timestep rather than a full
        %step: full seems to work better
        for t_inner = 1:nt_inner
            %important how update the nonlinear k in time
            %k(1) = p(1)^k_pow; %at the old timestep
            %k(end) = p(end)^k_pow;
            k_avg = kold; %can avergae old and new time
            %k_avg = p.^k_pow;
            %change to central difference (k(3:end) - k(1:end-2)) / 2h
            %rather than (k(2:end-1) - k(1:end-2)) / h for upwinding
            if (upwind) %check sign
                speed = k_avg(2:end-1) - k_avg(1:end-2);
                for j = 1:length(speed)
                    if (speed(j) < 0)
                        firstDeriv(j) = speed(j);
                    else
                        firstDeriv(j) = k_avg(j+2) - k_avg(j+1);
                    end
                end
                firstDeriv = k_avg(2:end-1) - k_avg(1:end-2);
            else %use second order central differencing scheme
                firstDeriv = 0.5* (k_avg(3:end) - k_avg(1:end-2));
                speed = firstDeriv;
            end
            
            %Solve PDE for k as predictor: %add in predicted pressure for
            %right k
            %k(2:end-1) = kold(2:end-1) + lambda *( (1/k_pow) * speed.*firstDeriv ...
                        %+ k_avg(2:end-1).* (k_avg(3:end) - 2*k_avg(2:end-1) + k_avg(1:end-2)));
            %k(2:end-1) = kold(2:end-1) + lambda *((1/k_pow)*speed.*firstDeriv ...
                        %+ k_avg(1:end-2).*(k_avg(3:end) - 2*k_avg(2:end-1) + k_avg(1:end-2)));
            k(2:end-1) = k_avg(2:end-1) + lambda/t_inner *(k_avg(3:end) - 2*k_avg(2:end-1) + k_avg(1:end-2));
                 
            %k(1) = bdry_left(t(i+1))^k_pow;
            %k(end) = bdry_right(t(i+1))^k_pow;
            
            %Compute exact solution for plotting and testing       
            if (use_exact)
                if (strcmp(timeScheme, 'BackwardEuler')) %need to use future k values- compute at next t
                    ind_exact = x_cent < t(i+1);
                    ind_comp = x_cent >= t(i+1);
                    t_coord = t(i+1);
                elseif (strcmp(timeScheme, 'ForwardEuler'))
                    ind_exact = x_cent < t(i);
                    ind_comp = x_cent >= t(i);
                    t_coord = t(i);
                elseif(strcmp(timeScheme, 'CrankNicolson'))
                    ind_exact = x_cent < 0.5*(t(i) + t(i+1));
                    ind_comp = x_cent >= 0.5*(t(i) + t(i+1));
                    t_coord = 0.5*(t(i) + t(i+1));
                end
                k(ind_exact) = k_pow * (t_coord - x_cent(ind_exact));
                k(ind_comp) = eps; %see if fixes harmonic problem
                %k = ones(N+1,1)'; use for heat eqtn test case
            end
            k_exact(x_cent < t(i+1)) = k_pow*(t(i+1)-x_cent(x_cent < t(i+1)));
            k_exact(x_cent >= t(i+1)) = eps; %see if fixes harmonic problem
            p_exact = (k_exact.^(1/k_pow))';
            %define Ti which are the tridiagonal entries of the matrix. Formula
            %simplifies to below for uniform spacing k
            %k = k_exact';
            k_plus = (k_i(k(2:end-1), k(3:end)).^2 ./ k_A(k(2:end-1), k(3:end)))';
            k_minus = (k_i(k(1:end-2), k(2:end-1)).^2 ./ k_A(k(2:end-1), k(1:end-2)))';
            
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

            if (strcmp(timeScheme, 'BackwardEuler')) %Baskward Euler
                M = alpha * speye(N-1,N-1) + lambda * D_k;
            elseif (strcmp(timeScheme, 'CrankNicolson')) %Crank Nicolson
                M = alpha * speye(N-1,N-1) + 0.5 * lambda * D_k;
            end

            %Add Dirichlet BC contribution to RHS-fix boundary by dividing by
            %h^2 as done with matrix term D_k
            bdry_term = zeros(2,1);
            if (strcmp(timeScheme, 'ForwardEuler'))
                bdry_term(1) = p(1); %explicit use at current timestep
                bdry_term(2) = p(end);
            elseif(strcmp(timeScheme, 'BackwardEuler'))
                bdry_term(1) = bdry_left(t(i+1)); %p(1) at next timestep
                bdry_term(2) = bdry_right(t(i+1));
            elseif(strcmp(timeScheme, 'CrankNicolson'))
                %average of bdry value at current and new timestep-one from
                %each matrix to add to first and last components of F
                bdry_term(1) = 0.5 * (p(1) + bdry_left(t(i+1)));
                bdry_term(2) = 0.5 * (p(end) + bdry_right(t(i+1))); 
            end
            F(1) = F(1) + k_minus(1)*bdry_term(1) / h^2;

            if (strcmp(BC_type, 'Dir'))
                F(end) = F(end) + k_plus(end)*bdry_term(2) / h^2; %Dirichlet RHS
            end

            %solve (N-1)x(N-1) linear system for interior pressures
            if (strcmp(timeScheme, 'ForwardEuler'))
                p(2:end-1) = pold(2:end-1) + (1 / alpha)*(delta_t * F - lambda*D_k*pold(2:end-1));
            elseif(strcmp(timeScheme, 'BackwardEuler'))
                %constant F in time otherwise need future value
                p(2:end-1) = M \ (alpha * pold(2:end-1) + delta_t * F); %check RHS correct
            elseif(strcmp(timeScheme,'CrankNicolson'))
                %constant F in time otherwise need tie midpoint value 1/2(f(xi,tn) + f(xi,tn+1));
                p(2:end-1) = M \ (alpha * pold(2:end-1) + delta_t * F - 0.5 * lambda * D_k*pold(2:end-1));
            end
            %enforce 0 Neumann BC
            if (strcmp(BC_type,'Neu'))
                p(end) = p(end-1);
            end
            F = F0; %set back at time step otherwise accumulating BC terms
            p(1) = bdry_left(t(i+1));
            if (strcmp(BC_type, 'Dir'))
                %set Dirichlet RHS
                p(end) = bdry_right(t(i+1));
            end
        end
    end
    k = (p.^k_pow)';
    %k's are 1 ts behind
    k_exact(x_cent < t(end)) = k_pow*(t(end)-x_cent(x_cent < t(end)));
    k_exact(x_cent >= t(end)) = eps; %see if fixes harmonic problem
    if (use_exact) %compute at final timestep for FE
        k = k_exact;
    end
    p_exact = (k_exact.^(1/k_pow))';
    p_t(end) = p(ind);
    p_t_exact(end) = p_exact(ind);
    
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
    
    filename = sprintf('~/Desktop/parabolic_conv_and_plots/numerical_k/expandKeqtn/%s/h_sqr/arithmetic',timeScheme);
    h1 = figure(1);
    hold on;
    plot(x_cent,p,x_cent,p_exact)
    legend('p','p_{exact}')
    title(sprintf('Plot of numerical pressure versus exact pressure'))
    xlabel('position x')
    ylabel('pressure p')
    %saveas(h1, strcat(filename,sprintf('/nx%i.jpg',N)));
    
    h2 = figure(2);
    hold on;
    plot(x_cent,p,'*',x_cent,p_exact,'*')
    legend('p','p_{exact}')
    title(sprintf('Plot of numerical pressure versus exact pressure'))
    xlabel('position x')
    ylabel('pressure p')
    %saveas(h2, strcat(filename, sprintf('/nx%i_star.jpg', N)));
    
    h3 = figure(3);
    plot(t,p_t, t, p_t_exact)
    xlabel('time t')
    ylabel(sprintf('pressure at position x = %.1f', x_coord))
    title(sprintf('Plot of pressure versus time at x = %.1f', x_coord))
    legend('numerical pressure','exact pressure')  
    %saveas(h3, strcat(filename, sprintf('/time%.1f_nx%i.jpg',x_coord,N))); 
end
