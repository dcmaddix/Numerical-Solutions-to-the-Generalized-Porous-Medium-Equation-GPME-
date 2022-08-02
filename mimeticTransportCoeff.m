%Mimetic 2016 paper implementing Scheme I for various ki functions using
%basic time implementaions: Solving alpha*dp/dt -  d/dx(kdp/dx) = F
function [p p_t t pos_5] = mimeticTransportCoeff(k_i, N, k_pow, bdry, ghostNodes,...
                                           F0, p, alpha, BC_type, timeScheme,...
                                           nt, dt, dx, use_exact, errorType,...
                                           x_coord, theta, x_cent, avg, ...
                                           eps_pxx, eps_px, eps_pxxx, eigen, ...
                                           eps, ramp, slowDiff)
%% INPUTS:
    %k_i - function to approx cell-centered diffusion coeff k_tilde at node
           %point usually depends on the arithmetic and harmonic cell-centered
           %averages
    
    %N - number of interior grid points
    
    %k_pow - Assuming k = p^(k_pow)
         
    %bdry - function handle for left bdry function which grows in time and
            %right bdry value which is either Neumann or Dir BC
            
    %ghostNodes - left and right ghost node for wider higher order stencil
            
    %F0 - forcing function function handle that depends on time
    
    %p - Initial Condition
    
    %alpha - compressibility constant alphadp/dt
    
    %BC_type - Dir or Neu for right bdry assuming Dirichlet on left bdry
    
    %timeScheme - ForwardEuler, BackwardEuler, CrankNicolson, RK2_TVD
    
    %nt - Final time
    
    %dt - Timestep
    
    %dx - Spatial Step
    
    %use_exact - Boolean Flag as whether or not to use exact k for
                %unlocking harmonic
    
    %errorType - l2norm, maxnorm
    
    %x_coord - x position to store pressure at every time for pressure
                %versus time plot to check for oscillations
     
    %theta - theta method for Crank Nicolson theta*p^n + (1-theta)*p^(n+1)
             %better behaved and less spurious oscillations than theta =
             %0.5 for CN
    
    %x_cent - position vector
    
    %avg - averaging of permeabilites either arith or fv harm
    
    %eps_pxx - parameter for turning on or off anti-diffusion
    
    %eps_px - parameter for turning on or off advective error term
    
    %eigen - boolean if turned on compute the eigenvalues for linear
             %stability analysis
%% OUTPUTS:
    %p - pressure for plotting
    
    %p_t - pressure at specific x-coordinate
    
    %t- time vector
     
%% Initialization         
    %Plot pressure at this x-coordinate for each time to get pressure vs
    %time plot
    ind = find(x_cent < x_coord+1e-6 & x_cent > x_coord - 1e-6, 1, 'first');
    if (isempty(ind)) %make sure xcoord exists
        fprintf('Ind warning\n')
        return
    end
    
    %Define arithmetic averaging
    k_A = @(k1,k2) 0.5 * (k1 + k2);
    %Compute time vector
    t = 0: dt : nt;
    p_t = zeros(1,length(t)); %time vector where initial condition is given at t=0
    pos_5 = zeros(1,length(t)- 1);
    %compute constant lambda
    lambda = dt / dx^2;
    TV = zeros(1,length(t));
    trun = zeros(1,length(t));
    
    %max/min eigenvalue storage
    if (eigen)
        eval_t = zeros(2,length(t));
    end
    
    %define for exact solution
    if (use_exact)
        eps = 1e-9;
        k_exact = zeros(N+1,1)';
        p_exact = zeros(N+1,1);
        p_t_exact = zeros(1,length(t)); %plotting p at specific x coord versus time
    end
    %% Vector of interior indices
    i = 2:length(x_cent) - 1;
    %Time loop: already have it at first time so less iteration
    %special Jakolein case
    if (k_pow == -1)
        p_star = 0.5;
        upper = p_star + eps; %band
        lower = p_star - eps;
        k_upper = 1.0;
        k_lower = 0.01;
        %linear case compute slope of line
        if (strcmp(ramp, 'linear'))
            m = (k_upper - k_lower) / (2 * eps); %upper - lower (p_star cancels)
            b = k_upper - m * upper;
        elseif(strcmp(ramp, 'arctan'))
            C = (k_upper + k_lower) / 2.0;
            A = (k_upper - C) * (4 / pi);
            B = 1 / eps;
        elseif(strcmp(ramp, 'sine'))
            C = (k_upper + k_lower) / 2.0;
            A = k_upper - C;
            B = pi / (2*eps);
        end
    end
    p_x = zeros(length(i), 1);
    p_xx = zeros(length(i), 1);
    p_xxx = zeros(length(i), 1);
    p_xxxx = zeros(length(i), 1);
    p_xx2 = zeros(length(i), 1);
    for n = 1:length(t)-1
        %write simple finite volume solver and separate BC    
        %this makes it explicit change to Newton Raphson for fully-implicit
        F = F0;
        %store at current time
        p_t(n) = p(ind);
        TV(n) = sum(abs(p(2:end) - p(1:end-1)));
        %update at each timestep depends on new p depends on which
        %timestepping scheme
        if (k_pow ~= -1)
            k = (p.^k_pow)';
        else %Jakolein discontinous nonlinear k case
            if (eps ~= 0)
                if (strcmp(ramp, 'linear'))
                    k = m * p' + b;
                elseif(strcmp(ramp, 'arctan'))
                    k = A * atan(B * (p - 0.5))' + C;
                elseif(strcmp(ramp, 'sine'))
                    k = A * sin(B * (p - 0.5) )' + C;
                end
            end
            k(p <= lower) = k_lower;
            k(p > upper) = k_upper;
            %Only define in ramp outside cap to values
        end
        if (slowDiff)
            k = exp(- 1 ./ p)';
        end
        %% Exact parameters
        if (use_exact)
            if (strcmp(timeScheme, 'BackwardEuler')) %need to use future k values- compute at next t
                ind_exact = x_cent < t(n+1);
                ind_comp = x_cent >= t(n+1);
                t_coord = t(n+1);
                %using explicit TVD RK Scheme
            elseif (or(strcmp(timeScheme, 'ForwardEuler'), strcmp(timeScheme, 'RK2_TVD')))
                ind_exact = x_cent < t(n);
                ind_comp = x_cent >= t(n);
                t_coord = t(n);
            elseif(strcmp(timeScheme, 'CrankNicolson'))
                t_coord = theta *t(n) + (1-theta)*t(n+1);
                ind_exact = x_cent < t_coord;
                ind_comp = x_cent >= t_coord;
            end
            k(ind_exact) = k_pow * (t_coord - x_cent(ind_exact));
            k(ind_comp) = 1e-3^(k_pow); %see if fixes harmonic problem
            
            p_exact(x_cent < t(n)) = (k_pow*(t(n)-x_cent(x_cent < t(n))) ).^(1/k_pow);
            p_exact(x_cent >= t(n)) = eps^(1/k_pow); %see if fixes harmonic problem
            p_t_exact(n) = p_exact(ind); %store at every timestep
        end
        %do taylor linearization of matrix-can't change for BE/CN
        %important how update the nonlinear k in time
        %% Compute truncation error term
        %% Averaging of k and matrix
        %define Ti which are the tridiagonal entries of the matrix. Formula
        %simplifies to below for uniform spacing k
        if (strcmp(avg, 'upwind'))
            %k is increasing with p
            k_plus = k(i)'; %equivalent to max(k(i), k(i+1))' since p(i) > p(i+1) k(p(i)) > k(p(i+1))
            k_minus = k(i-1)'; %equivalent to max(k(i-1), k(i))' since p(i-1) > p(i) k(p(i-1)) > k(p(i))
        else
            k_plus = (k_i(k(i), k(i+1)).^2 ./ k_A(k(i), k(i+1)))';
            k_minus = (k_i(k(i-1), k(i)).^2 ./ k_A(k(i-1), k(i)))';
        end
        if(~strcmp(timeScheme, 'ForwardEuler'))
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

            %Not fully implicit since D_k depends on p at the old timestep
            %rather than the new timestep
            if (strcmp(timeScheme, 'ForwardEuler')) 
                M = speye(N-1,N-1) - (lambda/alpha) * D_k;
            elseif (strcmp(timeScheme, 'BackwardEuler')) 
                M = speye(N-1,N-1) + (lambda/alpha)*D_k;
            elseif (strcmp(timeScheme, 'CrankNicolson')) %Crank Nicolson
                M  = speye(N-1,N-1) + (1-theta)*(lambda/alpha)*D_k;
            end
        end
        
        %Add Dirichlet BC contribution to RHS-fix boundary by dividing by
        %h^2 as done with matrix term D_k
        %fine to use old k_minus and p for F.E but B.E. needs to be added
        %at new timestep
        
        %% compute boundary terms
        bdry_term = zeros(2,1);
        if (or(strcmp(timeScheme, 'ForwardEuler'), strcmp(timeScheme, 'RK2_TVD')))
            bdry_term(1) = p(1); %explicit use at current timestep
            bdry_term(2) = p(end);
        elseif(strcmp(timeScheme, 'BackwardEuler'))
            bdry_term(1) = bdry{1}(t(n+1)); %p(1) at next timestep
            bdry_term(2) = bdry{2}(t(n+1));
        elseif(strcmp(timeScheme, 'CrankNicolson'))
            %average of bdry value at current and new timestep-one from
            %each matrix to add to first and last components of F
            bdry_term(1) = (theta*p(1) + (1-theta)*bdry{1}(t(n+1)));
            bdry_term(2) = (theta*p(end) + (1-theta)*bdry{2}(t(n+1))) ;
        end
        
        F(1) = F(1) + k_minus(1) * bdry_term(1) / dx^2; %at next timestep for Backward Euler bdry terms
        
        if (strcmp(BC_type, 'Dir')) %verified sign for forward euer
            F(end) = F(end) + k_plus(end)*bdry_term(2) / dx^2; %Dirichlet RHS correct for bdry term on RHS vector
        end
        %% fix oscillations in harmonic case positive for arithmetic case
        %use 4th order to cancel
        ghostNodes_t = [ghostNodes{1}(t(n)) ghostNodes{2}(t(n)) ...
                            ghostNodes{1}(t(n)) ghostNodes{2}(t(n))];
        if (strcmp(avg, 'fv harm') ||k_pow == -1) %fourth order to cancel
            %Fourth Order First Derivative
            if (slowDiff)
                dp_dx = (p(i) - p(i-1)) / dx;
            else
                dp_dx = central_4thOrder(p, ghostNodes_t, dx);
                %dp_dx = (p(i+1) - p(i-1)) / (2*dx);
                %dp_dx = (p(i) - p(i - 1)) / dx;
            end
            %Fourth Order Second Derivative
            d2p_dx2 = central_2ndDeriv_4thOrder(p, ghostNodes_t, dx);
            %Fourth Order Third Derivative
            d3p_dx3 = central_3rdDeriv_4thOrder(p, ghostNodes_t, dx);
            %Fourht Order Fourth Derivative
            d4p_dx4 = central_4thDeriv_4thOrder(p, ghostNodes_t, dx);
        elseif (strcmp(avg, 'arith') || strcmp(avg, 'upwind'))
            dp_dx = (p(i+1) - p(i-1)) / (2 * dx); %central
            dp_dx = central_4thOrder(p, ghostNodes_t, dx);
            d2p_dx2 = (p(i-1) - 2 * p(i) + p(i+1)) / dx^2; %to add create oscillations use second order
        end
        %% k = p^n
        if (k_pow ~= -1) %k = p^n case
            p_xx = dx^2 * dp_dx.^2 .* d2p_dx2;
            p_xx = dx^2 * dp_dx.^2 .* (p(i+1) - 2*p(i) + p(i-1)) / (dx^2);
            p_x = dx^2 * dp_dx.^4;
        end
        %have negative signs in truncation error from semi-discrete: Cancel 
        %to remove oscillations
        if (strcmp(avg, 'fv harm'))
            if (k_pow == 1)
                p_xx = 3 / 4 * p(i).^(-1) .* p_xx;
                p_x = -1/ 4 * p(i).^(-2) .* p_x;
            elseif (k_pow == 2)
                p_xx = 3 / 2 * p_xx;
                p_x = 0 * p_x;
            elseif (k_pow == 3)
                p_xx = 9 / 4 * p(i) .* p_xx;
                p_x = 5 / 4 * p_x;
            elseif (slowDiff) %confirmed correct derivatives!
                k_p = k(i)' .* p(i).^(-2);
                k_pp = k(i)' .* (p(i).^(-4) - 2*p(i).^(-3));
                k_ppp = k(i)'.* (p(i).^(-6) - 6*p(i).^(-5) +6*p(i).^(-4));
                p_xx = -3/4 * (k_pp - k_p.^2 .* k(i).^(-1)') .*p_xx;
                p_x = (-1/6 * k_ppp + 1/4 * (2*k_p .*k_pp .* k(i)'.^(-1)...
                      - k_p.^3 .*k(i)'.^(-2))) .* p_x;           
            end 
        %Have positive sign in truncation error from semi-discrete: Cancel
        %to create oscillations
        elseif(strcmp(avg, 'arith'))
            if (k_pow == 1)
                p_xx = 0 * p_xx;
                p_x = 0 * p_x;
            elseif (k_pow == 2)
                %Add antidiffusion 3/2 in term and not strong enough for
                %oscillations: This is negative quantity that harmonic adds
                %on
                p_xx = -3 * p_xx;
                p_x = 0 * p_x;
            elseif (k_pow == 3)
                %Canceling from arithmetic trunc error
                p_xx = -9 / 2 * p(i) .* p_xx;
                p_x  = -p_x;
            end
        end   
        %Compute coefficient of pxx = (k + (A+B+deltaB)) negative at
        %specific gridpt for harmonic
        if (strcmp(avg, 'fv harm')) %applicable to k = p^n case
            if (k_pow == 1)
                trun = k(i)' - (dp_dx.^2 .* (7/2 * dt + 3/4 * dx^2 * p(i).^(-1)));
            elseif (k_pow == 2)
                trun = k(i)' - (dp_dx.^2 .* (21 * dt * p(i).^2 + 3/2 * dx^2) );
            elseif (k_pow == 3)
                trun = k(i)' - (dp_dx.^2 .* (105/2 * dt * p(i).^4 + 9/4 * dx^2 * p(i)));
            end
        elseif(strcmp(avg,'arith'))
             if (k_pow == 1)
                trun = k(i)' - (dp_dx.^2 .* 7/2 * dt); 
             elseif (k_pow == 2)
                trun = k(i)' + dp_dx.^2 .* (-21 * dt * p(i).^2 + 3/2 * dx^2);
             elseif (k_pow == 3)
                trun = k(i)' + dp_dx.^2 .* (-105/2 * dt * p(i).^4 + 9/2 * dx^2 * p(i));
             end
        elseif (~strcmp(avg, 'upwind'))
            fprintf('Undefined Averaging Type\n')
            return;
        end
        if (any(trun < 0))
            min(trun);
        end
        
        %% Nonlinear discontinuous permeability case
        %Express in terms of p derivatives of k for general nonlinear k(p)
        if(k_pow == -1)
            k_p = zeros(length(i), 1);
            k_pp = zeros(length(i), 1);
            k_ppp = zeros(length(i), 1);
            if (eps ~= 0)
                if (strcmp(ramp, 'linear'))
                    k_p = m * ones(length(i), 1);
                elseif(strcmp(ramp, 'sine'))
                    k_p = A * B * cos(B * (p(i) - 0.5));
                    k_pp = -A * B^2 * sin(B * (p(i) - 0.5));
                    k_ppp = -A * B^3 * cos(B * (p(i) - 0.5));
                elseif (strcmp(ramp, 'arctan')) %cancel to avoid large inv eps
                    inv = (eps^2 + (p(i) - 0.5).^2).^(-1);
                    k_p = A * eps * inv;
                    k_pp = -2 * A * eps * (p(i) - 0.5) .* inv.^2;
                    k_ppp = 8 * A * eps * (p(i) - 0.5).^2 .* inv.^3 ...
                            - 2 * A * eps *inv.^2;
                end
            end
            const_ind = p(i) < lower | p(i) > upper; %constant outside range
            k_p(const_ind) = 0.0; %not working to for p_xx
            k_pp(const_ind) = 0.0; %not working to for p_xx
            k_ppp(const_ind) = 0.0;%otherwise travel too far only apply locally
            %Use to calc deriv in modified eqtn- cancelling everything
            %Anti-diffusion in this term can lead to oscillations
            p_xx = -dp_dx.^2 .* d2p_dx2 .* (-7 * dt / 2 * (k_p.^2 + k(i)' .* k_pp)...
                                           + 3 * dx^2 / 4 * k_pp);
            %Second derivative squared term can impact the frequencies
            p_xx2 = -d2p_dx2.^2 .* k_p .* (-2 * k(i)' * dt + dx^2 / 4);

            %Dissipative Third Derivative Term
            p_xxx = -d3p_dx3 .* dp_dx .* k_p .* (-3 * k(i)' * dt + dx^2 / 3);

            %Advective First Derivative to the 4th term
            p_x = -dp_dx.^4 .*(-dt / 2 * (3 * k_p .* k_pp + k(i)' .* k_ppp)...
                              + dx^2 / 6 * k_ppp);
            %Fourth Derivative term
            p_xxxx = -d4p_dx4 .* k(i)' .* (-k(i)' * dt/ 2 + dx^2 /12 );
            
            if (strcmp(avg, 'fv harm'))
                %Harmonic adds on additional antidiffusive term to cancel form
                %semidiscrete equation
                %p_xx = p_xx + (3 * dx ^ 2 / 4 * k_p.^2 .* k(i)'.^(-1)...
                          %  .* dp_dx.^2 .* d2p_dx2);
                %Helps correct propagation speed and unlock harmonic
                %p_x = p_x + dx^2 / 4 * (2 * k_p .* k_pp .* k(i)'.^(-1) ...
                                   %   - k_p.^3 .*k(i)'.^(-2)) .* dp_dx.^4;
            end
        end
        %% Additional explicit update in time for p: k = p^n only nonzero
        % contributions from the first two term
         %% Update pressure step
        %solve (N-1)x(N-1) linear system for interior pressures
        if (strcmp(timeScheme, 'ForwardEuler'))
            %need to add artificial diffusion to cancel dx^2 term!
            %p(i) = M*p(i) + dt/alpha * F;
            %p(i) = p(i) + (1 / alpha)*(dt * F - lambda * D_k * p(i));
            %Don't need to add to F if doing this case below
            F = F0;
            p(i) = p(i) + (1 / alpha) * (dt * F + lambda * (k_minus.*p(i-1) ...
                - (k_minus + k_plus) .* p(i) + k_plus .* p(i+1)));
        elseif (strcmp(timeScheme, 'RK2_TVD'))
            %p^(1) = p^n+ deltat_t * D_k (FE step) m = 1 first stage- F is
            %containing bdry info
            p1 = p;
            %p1(i) = p(i) + (1 / alpha) * (dt * F + lambda * (k_minus.*p(i-1) ...
                 %- (k_minus + k_plus) .* p(i) + k_plus .* p(i+1)));
            p1(i) = p(i) + (1/alpha) * (dt * F - lambda * D_k * p(i));
            p1(1) = bdry{1}(t(n+1));
            p1(end) = bdry{2}(t(n+1));
            %p(i) = 0.5 * (p(i) + p1(i)) + 0.5 * (1 / alpha) * (dt * F + lambda * (k_minus.*p1(i-1) ...
                 %- (k_minus + k_plus) .* p1(i) + k_plus .* p1(i+1)));
            F = F0;
            F(1) = F(1) + k_minus(1)*p1(1) / dx^2;
            if (strcmp(BC_type, 'Dir')) %verified sign for forward euer
                F(end) = F(end) + k_plus(end)*p1(end) / dx^2; %Dirichlet RHS correct for bdry term on RHS vector
            end
            %do second stage and add on first and last component to F from
            %tridiagonal mat vec with
            p(i) = 0.5*(p(i) + p1(i)) + 0.5/alpha*(dt * F - lambda*D_k*p1(i));
        elseif(strcmp(timeScheme, 'BackwardEuler'))
            %constant F in time otherwise 
            p(i) = M \ (p(i) + dt/alpha * F); 
        elseif(strcmp(timeScheme,'CrankNicolson'))
            %F is constant in time if had time dependency we want
            %1/2(f(xi,tn) + f(xi,tn+1));
            p(i) = M \ (p(i) - theta * lambda * D_k*p(i) + dt/alpha * F);
            %p(i) = M * p(i) + M_F \ ((dt / alpha) * F); 
        end
        p(i) = p(i) + dt / alpha * (eps_pxx*p_xx + eps_px * p_x + ...
                                   eps_pxxx*p_xxx + 0.0*p_xxxx + 0.0*p_xx2);
        %% Interpolate to find value at 0.5
        shift = 0.1;
        p_int = p(shift < p  & p < max(p) - shift);
        x_int = x_cent(shift < p  & p < max(p) - shift);
        %find pos x at corresponding pressure points reverse interp
        if (length(p_int) >= 2)
            pos_5(n) = interp1(p_int, x_int, 0.5);
        end
        %% Update boundaries
        %enforce 0 Neumann BC
        if (strcmp(BC_type,'Neu'))
            p(end) = p(end-1);
        end
        p(1) = bdry{1}(t(n+1));
        if (strcmp(BC_type, 'Dir'))
            %set Dirichlet RHS
            p(end) = bdry{2}(t(n+1));
        end
        %% Eigenvalue Stability analsyis
        %expensive to invert stability analysis ok for these two since smae
        %oscillation results reglardless of tiem scheme
        if (eigen)
            if (or(strcmp(timeScheme, 'ForwardEuler'),strcmp(timeScheme, 'BackwardEuler')))
                evals = eig(M);
                if(strcmp(timeScheme, 'BackwardEuler'))
                    %For implicit need to invert eignvalues to calculate for inversematrix
                    evals = 1 ./ evals; %p^(n+1) = (I+lambdaD_k)^(-1)p^n
                end
                eval_t(1,n) = max(evals);
                eval_t(2,n) = min(evals);
                if (~isreal(evals))
                    fprintf('Warning: Complex eigenvalues!')
                    return;
                end
                if (max(abs(evals)) >= 1) %for stability make sure less than 1 in abs value
                    max(evals);
                    min(evals);
                    fprintf('Warning: Eigenvalues larger in magnitude than 1!')
                    return;
                end
            end
        end
    end
    %% Final trun error
    %copy same boundary node to both ghost nodes on either side using Dir
    %BC alteratnes 1st left, 1st right, 2nd lef,t 2nd right in function
    ghostNodes_t = [ghostNodes{1}(t(n)) ghostNodes{2}(t(n)) ...
                    ghostNodes{1}(t(n)) ghostNodes{2}(t(n))];
    dp_dx = central_4thOrder(p, ghostNodes_t, dx);
    if (strcmp(avg, 'fv harm'))
        if (k_pow == 1)
            trun = k(i)' - (dp_dx.^2 .* (7/2 * dt + 3/4 * dx^2 * p(i).^(-1)) );
        elseif (k_pow == 2)
            trun = k(i)' - (dp_dx.^2 .* (21 * dt * p(i).^2 + 3/2 * dx^2) );
        elseif (k_pow == 3)
            trun = k(i)' - (dp_dx.^2 .* (105/2 * dt * p(i).^4 + 9/4 * dx^2 * p(i)));
        end
    elseif(strcmp(avg,'arith'))
         if (k_pow == 1)
            trun = k(i)' - (dp_dx.^2 .* 7/2 * dt);
         elseif (k_pow == 2)
            trun = k(i)' + dp_dx.^2 .* (-21 * dt * p(i).^2 + 3/2 * dx^2);
         elseif (k_pow == 3)
            trun = k(i)' + dp_dx.^2 .* (-105/2 * dt * p(i).^4 + 9/2 * dx^2 * p(i));
         end
    end
    if (k_pow ~= -1)
        figure(100)
%        plot(x_cent(2:end-1), trun, '-o')
        hold on
        ylabel('k + p_{xx} truncation coeff.')
        xlabel('x')
        title('Plot of coefficient of p _{xx}')
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca, 'FontSize',16,'FontWeight','bold');
        save('trun_harm_m1Nx100','trun', 'x_cent')
    end
    %% Final k and p
    %Compute at final time
    if (k_pow ~= -1)
        k = (p.^k_pow)';
    end
    figure(20)
    plot(p,k, '-o')
    title('k versus p')
    xlabel('p')
    ylabel('k')
    hold on;
    TV(end) = sum(abs(p(2:end) - p(1:end-1)));
    %k's are 1 ts behind
    if (use_exact)
        k_exact(x_cent < t(end)) = k_pow*(t(end)-x_cent(x_cent < t(end)));
        k_exact(x_cent >= t(end)) = eps; %see if fixes harmonic problem
        k = k_exact;
        p_exact = (k_exact.^(1/k_pow))';
        p_t_exact(end) = p_exact(ind);
    end
    p_t(end) = p(ind);
    figure(11)
    plot(t, TV)
    figure(12)
    plot(x_cent, p)
end
