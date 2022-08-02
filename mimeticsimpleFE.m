%Mimetic 2016 paper implementing Scheme I for various ki functions using
%basic time implementaions: Solving alpha*dp/dt -  d/dx(kdp/dx) = F
function [p p_t p_exact p_t_exact t alpha_star xi_plot] = mimeticsimpleFE(k_i, k_pow, ...
                                                                  bdry, p, alpha, ...
                                                                  nt, dt, dx, x_coord, ...
                                                                  x_cent, avg, eps, ...
                                                                  ramp, F, k_lower,...
                                                                  modEqtn, globalAtan, noAvg,...
                                                                  fluxEqual, t0, useOldv, ...
                                                                  artDiffTime, artDiffSpace, ...
                                                                  newUpdateStep, weightedTimeAvg, ...
                                                                  weightedTimeAvgOld, IC_evolve, ...
                                                                  exact, levelset, linear)
%% INPUTS:
    %k_i - function to approx cell-centered diffusion coeff k_tilde at node
           %point usually depends on the arithmetic and harmonic cell-centered
           %averages
    
    %k_pow - Assuming k = p^(k_pow)
         
    %bdry - function handle for left bdry function which grows in time and
            %right bdry value which is either Neumann or Dir BC
    
    %p - Initial Condition
    
    %alpha - compressibility constant alphadp/dt
    
    %nt - Final time
    
    %dt - Timestep
    
    %dx - Spatial Step
    
    %x_coord - x position to store pressure at every time for pressure
                %versus time plot to check for oscillations
     
    %x_cent - position vector
    
    %avg - averaging of permeabilites either arith or fv harm
    
    %eps - Width of gap region in discontinuity in nonlinear permeability
          %function
    %ramp - Type of ramp to add in: linear, sine or arctan
    
    %F - RHS forcing function
    
    %k_lower - minimum piecewise constant k value assuming k_upper = 1
    
    %modEqtn - flag as to whether to add on truncation error terms
    
    %globalAtan - flag as to whether use global arctan or arctan only
                 %defined in ramp
%% OUTPUTS:
    %p - numerical pressure for plotting
    
    %p_t - numerical pressure at specific x-coordinate as stored at every
            %timestep
    
    %p_exact - exact solution at final timestep
    
    %p_t_exact - exact pressure at specific x-coordinate as stored at every
                %timestep
    
    %t- time vector
    
    %alpha_star - exact shock location is alpha_star*sqrt(t)
     
%% Initialization         
    %Plot pressure at this x-coordinate for each time to get pressure vs
    %time plot
    cell = find(x_cent < x_coord + 1e-6 & x_cent > x_coord - 1e-6, 1, 'first');
    probe = cell; %for shfiting and plotting next grid cell
    if (isempty(cell)) %make sure xcoord exists
        fprintf('Ind warning\n')
        return
    end
    %Define arithmetic averaging
    k_A = @(k1,k2) 0.5 * (k1 + k2);
    t_shock = [];
    %Compute time vector
    %start at later time
    t = t0 : dt : t0 + nt;
    %start with old intial condition
    p_t = zeros(length(x_cent),length(t)); %time vector where initial condition is given at t=0
    p_t = zeros(1,length(t));
    tvd = zeros(1,length(t));
    k_t = tvd;
    p_t_exact = p_t;
    xi_plot = zeros(1,length(t));
    speed_plot = zeros(1,length(t));
    xi_error = zeros(1,length(t));
    count = 0;
    %Initialize exact solution to 0 doesn't exist in p^n case
    p_exact = zeros(length(x_cent), 1);
    alpha_star = 0.0;
    %compute constant lambda
    lambda = dt / dx^2;
    %% Vector of interior indices
    i = 2 : length(x_cent) - 1;
    %Initialize to 0
    k_minus_exact = zeros(length(i),1);
    num_speed = zeros(length(t),1);
    k_plus_exact = k_minus_exact;
    pold = zeros(length(p), 3);
    %Time loop: already have it at first time so less iteration
    %special Jakolein case
    if (k_pow == -1)
        p_star = 0.5;
        upper = p_star + eps; %band
        lower = p_star - eps;
        k_upper = 1.0;
        %linear case compute slope of line
        if (strcmp(ramp, 'linear'))
            m = (k_upper - k_lower) / (2 * eps); %upper - lower (p_star cancels)
            b = k_upper - m * upper;
        elseif(strcmp(ramp, 'arctan'))
            B = 1 / eps;
            %global arctan
            if (globalAtan)
                A = (k_upper - k_lower) / (2*atan(0.5*B));
                C = 1 - A*atan(0.5 * B);
            else %ramp arctan
                C = (k_upper + k_lower) / 2.0;
                A = (k_upper - C) * (4 / pi);
            end
        elseif(strcmp(ramp, 'sine'))
            C = (k_upper + k_lower) / 2.0;
            A = k_upper - C;
            B = pi / (2*eps);
        end
        %% exact-solution solve nonlinear for exact solution
        if (k_lower < 1e-6)
            lim = true; %use limiting solution where RHS solution is 0
        else 
            lim = false; %use solution for exact eps solving for
        end
        [alpha_star A_exact B_exact] = exact_sol_speed(lim, p_star, k_lower);
    end
    %Plot flux at specific x coord within same cell at every timestep
    flux_in = zeros(length(t), 1);
    flux_out = zeros(length(t), 1);
    flux_in_exact = zeros(length(t), 1);
    flux_out_exact = zeros(length(t), 1);
    k_eff = zeros(length(t), 1);
    k_eff_exact = zeros(length(t), 1);
    dx_star_t = zeros(length(t), 1);
    dx_star_exact = zeros(length(t), 1);
    for n = 1:length(t)-1
        %tvd(n) = sum(abs(p(2:end) - p(1:end-1)));
        %store at current time
        p_t(n) = p(probe);
        %p_t(:,n) = p;
        if (k_pow == -1)
              [p_t_exact(n) p_exact] = exact_solution(t(n), x_cent, cell, alpha_star,...
                                              A_exact, B_exact, k_lower, lim);
              tvd(n) = sum(abs(p_exact(2:end) - p_exact(1:end-1)));
              if (n == 1)
                  p_exact = p;
              end
        end
        %update at each timestep depends on new p depends on which
        %timestepping scheme
        if (k_pow ~= -1)
            k = (p.^k_pow)';
        else %Jakolein discontinous nonlinear k case
            if (eps ~= 0)
                if (strcmp(ramp, 'linear'))
                    k = m * p' + b;
                elseif(strcmp(ramp, 'arctan'))
                    k = A * atan(B * (p - 0.5))' + C; %global arctan
                elseif(strcmp(ramp, 'sine'))
                    k = A * sin(B * (p - 0.5) )' + C;
                end
            end
            if (~globalAtan || eps == 0)
                k(p <= lower) = k_lower;
                k(p > upper) = k_upper;
            end
        end
        %Added factor since too fast
        vel_exact = 0.1 * alpha_star * 0.5 * t.^(-0.5);
        %% Averaging of k and matrix
        %define Ti which are the tridiagonal entries of the matrix. Formula
        %simplifies to below for uniform spacing k
        if (strcmp(avg, 'upwind'))
            %k is increasing with p
            k_plus = k(i)'; %equivalent to max(k(i), k(i+1))' since p(i) > p(i+1) k(p(i)) > k(p(i+1))
            k_minus = k(i-1)'; %equivalent to max(k(i-1), k(i))' since p(i-1) > p(i) k(p(i-1)) > k(p(i))
        elseif(k_pow == -1 && strcmp(avg, 'int')) %This case is only defined for discontinuous nonlinear k case
            %k_i-1/2 averaged on cell face
            k_minus = integral_average(p, i - 1, p_star, k_lower, k_upper);
            %k_i+1/2 averaged on right cell face
            k_plus = [k_minus(2:end); ...
                      integral_average(p, length(p) - 1, p_star, k_lower, k_upper)];
            %compute on last cell p is length N+1 so on(N,N+1)
            ind_shock = find(p(i - 1) >= p_star & p_star > p(i));
            if(ind_shock == cell - 1) %cell is i+1
                t_shock = [t_shock t(n)];
            end
            k_t(n) = k_plus(probe-2);
        elseif(k_pow == -1 && strcmp(avg, 'k_eff'))
            %k_i-1/2_exact
            if (fluxEqual || newUpdateStep)
                [k_minus k_plot ratio] = computek_eff(p, i - 1, p_star, ...
                                                        k_lower, k_upper);
                k_plus = [k_minus(2:end); ...
                      computek_eff(p, length(p) - 1, p_star, k_lower, k_upper)];
            elseif(useOldv) %locks
                %k_i-1/2 averaged on cell face
                if (n == 1)
                    vold = Inf; %take the computed on the first timestep
                else
                    vold = v;
                end
                [k_minus k_plot v] = computek_eff_v(p, i - 1, p_star, k_lower, ...
                                                    k_upper, vold, dx);
                %k_i+1/2 averaged on right cell face
                k_plus = [k_minus(2:end); ...
                            computek_eff_v(p, length(p) - 1, p_star, k_lower, ...
                                                    k_upper, vold, dx)];            
            else %works using exact v
                [k_minus k_plot ratio] = computek_eff_exactv...
                                         (p, i - 1, p_star, k_lower, k_upper, ...
                                         vel_exact(n), dx);
                k_plus = [k_minus(2:end); ...
                      computek_eff_exactv(p, length(p) - 1, p_star, k_lower, ...
                                         k_upper, vel_exact(n), dx)];
            end
            %Compute exact using the equations for dx* and the exact
            %solution
            [k_minus_exact k_plot_exact exact_ratio] = computek_eff...
                                                        (p_exact, i - 1, p_star, ...
                                                        k_lower, k_upper);
            k_plus_exact = [k_minus_exact(2:end); ...
                       computek_eff(p_exact, length(p) - 1, p_star, k_lower, k_upper)];
            %compute on last cell p is length N+1 so on(N,N+1)  
            k_eff(n) = k_plot;
            k_eff_exact(n) = k_plot_exact;
            if (~useOldv) %not defined when use oldV
                dx_star_t(n) = ratio;
            end
            dx_star_exact(n) = exact_ratio;
        else
            k_minus = (k_i(k(i-1), k(i)).^2 ./ k_A(k(i-1), k(i)))';
            %Find all indices where 0 implies k(i) = k(i-1) = 0
            %Both arith and harm should give 0
            k_minus(k_A(k(i), k(i-1)) == 0) = 0.0; 
            k_plus = [k_minus(2:end); ...
                     (k_i(k(end-1), k(end)).^2 ./ k_A(k(end-1), k(end)))];
            %Check last component
            if (k_A(k(end-1), k(end)) == 0.0)
                k_plus(end) = 0.0;
            end
            k_t(n) = k_minus(probe-1); %k_minus(cell-1)
            ind_shock = find(p(i - 1) >= p_star & p_star > p(i));
            if(ind_shock == cell - 1)
                t_shock = [t_shock t(n)];
            end
        end
        %different for shock cell!
        flux_in(n) = k_minus(cell-1) .* (p(cell-1) - p(cell)) / dx;
        flux_out(n) = k_plus(cell-1) .* (p(cell) - p(cell + 1)) / dx;
        %Compute F_i^+
        flux_in(n) = k_plus(cell-1) .* (pold(cell) - pold(cell + 1)) / dx;
        %Compute F_{i+1}^-
        flux_out(n) = k_minus(cell) .* (pold(cell) - pold(cell + 1)) / dx;

        flux_in_exact(n) = k_minus_exact(cell-1) .* (p_exact(cell-1) - p_exact(cell)) / dx;
        flux_out_exact(n) = k_plus_exact(cell-1) .* (p_exact(cell) - p_exact(cell+1)) / dx;
        %Compute F_i^+
        flux_in_exact(n) = k_plus_exact(cell-1) .* (p_exact(cell) - p_exact(cell+1)) / dx;
        %Compute F_{i+1}^-
        flux_out_exact(n) = k_minus_exact(cell) .* (p_exact(cell) - p_exact(cell+1)) / dx;
        if (modEqtn) %global arctan case
            k_p = A * eps .* (eps^2 + (p(i) - 0.5).^2).^(-1);
            k_pp = -2 * A * eps .* (eps^2 + (p(i) - 0.5).^2).^(-2) .*(p(i) - 0.5);
            dp_dx = (p(i+1) - p(i)) / dx;
            d2p_dx2 = (p(i+1) - 2*p(i) + p(i-1)) / dx^2;
            p_xx_coeff = 3*dx^2 / 4 * k_p.^2 .*k(i)'.^(-1) .* dp_dx.^2;
            p_x_coeff = dx^2/4 *(2*k_p .* k_pp .*k(i)'.^(-1) -...
                        k_p.^3.*k(i)'.^(-2)).*dp_dx.^3;
        end
        
        %% Update pressure step
        %solve (N-1)x(N-1) linear system for interior pressures
        if (~noAvg)
                %old identify shock for exactv calcul
            if (k_pow == -1 && ~fluxEqual)
                ind_shock = find(p(i - 1) >= p_star & p_star > p(i));%Normally ind_shock + 1 but since
                %shock between i - 1 and i need to subtract 1
                %% Add artificial diffusion in space to shock cells
                %centered at shock cell i
                diff_coeff = 1e-1 * artDiffSpace;
                pold = p;
                if (ind_shock > 1)
                    d2p_dx2 = (p(ind_shock + 1) - 2 * p(ind_shock) + p(ind_shock - 1)) / dx^2;
                end
                if (~newUpdateStep)
                    p(ind_shock) = p(ind_shock) + (1 / alpha) * (dt * F(ind_shock - 1) + lambda * (k_minus(ind_shock - 1) * p(ind_shock - 1) ...
                        - (k_minus(ind_shock - 1) + k_plus(ind_shock - 1)) .* p(ind_shock) + k_plus(ind_shock - 1) .* p(ind_shock+1))) + dt * diff_coeff * d2p_dx2;
                    %Calculate next permeability in shock
                    k_minus_iplus1 = k_plus(ind_shock - 1) - vel_exact(n)*dx / 8;
                    p(ind_shock + 1) = pold(ind_shock + 1) + (1 / alpha) * (dt * F(ind_shock) + lambda * (k_minus(ind_shock) .* pold(ind_shock) ...
                                - (k_minus_iplus1 + k_plus(ind_shock)) .* pold(ind_shock + 1) + k_plus(ind_shock) .* pold(ind_shock+2))) + dt * diff_coeff * d2p_dx2;
                    j = setdiff(i, [ind_shock ind_shock + 1]);
                else
                    %% Testing perturbations
                    epsilon = 5e-3;
                    pow = rand(1);
                    if (pow < 0.5)
                        pow = 0;
                    else
                        pow = 1;
                    end
                    wave_num = (2 * pi) / (dx / 2);
                    shock_loc = alpha_star * t(n)^(0.5);
                    perturb = epsilon * sin(shock_loc * wave_num);
                    %% Adding pertubation calc
                    %shock_loc = alpha_star * t(n)^(0.5);% + (-1)^pow * epsilon * dx;
                    %shock_loc = shock_loc * (1 + perturb);
                    %shock_loc = shock_loc * (0.9);
                    tol = 2e-4; %make proportional to dx 2e-4 case that works for exact sol
                    if ((length(x_cent) - 1) >= 200) %Need to increase the tolerance as the grid size is decreased
                        tol = 1e-4; %other tolerance doesn't work for 200 grid points
                    end
                    skip_rightcell = false;
                    skip_leftcell = false;
                    %Both conditions could happen
                    %compare to numerical speed
                    fR = -(k_lower .* (p(ind_shock + 3) - p(ind_shock + 2)) / dx);
                    fL = -(k_upper .* (p(ind_shock) - p(ind_shock - 1)) / dx);
                    %Correct with p(ind_shock + 1) = 0 in denominator!
                    if (linear)
                        %can use ind_shock +2 but not ind_shock +1 in RH condtion
                        speed = (fR - fL) / (p(ind_shock+2)-p(ind_shock));
                    else
                        speed = (p(ind_shock - 1) - p(ind_shock)) / (dx * (p(ind_shock)));
                        speed = (log(p(ind_shock - 1)) - log(p(ind_shock)))/ dx;
                    end
                    if (n == 1)
                        %Do a bisection here to get close
                        if (linear)
                            xi = 0.25; %linear ic
                            temp = xi;
                            phi = x_cent - 0.25;
                        else
                           xi = shock_loc - x_cent(ind_shock);
                           temp = shock_loc(1);
                           %initial signed distance function
                           phi = x_cent - shock_loc(1);
                        end
                        xi_plot(n) = temp;
                    end
                    %Crossed gridcell reset for next timestep
                    if (exact)
                       xi = shock_loc - x_cent(ind_shock);
                    %This equation is causing the lag: when we use exatc
                    %perfect!
                    %Only crosses one cell at a time
                    elseif (n > 1 && (ind_shockold ~= ind_shock)) %new cell
                       fprintf('%i %i\n', ind_shockold, ind_shock)
                       %perfect results
                       %xi = mod(xi, dx);
                       count = count + 1; %number times change cell
                    elseif (n > 1)
                       %xi = xi + speed * dt;
                    end
                    %xi_plot(n) = xi + count * dx; %counteract mod operator
                    if (levelset)
                        phi(2:end) = phi(2:end) - dt*speed*(phi(2:end) - phi(1:end-1)) / dx;
                        xi = abs(phi(ind_shock)); %signed distance take abs
                        xi_plot(n+1) = interp1(phi, x_cent, 0); %shock position in ls
                    else
                        temp = temp + speed * dt; %never reset as grid cell
                        xi = temp - x_cent(ind_shock);
                        xi_plot(n+1) = temp;
                    end
                    %value
                    speed_plot(n) = speed;
                    %Store old ind_shock
                    ind_shockold = ind_shock;
                    %Check when to skip cell
                    if (dx - xi <= tol)
                       %affects i cell update that is not skipped
                       skip_rightcell = true;
                    end
                    if (xi <= tol)
                       skip_leftcell = true;
                    end
                    %how to reset when shock crossed grid cell
                    %% Other initial condition
                    if (~IC_evolve)
                       xi = xi + dx; 
                    end
                  
                    %use regular stencil define averages
                    j = setdiff(i, []);
                    %Define k_minus, with special evolving value at the
                    %shock for ind_shock + 1
                    %lambda = k_minus;
                    k_minus(1:ind_shock - 1) = k_upper;
                    %add in correct volumes
                    %lambda(1:ind_shock - 2) = dt / dx^2;
                    %lambda(ind_shock - 1) = dt / (dx*(dx+xi)/2);
                    %lambda(ind_shock) = dt / (dx * (dx - xi/2));
                    %lambda(ind_shock+1:end) = dt / dx^2;
                    k_minus(ind_shock) = k_lower * (dx / (dx - xi)) * (p_star - p(ind_shock +1)) ...
                                            / (p(ind_shock) - p(ind_shock + 1));
                    k_minus(ind_shock + 1:end) = k_lower;
                    %Define k_plus with special evolving value at the shock
                    %for ind_shock
                    k_plus(1:ind_shock - 2) = k_upper;
                    k_plus(ind_shock - 1) = k_upper * (dx / xi) * (p(ind_shock)- p_star) ...
                                       / (p(ind_shock) - p(ind_shock + 1));
                    k_plus(ind_shock:end) = k_lower;
                    %define tolerance parameter
                    %If close to x_i: area negligble use standard FV use
                    %k_lower other portion neglible p_i close to p_star

                    if (skip_leftcell)
                        k_plus(ind_shock - 1) = k_lower;
                    %If shock close to x_{i+1}, use k_upper other area
                    %negligble
                    end
                    if (skip_rightcell)
                        k_minus(ind_shock) = k_upper;
                    end
                end
                %define smaller cell area
                %Artificial diffusion in space to chock cell i + 1
                d2p_dx2 = (p(ind_shock + 2) - 2 * p(ind_shock + 1) + p(ind_shock)) / dx^2;
                d2p_dx2 = (pold(j+1) - 2 * pold(j) + pold(j-1)) / dx^2;
                k_t(n) = k_plus(cell-2);
                if(ind_shock == cell - 1) %cell is i+1
                    t_shock = [t_shock t(n)];
                end
                %k_t(n) = k_minus(cell-1);
                %Calcualte updated fluxes
                flux_in(n) = k_minus(cell-1) .* (pold(cell-1) - pold(cell)) / dx;
                flux_out(n) = k_plus(cell-1) .* (pold(cell) - pold(cell + 1)) / dx;
                %Compute F_i^+
                flux_in(n) = k_plus(cell-1) .* (pold(cell) - pold(cell + 1)) / dx;
                %Compute F_{i+1}^-
                flux_out(n) = k_minus(cell) .* (pold(cell) - pold(cell + 1)) / dx;
       
                p(j) = pold(j) + (1 / alpha) * (dt * F(j-1) + lambda .* (k_minus(j-1).*pold(j-1) ...
                        - (k_minus(j-1) + k_plus(j-1)) .* pold(j) + k_plus(j-1) .* pold(j+1))) + dt * d2p_dx2*diff_coeff;
                %investigate relation between flux and speed for exatc sol.
                fR = -(k_lower .* (p(ind_shock + 2) - p(ind_shock + 1)) / dx);
                fL = -(k_upper .* (p(ind_shock) - p(ind_shock - 1)) / dx);
                if (~isempty(ind_shock))
                    num_speed(n) = (fR - fL) / (p(ind_shock + 1) - p(ind_shock));
                end
            else
                %Add artificial diffusion in time
                if (n == 1)
                    pold(:,1) = p; %first time step
                end
                if (n == 2)
                    pold(:,2) = p;
                end
                if (n == 3)
                    pold(:,3) = p;
                end
                if (n >= 4) %standard switch in time
                    pold(:,1) = pold(:,2);
                    pold(:,2) = pold(:,3);
                    pold(:,3) = p;
                end
                if (n <= 3)
                    d2p_dt2 = 0;
                else
                    d2p_dt2 = (pold(i,3) - 2 * pold(i,2) + pold(i,1)) / dt^2;
                end
                if (artDiffSpace)
                    d2p_dx2 = (p(i+1) - 2 * p(i) + p(i-1)) / dx^2;
                end
                %smoothes solution removes oscillations in evolved IC case
                diffCoeff = 0.5 * 1e-1 * artDiffSpace;
                %Add artificial diffusion in space
                d2p_dx2 = (p(i+1) - 2 * p(i) + p(i-1)) / dx^2;
                if (weightedTimeAvg)
                    %Ignore at first timestep need to store initial
                    %condition first
                    if (weightedTimeAvgOld)
                        if (n > 2)
                            Dp_old = (k_minusold .* p_old(i-1) - (k_minusold + k_plusold) .* p_old(i)...
                                + k_plusold .* p_old(i+1)) / dx^2;
                            Dp_old2 = (k_minusold2 .* p_old2(i-1) - (k_minusold2 + k_plusold2) .* p_old2(i)...
                                + k_plusold2 .* p_old2(i+1)) / dx^2;
                        else
                            Dp_old = 0; Dp_old2 = 0;
                        end
                    else
                        k_minusold = k_minus;
                        k_plusold = k_plus;
                        if (n == 1)
                            c = 0;
                            Dp_old = 0;
                        else %begin at second timestep
                            c = 2 * 1e-5; %artficial diffusion coeff
                            Dp_old = (k_minusold .* p_old(i-1) - (k_minusold + k_plusold) .* p_old(i)...
                                    + k_plusold .* p_old(i+1)) / dx^2;
                        end
                    end
                    Dp = (k_minus .* p(i-1) - (k_minus + k_plus) .* p(i) + ...
                                 k_plus .* p(i+1)) / dx^2;
                    if (weightedTimeAvgOld)
                        if (n == 1)
                            p_old2 = p;
                            k_minusold2 = k_minus;
                            k_plusold2 = k_plus;
                            c = 0;
                        elseif (n == 2)
                            p_old = p;
                            k_minusold = k_minus;
                            k_plusold = k_plus;
                            c = 0;
                        else
                            p_old2 = p_old;
                            p_old = p;
                            k_minusold2 = k_minusold;
                            k_minusold = k_minus;
                            c = 1e-4;
                        end
                        p(i) = p(i) + (1 / alpha) * (dt * F + dt * Dp...
                                        + c * (Dp_old - Dp_old2));
                     else
                         p_old = p;
                         p(i) = p(i) + (1 / alpha) * (dt * F + ...
                                         (dt + c) * Dp - c * Dp_old);
                     end
                else
                    diffTimeCoeff = 1e-6 * artDiffTime;
                    p(i) = p(i) + (1 / alpha) * (dt * F + lambda * (k_minus .*p(i-1) ...
                       -(k_minus + k_plus) .* p(i) + k_plus .* p(i+1))) + ...
                        + dt * diffTimeCoeff * d2p_dt2 + dt * diffCoeff * d2p_dx2; %if false will be adding on 0;
                end
            end
        elseif(k_pow == -1) %No avg conserv form
            p(i) = compute_newP_noAvg(i, p, F, dt, alpha, lambda, k_lower, k_upper, p_star);
        else %No avg conserv case for k = p^n
            %Test artificila diffusion in time
            p(i) = p(i) + (1 / alpha) * (dt * F + lambda / (k_pow + 1) * (p(i+1).^(k_pow + 1) ...
                                        -2 * p(i).^(k_pow + 1) + p(i-1).^(k_pow +1)));      
        end
        if (modEqtn)
            p(i) = p(i) + dt  * (p_x_coeff .*dp_dx);
            p(i) = p(i) + dt  * (p_xx_coeff .*d2p_dx2 +p_x_coeff .*dp_dx);
        end
        %% Update boundaries
        p(1) = bdry{1}(t(n+1));
        p(end) = bdry{2}(t(n+1));
        %hold off
        %plot(x_cent,p,'-o')
        %        pause(1e-6)
    end
    %Consturct exact solution as well at final time
    p_t(end) = p(probe);
    if (k_pow == -1)
        [p_t_exact(end) p_exact] = exact_solution(t(end), x_cent, cell, alpha_star,...
                                              A_exact, B_exact, k_lower, lim);
        if (fluxEqual)
                [k_minus k_plot ratio] = computek_eff(p, i - 1, p_star, ...
                                                        k_lower, k_upper);
                k_plus = [k_minus(2:end); ...
                      computek_eff(p, length(p) - 1, p_star, k_lower, k_upper)];
            elseif(useOldv)
                %k_i-1/2 averaged on cell face
                [k_minus k_plot] = computek_eff_v(p, i - 1, p_star, k_lower, ...
                                                    k_upper, vold, dx);

                %k_i+1/2 averaged on right cell face
                k_plus = [k_minus(2:end); ...
                            computek_eff_v(p, length(p) - 1, p_star, k_lower, ...
                                                    k_upper, vold, dx)];            
            else
                [k_minus k_plot ratio] = computek_eff_exactv...
                                         (p, i - 1, p_star, k_lower, k_upper, ...
                                         vel_exact(end), dx);
                k_plus = [k_minus(2:end); ...
                      computek_eff_exactv(p, length(p) - 1, p_star, k_lower, ...
                                         k_upper, vel_exact(end), dx)];
        end
        [k_minus_exact k_plot_exact exact_ratio] = computek_eff(p_exact, i - 1, p_star, k_lower, k_upper);
        k_plus_exact = [k_minus_exact(2:end); ...
              computek_eff(p_exact, length(p) - 1, p_star, k_lower, k_upper)];
        %compute on last cell p is length N+1 so on(N,N+1)  
        k_eff(end) = k_plot;
        k_eff_exact(end) = k_plot_exact;
        if (~useOldv)
            dx_star_t(end) = ratio;
        end
        dx_star_exact(end) = exact_ratio;
        %Fix using k at next time using final pressurecheck tiemk_minus
        flux_in(end) = k_minus(cell - 1) .* (p(cell - 1) - p(cell)) / dx;
        flux_out(end) = k_plus(cell - 1) .* (p(cell) - p(cell + 1)) / dx;
        flux_in_exact(end) = k_minus_exact(cell - 1) .* (p_exact(cell-1) - p_exact(cell)) / dx;
        flux_out_exact(end) = k_plus_exact(cell - 1) .* (p_exact(cell) - p_exact(cell + 1)) / dx;

        figure(4);
        plot(t, flux_in, t, flux_out, t, flux_in - flux_out, 'Linewidth', 2);
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca, 'FontSize',16,'FontWeight','bold');
        xlabel('time t')
        ylabel('Fluxes into and out of cell')
        title(sprintf('Flux in and flux out versus time at grid cell with cell-center %.2f', ...
                        x_cent(cell)))
        legend('flux_{in}', 'flux_{out}', 'Location', 'Best')
        %saveas(gcf, '~/Desktop/NumFluxes.jpg')

%         figure(5)
%         plot(p,k, '-o', 'Linewidth', 2)
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         set(gca, 'FontSize',16,'FontWeight','bold');
%         title('Plot of k vs p for uniform p')
%         xlabel('pressure p');
%         ylabel('permeability k');
        
%         figure(6)
%         plot(t, dx_star_t, '-o', t, dx_star_exact, '-o', 'Linewidth', 2)
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         set(gca, 'FontSize',16,'FontWeight','bold');
%         title('Plot of dx*/dx vs t')
%         ylabel('dx*/dx'); %ratio 
%         xlabel('time t');
%         legend('num', 'exact')

%         figure(7)
%         plot(t, k_eff, '-o', t, k_eff_exact, '-o', 'Linewidth', 2)
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         set(gca, 'FontSize',16,'FontWeight','bold');
%         title('Plot of k_{eff} vs t')
%         ylabel('k_{eff}');
%         xlabel('time t');
%         legend('k_{eff}', 'exact k_{eff}')
        %saveas(gcf, '~/Desktop/k_eff.jpg')

%         figure(8);
%         plot(t, flux_in_exact, t, flux_out_exact, 'Linewidth', 2);
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         set(gca, 'FontSize',16,'FontWeight','bold');
%         xlabel('time t')
%         ylabel('exact Fluxes into and out of cell')
%         title(sprintf('exact Flux in and flux out versus time at grid cell with cell-center %.2f', ...
%                         x_cent(cell)))
%         legend('exact flux_{in}', 'exact flux_{out}', 'Location', 'Best')

        %Plot velocity using dx* as approximately computed using constraint
        %that dp*/dt = 0
        %dx_star = dx* / dx. To get dx* need to multiply by dx
        %find grid cell transitions
        vel = diff(dx_star_exact);
        grid_cells = find(vel < 0);
        vel_num = diff(dx_star_t);
        grid_cells_num = find(vel_num < 0);
        n_cells_num = length(grid_cells_num);
        n_cells = length(grid_cells);
        dx_star_t = dx_star_t * dx;
        dx_star_exact = dx_star_exact * dx;
        %Don't add dx at junction dx* = 1 for celli-1 is point i and dx*=0 in
        %next cell is dx* = 0 is point i at grid cell transition
        for i = 1:n_cells - 1
            dx_star_exact(grid_cells(i) + 1:grid_cells(i+1)) = ...
                    dx_star_exact(grid_cells(i) + 1:grid_cells(i+1)) + i*dx;
        end
        for i = 1:n_cells_num - 1
            dx_star_t(grid_cells_num(i) + 1:grid_cells_num(i+1)) = ...
                    dx_star_t(grid_cells_num(i) + 1:grid_cells_num(i+1)) + i*dx;
        end
        if (~isempty(grid_cells))
            dx_star_exact(grid_cells(end)+1:end) = dx_star_exact(grid_cells(end)+1:end) ...
                                            + n_cells * dx;
        end
        if (~isempty(grid_cells_num))
            dx_star_t(grid_cells_num(end)+1:end) = dx_star_t(grid_cells_num(end)+1:end) ...
                                            + n_cells_num * dx;
        end
        shock_loc = alpha_star*t.^(0.5);
        dx_star_exact = dx_star_exact + shock_loc(1);
        dx_star_t = dx_star_t + shock_loc(1);

        %Plot exact shock location alpha_star*sqrt(t) versus time
%         figure(9)
%         plot(t, shock_loc, '-o', t, dx_star_exact, '-o', t, dx_star_t, '-o','Linewidth', 2)
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         set(gca, 'FontSize',16,'FontWeight','bold');
%         title('Plot of shock location position vs time')
%         ylabel('shock location');
%         xlabel('time t');
%         legend('Exact','Exact Constraint', 'Num Constraint')
%         legend('Location','Southeast')

        %saveas(gcf, '~/Desktop/ExactPosAndVel/New/pos.jpg')

        %Plot velocities
%         figure(10)
%         plot(t(2:end), diff(shock_loc) / dt, '-o', t(2:end), diff(dx_star_exact) / dt, '-o', t(2:end), diff(dx_star_t) / dt,'-o', 'Linewidth', 2)
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         set(gca, 'FontSize',16,'FontWeight','bold');
%         title('Plot of velocity vs time')
%         ylabel('velocity');
%         xlabel('time t');
%         legend('Exact', 'Exact Constraint', 'Num Constraint')
        %saveas(gcf, '~/Desktop/ExactPosAndVel/New/vel.jpg')

        %Plot exact velocities without large numerical
%         figure(11)
%         plot(t(2:end), diff(shock_loc) / dt, '-o', t(2:end), diff(dx_star_exact) / dt,'-o', 'Linewidth', 2)
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         set(gca, 'FontSize',16,'FontWeight','bold');
%         title('Plot of velocity vs time')
%         ylabel('velocity');
%         xlabel('time t');
%         legend('Exact', 'Exact Constraint')
        %saveas(gcf, '~/Desktop/ExactPosAndVel/New/vel_exact.jpg') 

        %figure(12)
        %plot(t, num_speed, '-o');
        if (~fluxEqual)
            %xi_plot(end) = temp + speed * dt;
            %xi_plot = xi_plot + shock_loc(1); %integrate in time need to add back initial constant
            %xi_error(end) = xi - xi_exact;
            speed_plot(end) = speed;
            figure(100)
            plot(t-0.0479, shock_loc, t-0.0479, xi_plot, '--r')
            legend('Exact','Numerical', 'Location','Southeast')
            set(findall(gca, 'Type', 'Line'),'LineWidth',2);
            set(gca, 'FontSize',16,'FontWeight','bold');
            title('Shock Position')
            xlabel('t')
            ylabel('x^*')
            figure(101)
            %plot(t(1:end-1), diff(shock_loc) / dt, t(1:end-1), diff(xi_plot) / dt)
            vel_exact = alpha_star * 0.5./ sqrt(t);
            plot(t, vel_exact, t, speed_plot)
            title('Speed')
            legend('Exact', 'Num')
            figure(102)
            plot(t-0.0479, xi_plot - shock_loc)
            set(findall(gca, 'Type', 'Line'),'LineWidth',2);
            set(gca, 'FontSize',16,'FontWeight','bold');
            xlabel('t')
            ylabel('Shock Position Error')
            title('Shock Position Error as a function of time')
            figure(103)
            plot(t, abs(vel_exact - speed_plot))
            title('Error in the Speed')
            %save('xi_longtime02.mat','p','p_t','x_cent','t','xi_plot');
        end
    end
    save('ktplus_arith50left.mat','t','k_t', 'p_t');
end
