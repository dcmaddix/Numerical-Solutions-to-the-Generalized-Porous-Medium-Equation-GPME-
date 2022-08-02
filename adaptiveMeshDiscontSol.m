%Only for the discont case k_pow == -1
%Mimetic 2016 paper implementing Scheme I for various ki functions using
%basic time implementaions: Solving alpha*dp/dt -  d/dx(kdp/dx) = F
function [p p_t p_exact p_t_exact t alpha_star] = adaptiveMeshDiscontSol(bdry, p, p_shock, ...
                                                                  nt, dt, dx, dx_shock, x_coord, ...
                                                                  x_cent, x_shock, k_lower, t0,...
                                                                  IC_evolve, exact, int)
%% INPUTS:
 
%% OUTPUTS:
     
%% Initialization    
    fid = fopen('xi.txt', 'w');
    cell = find(x_cent < x_coord + 1e-6 & x_cent > x_coord - 1e-6, 1, 'first');
    if (isempty(cell)) %make sure xcoord exists
        fprintf('Ind warning\n')
        return
    end
    %Compute time vector
    %start at later time
    t = t0 : dt : t0 + nt;
    p_t = zeros(1,length(t)); %time vector where initial condition is given at t=0
    p_t_exact = zeros(1,length(t));
    %Initialize exact solution to 0 doesn't exist in p^n case
    %compute constant lambda
    lambda = dt / dx^2;
    %% Vector of interior indices
    i = 2 : length(x_cent) - 1;
    %Initialize to 0
    k_minus = zeros(length(i),1);
    k_plus = k_minus;
    %Time loop: already have it at first time so less iteration
    %special Jakolein case
    p_star = 0.5;
    k_upper = 1.0;
    %% exact-solution solve nonlinear for exact solution
    if (k_lower < 1e-6)
        lim = true; %use limiting solution where RHS solution is 0
    else 
        lim = false; %use solution for exact eps solving for
    end
    [alpha_star A_exact B_exact] = exact_sol_speed(lim, p_star, k_lower);
    for n = 1:length(t)-1
        %store at current time
        p_t(n) = p(cell);
        if (n > 1)
            old_indshock = ind_shock;
        end
        ind_shock = find(p(i - 1) >= p_star & p_star > p(i));
        [p_t_exact(n) ~] = exact_solution(t(n), x_cent, cell, alpha_star,...
                                          A_exact, B_exact, k_lower, lim);
       
        %% Update pressure step
        %solve (N-1)x(N-1) linear system for interior pressures
        shock_loc = alpha_star * t(n)^(0.5);
        x_shock = x_cent(ind_shock): dx_shock: x_cent(ind_shock + 1);
        x_shock = x_shock(2:end-1);
        %compare shock pos and interpolate when change cells
        if (n > 1 && ind_shock ~= old_indshock)
            p_shock = interp1(x_cent, p, x_shock, 'spline');
        end
        %Search p_shock before defining the k values and grid spacing
        %Try integral average in shock cells
        shock_left = false;
        shock_right = false;
        if (p(ind_shock) >= p_star && p_star > p_shock(1))
            shock_left = true;
        else
            ind_ref = find(p_shock(1:end-1) >= p_star & p_star > p_shock(2:end));
            if (isempty(ind_ref)) %between p_shock(end) and p(ind_shock + 1)
                shock_right = true;
            end
        end
        pold = p;
        p_shockold = p_shock;
        ind_left = 2:ind_shock - 1;
        %% Standard stencil away from shock
        p(ind_left) = lapl_update(pold, dt, dx, ind_left, k_upper);
        %% Standard stencil away from shock
        ind_right = ind_shock + 2: length(p) - 1;
        p(ind_right) = lapl_update(pold, dt, dx, ind_right, k_lower);
        
        %% ind_shock
        %if shock is present
        if (n == 1) %change for other IC
            %use exact for initial timestep
            xiold = shock_loc - x_cent(ind_shock);
        end
        %% shock_left case simpler betwen p(ind_shock) and p_shock(1)
        if (shock_left)
            xi = computeXi(xiold, dt, dx, dx_shock, dx_shock, ...
                           pold(ind_shock - 1), pold(ind_shock), p_shockold(1), ...
                           p_shockold(2), k_upper, k_lower, p_star);
            if (exact)
                xi = shock_loc - x_cent(ind_shock);
            end
            if (int)
                F_min = calcFluxes(pold(ind_shock), p_shockold(1), dx_shock, ...
                               true, p_star, k_lower, k_upper);
            else
                F_min = calcFluxesXi(pold(ind_shock), p_shockold(1), dx_shock, ...
                                true, p_star, k_lower, k_upper, xi, false);
            end
                            %fix all xi ones with right function
            F_plus = calcFluxes(p_shockold(1), p_shockold(2), dx_shock, ...
                                false, p_star, k_lower, k_upper);
            p_shock(1) = update(p_shockold(1), dt, dx_shock, F_min, F_plus); 
            p_shock(2:end-1) =  lapl_update(p_shockold, dt, dx_shock, ...
                                             2:length(p_shock) - 1, k_lower);
            p_shock(end) = p_shockold(end) + dt / dx_shock^2 * k_lower ...
                                * (p_shockold(end-1) - 2 * p_shockold(end) ...
                                + pold(ind_shock + 1));
        %% shock_right case between p_shock(end) and p(ind_shock + 1)
        elseif (shock_right)
            %update on shock grid using same xi
            xi = computeXi(xiold, dt, dx_shock, dx_shock, dx, ...
                           p_shockold(end - 1), p_shockold(end), ...
                           pold(ind_shock + 1), pold(ind_shock + 2), ...
                           k_upper, k_lower, p_star);
            if (exact)
                xi = shock_loc - x_shock(end);
            end
            %store for later use on the plus side of shock
             p_shock(1) = p_shockold(1) + dt / dx_shock^2 * k_upper ...
                                * (pold(ind_shock) - 2 * p_shockold(1) ...
                                + p_shockold(2));
             p_shock(2:end-1) =  lapl_update(p_shockold, dt, dx_shock, ...
                                              2:length(p_shock) - 1, k_upper);
             F_min = calcFluxes(p_shockold(end-1), p_shockold(end), dx_shock, ...
                                 false, p_star, k_lower, k_upper);
             F_plus = calcFluxes(p_shockold(end), pold(ind_shock + 1), dx_shock, ...
                                 true, p_star, k_lower, k_upper);
             p_shock(end) = update(p_shockold(end), dt, dx_shock, F_min, F_plus);
        end
        %% ind_shock
        %xi = shock_loc - x_cent(ind_shock);
        F_min = calcFluxes(pold(ind_shock - 1), pold(ind_shock), dx, ...
                          false, p_star, k_lower, k_upper);
        F_plus = calcFluxes(pold(ind_shock), p_shockold(1), dx_shock, ...
                              shock_left, p_star, k_lower, k_upper);
        h = (dx + dx_shock) / 2;
        p(ind_shock) = update(pold(ind_shock), dt, h, F_min, F_plus);
        %% ind_shock + 1
        if (int)
            F_min = calcFluxes(p_shockold(end), pold(ind_shock + 1), dx_shock, ...
                               shock_right, p_star, k_lower, k_upper);
        else
            F_min = calcFluxesXi(p_shockold(end), pold(ind_shock + 1), dx_shock, ...
                                shock_right, p_star, k_lower, k_upper, xi, false);
        end
        F_plus = calcFluxes(pold(ind_shock+1), pold(ind_shock+2), dx, ...
                                false, p_star, k_lower, k_upper);
        p(ind_shock + 1) = update(pold(ind_shock + 1), dt, h, F_min, F_plus);
        %% Fine grid
        if (~shock_left && ~shock_right)
            %% p_shock(1)
            F_min =  calcFluxes(pold(ind_shock), p_shockold(1), dx_shock, ...
                                false, p_star, k_lower, k_upper);
            %ind_ref determines if shock is between 1 and 2
            if (ind_ref == 1)
                xi = computeXi(xiold, dt, dx_shock, dx_shock, dx_shock, ...
                               pold(ind_shock), p_shockold(ind_ref), ...
                               p_shockold(ind_ref + 1),p_shockold(ind_ref + 2), ...
                               k_upper, k_lower, p_star);
            end
            if (exact)
                xi = shock_loc - x_shock(1);
            end
            F_plus = calcFluxes(p_shockold(1), p_shockold(2), dx_shock, ...
                              ind_ref == 1, p_star, k_lower, k_upper);
            p_shock(1) = update(p_shockold(1), dt, dx_shock, F_min, F_plus);
            %% 2:ind_ref - 1: all k_upper
            %if ind_ref = 1, 2 will be overwritten in ind_ref + 1 ok
            ind_left = 2:ind_ref - 1;
            p_shock(ind_left) = lapl_update(p_shockold, dt, dx_shock, ind_left, ...
                                                 k_upper);
            %Already took care of case when ind_ref is 1 so only defined on
            %interior
            %% ind_ref + 2: end-1
            ind_right = ind_ref + 2:length(p_shock) - 1;
            p_shock(ind_right) = lapl_update(p_shockold, dt, dx_shock, ...
                                            ind_right, k_lower);
            %if ind_ref will be corrected later
            F_min = calcFluxes(p_shockold(end - 1), p_shockold(end), dx_shock, ...
                                false, p_star, k_lower, k_upper);
            F_plus = calcFluxes(p_shockold(end), pold(ind_shock + 1), dx_shock, ...
                                false, p_star, k_lower, k_upper);
            p_shock(end) = update(p_shockold(end), dt, dx_shock, F_min, F_plus);
            %% ind_ref 
            if (ind_ref > 1)
                F_min =  calcFluxes(p_shockold(ind_ref - 1), p_shockold(ind_ref), dx_shock, ...
                                false, p_star, k_lower, k_upper);
                if (ind_ref <= length(p_shock) - 2)
                    xi = computeXi(xiold, dt, dx_shock, dx_shock, dx_shock, ...
                                   p_shockold(ind_ref - 1), p_shockold(ind_ref), ...
                                   p_shockold(ind_ref + 1), p_shockold(ind_ref + 2), ...
                                   k_upper, k_lower, p_star);
                else
                    xi = computeXi(xiold, dt, dx_shock, dx_shock, dx_shock, ...
                                   p_shockold(ind_ref - 1), p_shockold(ind_ref), ...
                                   p_shockold(ind_ref + 1), pold(ind_shock + 1), ...
                                   k_upper, k_lower, p_star);
                end
                if (exact)
                    xi = shock_loc - x_shock(ind_ref);
                end
                F_plus = calcFluxes(p_shockold(ind_ref), p_shockold(ind_ref + 1), ...
                                      dx_shock, true, p_star, k_lower, k_upper);
                %always defined since ind_ref is in 1 to length(p) - 1
                p_shock(ind_ref) = update(p_shockold(ind_ref), dt, dx_shock, F_min, F_plus);
            end
            %% ind_ref + 1- out of bounds for ind_ref = length(p) - 1
            if (int)
                F_min = calcFluxes(p_shockold(ind_ref), p_shockold(ind_ref + 1), dx_shock, ...
                                true, p_star, k_lower, k_upper);
            else
                F_min = calcFluxesXi(p_shockold(ind_ref), p_shockold(ind_ref + 1), dx_shock, ...
                                true, p_star, k_lower, k_upper, xi, false);
            end
            if(ind_ref <= length(p_shock) - 2)
                F_plus = calcFluxes(p_shockold(ind_ref + 1), p_shockold(ind_ref + 2), dx_shock, ...
                                false, p_star, k_lower, k_upper);
            elseif (ind_ref == length(p_shock) - 1) %F_min will stay the same
                F_plus = calcFluxes(p_shockold(ind_ref + 1), pold(ind_shock + 1), dx_shock, ...
                                false, p_star, k_lower, k_upper);%length(p_shock) is rightShock case
            end
            %if ind_ref = length(p) - 1 get last point covered if ind_ref =
            %1, 3 gets correctly overwritten here
            p_shock(ind_ref + 1) = update(p_shockold(ind_ref + 1), dt, dx_shock, F_min, F_plus); 
        end
        xiold = xi;
        %% Update boundaries
        p(1) = bdry{1}(t(n+1));
        p(end) = bdry{2}(t(n+1));
        %fprintf(fid, '%f\n',xi);
    end
    %Construct exact solution as well at final time
    p_t(end) = p(cell);
    [p_t_exact(end) p_exact] = exact_solution(t(end), x_cent, cell, alpha_star,...
                                          A_exact, B_exact, k_lower, lim);     
end

       %% Stencil for left of shock: grid spacing all dx fast code
        %Define k_minus, with special evolving value at the
        %shock for ind_shock + 1
%         k_minus(1:ind_shock - 1) = k_upper;
%         k_minus(ind_shock) = k_lower * (dx / (dx - xi)) * (p_star - p(ind_shock +1)) ...
%                                 / (p(ind_shock) - p(ind_shock + 1));
%         k_minus(ind_shock + 1:end) = k_lower;
%         %Define k_plus with special evolving value at the shock
%         %for ind_shock
%         k_plus(1:ind_shock - 2) = k_upper;
%         k_plus(ind_shock - 1) = k_upper * (dx / xi) * (p(ind_shock)- p_star) ...
%                            / (p(ind_shock) - p(ind_shock + 1));
%         k_plus(ind_shock:end) = k_lower;
%         %define tolerance parameter
%         %If close to x_i: area negligble use standard FV use
%         %k_lower other portion neglible p_i close to p_star
% 
%         if (xi < tol)
%             k_plus(ind_shock - 1) = k_lower;
%         %If shock close to x_{i+1}, use k_upper other area
%         %negligble
%         end
%         if (dx - xi < tol)
%             k_minus(ind_shock) = k_upper;
%         end
%         p(i) = p(i) + lambda * (k_minus(i-1).*p(i-1) ...
%             - (k_minus(i-1) + k_plus(i-1)) .* p(i) + k_plus(i-1) .* p(i+1));
