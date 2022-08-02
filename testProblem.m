function [p p_t x_cent t p_t_exact p_exact pos_5 alpha_star xi_plot] = testProblem(N,alpha, BC_type, timeScheme, ...
                     k_pow, nt_inner, timeMethod, avg, x_coord, errorType, ...
                     use_exact, upwind, weight, theta, nt, k_deriv, eps_pxx, ...
                     eps_px, eps_pxxx, eps_pxx2, conserv, eigen, eps, ramp, ...
                     k_lower, FDgrid, modEqtn, globalAtan, simple, noAvg, ...
                     Jameson, Jsolver, limiter, upwindFlux, fluxEqual, IC_evolve, ...
                     t0, useOldv, artDiffTime, artDiffSpace, newUpdateStep, ...
                     weightedTimeAvg, weightedTimeAvgOld, adaptiveMesh, exact, ...
                     int, N_IC, slowDiff, mimeticLock, lb, ub, levelset, linear)
    %% HEAT EQUATION TEST CASE: CONSTANT COEFF
    %lb = 0.0; %for bdry heat eqtn test case
    %ub = 1.0; %for bdry heat eqtn test case
    alpha_star = 0.0;
    %uniform grid spacing
    %dx = (ub - lb) / N; %pressure is size N+1
    dx = 1 / N;
    xi_plot = [];
    if (conserv)
        dt = alpha*dx^2 / 2;
        if (k_pow == 1)
            dt = dt / 2;
        elseif(k_pow == 2)
            dt = dt / 4;
        elseif(k_pow == 3)
            %for cancelation
            dt = dt / 8
            %for simple version
            %dt = dx^2 / 4; needs to be changed with cancellation
        elseif (k_pow == -1)
            %discont case
            dt = 1.8e-4;
            dt = dx^2 / 16; %need this to cancel terms in arctan for eps = 0.09
            dt = dx^2 / 4;
            %dt = 1.8e-4 %Jak ts for 50 points
        end
    else
        dt = alpha * dx^2 / 25; %need to divide by 4 or 8 for explicit
        dt = alpha * dx^2 / 4;
    end
    if (slowDiff)
        dt = dt / 2;
    end
    %dt = dt /2;
    %removes oscillations at shock
    if (k_pow == -1)
        dt = dt / 8; % for discont case 8 for N = 100, 2 for N = 50, N = 25
    end
    %dt = dt / 4; %removed for p^3 case for N = 25 exp case don't need it otherwise
    if (adaptiveMesh)
        n_shock = 10;
        dt = dt / n_shock;
    end
    %dt = 8*dt;
    %dt = dx / 2; %.1dx has oscillations use 0.5 for BackwardEuler
    %define spatial vector at cell faces to plot against velocity
    if (FDgrid) %true
        %FD grid pressure defined at N+1 points [0,dx,..,1-dx,1]
        %FV grid where first cell center is bdry point 0 with infinitesmal
        %cell volume around it
        x = linspace(lb + dx/2, ub - dx/2, N); %flux kp_x = u defiend on faces
        x_cent = [lb 0.5 * (x(1:end-1) + x(2:end)) ub]; %pressure cell center
    else
        %This case would require more complicated bdry condition treatment
        %by enforcing the condition on the left and right flux or ghost
        %cells. If this case is turned on we get Jakolein's results where
        %the curve is above the exact solution and mismatched
        x_cent = linspace(lb + dx/2, ub - dx/2, N); %pressure cell center
        x = [lb 0.5 * (x_cent(1:end-1) + x_cent(2:end)) ub]; %flux kp_x = u defiend on faces
    end
    x_cent = lb:dx:ub;
    %Define harmonic and arithmetic averaging functions
    k_H = @(k1,k2) (2*k1.*k2) ./ (k1+k2);
    k_A = @(k1,k2) (k1+k2) / 2;
    k_FV = @(k1,k2) (k_H(k1,k2).*k_A(k1,k2)).^(0.5);
    %Determine k_i function
    if (strcmp(avg, 'mfd harm'))
        k_i = @(k1,k2) k_H(k1,k2);
    elseif(strcmp(avg, 'arith'))
        k_i = @(k1,k2) k_A(k1,k2); %arithmetic averaging causes smoothing
    elseif(strcmp(avg, 'fv harm'))
        k_i = k_FV;%input to get simple FV scheme
    elseif(strcmp(avg, 'upwind') || strcmp(avg, 'int') || strcmp(avg, 'k_eff'))
        k_i = 0; %defined inside function
    else
        fprintf('ERROR: Undefined averaging type!\n')
        return
    end
    %Jak case
    p_t_exact = zeros(length(0:dt:nt),1);
    p_exact = zeros(length(x_cent),1);
    F = zeros(length(x_cent)-2,1); %constant zero forcing term
    %call scheme I with homogenous Dirichlet BC at both endpoints to make
    %consistent
    %% Initial Condition for k = p^n
    if (k_pow ~= -1)
        %works in monotonically increasing case as well 0.1 and 2.0
        if (mimeticLock)
            bdry{1} = @(t) (3*t)^(1/3);
            bdry{2} = @(t) 1e-3;
        else
            bdry{1} = @(t) 2.0;
            bdry{2} = @(t) 0.1;
        end
        init_cond = @(x) bdry{1}(0)*exp(log(bdry{2}(0)/bdry{1}(0)) * x  / 0.1);
    else %flag to use Jakolein's discontinuous perm case
        bdry{1} = @(t) 1.0;
        bdry{2} = @(t) 0.0;
        init_cond = @(x) -((bdry{1}(0) - bdry{2}(0)) / 0.1) * x + bdry{1}(0);
    end
    
    p = init_cond(x_cent)'; %3 points between 0 and 0.1: coarsest grid is N = 40
    p(x_cent > 0.1) = bdry{2}(0);
    if(k_pow == -1)
        p = zeros(length(x_cent),1);
    end
    if (mimeticLock)
        p(1) =  1e-3;
        p(2:end) = 1e-3;
    end
    %tried linear smooth IC and analysis still holds!
    %p = 2 - 1.9* x_cent';
    %Added to work with case where shock at index 1 and goes out of bounds!
    %p(2) = p(1);
    ghostNodes{1} = @(t) bdry{1}(t); %exponential left
    ghostNodes{2} = @(t) bdry{2}(t); %constant right
    %% Test later initial condition and plot it
    if (IC_evolve && N  <= 200)
        if(N_IC == 50)
            p = load('IC_0.0479/discont_50.mat');
        elseif(N_IC == 100)
            p = load('IC_0.0479/discont_100.mat');
        elseif (N_IC == 200)
            p = load('IC_0.0479/discont_200.mat');
        elseif(N_IC == 25 || N_IC == 12.5)
            p = load('IC_0.0479/discont_25.mat');
        end
        %p = load(sprintf('IC_0.0479/discont_%i.mat',N));
     p = p.p_exact;
     if (N_IC == 12.5)
        p = p(1:2:end);
     end
    end
    %linear case foam initial linear condition
    if (linear)
        p(x_cent <= 0.5) = 1 - 2*x_cent(x_cent <= 0.5)';
    end
    diffN = length(x_cent) - length(p);
    p = [p; zeros(diffN,1)];
    figure(1)
    hold on;
    %plot([lb - dx x_cent ub + dx],[ghostNodes{1}(0); p; ghostNodes{2}(0)],'-ob', 'Linewidth', 2)
    plot(x_cent,p,'-b', 'Linewidth', 2)
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca, 'FontSize',16,'FontWeight','bold');
    title('Initial Condition')
    xlabel('x')
    ylabel('p')
    hold on
    p_star = 0.5;
    k_upper = 1.0;
    %dt = 4*dt;
    if (adaptiveMesh)
        %%New section
        i = 2:length(p) - 1;
        ind_shock = find(p(i - 1) >= p_star & p_star > p(i));
        dx_inner = dx / n_shock;
        x_shock = x_cent(ind_shock): dx_inner :x_cent(ind_shock + 1);
        x_shock = x_shock(2:end-1);
        %%How to generate exact solution as IC
        lim = false; %old IC generated for limit = true 
        k_IC = 0.01;
        [alpha_star A_exact B_exact] = exact_sol_speed(lim, p_star, k_IC);
        [~, p_shock] = exact_solution(t0, x_shock, [], alpha_star,...
                                                  A_exact, B_exact, k_IC, lim);
        full_grid = [x_cent(1:ind_shock) x_shock x_cent(ind_shock+1:end)];
        full_p = [p(1:ind_shock); p_shock; p(ind_shock + 1:end)];
        figure(100)
        hold on
        plot(full_grid, full_p, '-or') 
    end
%% Call Solvers
    if (adaptiveMesh)
        [p p_t p_exact p_t_exact t alpha_star] = adaptiveMeshDiscontSol(bdry, p, p_shock,  ...
                                                                  nt, dt, dx, dx_inner, x_coord, ...
                                                                  x_cent, x_shock, k_lower, t0,...
                                                                  IC_evolve, exact, int);
        pos_5 = zeros(1, length(p_t) - 1);
    elseif (strcmp(timeMethod, 'regular'))
        if (conserv)
            %Simple implementation without cancellation combine together
            %with separate k_p deriv routine
            if (Jameson)
                %define integral flux for discont case break up at p*
                if (k_pow ~= -1)
                    intK = @(p) 1 / (k_pow + 1) * p.^(k_pow+1);
                else
                        
                    intK = @(p) double(p >= p_star) .* (k_lower * p_star + ...
                                k_upper * (p - p_star)) + ...
                                k_lower * double(p < p_star) .* p;
                end
                %Need to add integral flux ofr discon case above only holds
                %for k = p^n
                if (strcmp(Jsolver, 'JST'))
                    [p p_t t] = JST(dx,nt, dt, x_cent, p, x_coord, ...
                                    upwindFlux, intK, lb, ub, t0);
                elseif (strcmp(Jsolver, 'SLIP'))
                    [p p_t t] = SLIP(dx,nt,dt,x_cent,p, x_coord, limiter, ...
                                    upwindFlux, intK, lb, ub, t0);
                elseif (strcmp(Jsolver, 'USLIP'))
                    [p p_t t] = USLIP(dx,nt,dt,x_cent,p, x_coord, limiter, ...
                                      upwindFlux, intK, lb, ub, t0);
                end
                pos_5 = zeros(1, length(p_t) - 1);
            elseif (simple)
                [p p_t p_exact p_t_exact t alpha_star xi_plot] = mimeticsimpleFE(k_i, k_pow, bdry, p, alpha, nt, dt,...
                                     dx, x_coord, x_cent, avg, eps, ramp, F, k_lower, modEqtn, globalAtan, noAvg,...
                                     fluxEqual, t0, useOldv, artDiffTime, artDiffSpace, newUpdateStep, ...
                                     weightedTimeAvg, weightedTimeAvgOld, IC_evolve, exact, levelset,linear);
                pos_5 = zeros(1, length(p_t_exact) - 1);
            else %This works for cancellation in p^n case
                [p, p_t t pos_5] = mimeticTransportCoeff(k_i, N, k_pow, bdry, ghostNodes, F, p, ...
                                               alpha,BC_type,timeScheme,nt, dt, dx,...
                                               use_exact, errorType, x_coord, theta, x_cent, ...
                                               avg, eps_pxx, eps_px, eps_pxxx, eigen, eps, ramp, slowDiff);
            end
        else
             [p p_t t] = nonConservForm(bdry, ghostNodes, nt, p, dt, dx, ... 
                                        x_coord, x_cent, k_deriv,eps_pxx,...
                                        eps_px, eps_pxx2, eps_pxxx, k_pow, ramp, eps);
        end
    elseif (strcmp(timeMethod, 'nonlinear'))
            k_p = @(p) p.^(k_pow)'; %as input to nonlinear solve
            p = mimeticTransportCoeffNonlinear(k_i, N,k_p, bdry_left, bdry_right, ...
                                               F,alpha,BC_type,timeScheme,nt,p, dt);
    elseif (strcmp(timeMethod, 'expandKeqtn'))
            [p, error]  = mimeticTransportCoeffPredCorr(k_i, N,k_pow, bdry_left, bdry_right, ...
                                                        F,alpha,BC_type,timeScheme,nt,p,...
                                                        dt,dx,nt_inner, x_coord, upwind, errorType, ...
                                                        use_exact);
    elseif (strcmp(timeMethod, 'kFVeqtn'))
            [p, error] = mimeticTransportCoeffPredCorrFVk(k_i, k_FV, N,k_pow, bdry_left, bdry_right, ...
                                                          F,alpha,BC_type,timeScheme,nt,p, dt, ...
                                                            x_coord, errorType, upwind, weight);
                                                        elsefif (strcmp(timeMethod, 'p_predcorr'))
            [p error] = mimeticTransportCoeffPredCorrPress(N,k_pow, bdry_left, bdry_right, ...
                                         F,alpha,BC_type,timeScheme,nt,p,dt);

    end
end
%% Heat Equation Test Problem
%     p = (sin(pi*x_cent) - sin(3*pi*x_cent))';
%     k = @(p)ones(N+1,1)';
%     bdry_left = @(t) exp(-pi^2*t)+exp(-9*pi^2*t); lb = 0.5 ub = 1.5
%     otherwise (0,0) for lb = 0 ub = 1
%     bdry_right = @(t)-exp(-pi^2*t)-exp(-9*pi^2*t);
%     ghostNode_left = @(t) exp(-pi^2*t).*sin(pi*(lb-h)) - exp(-9*pi^2*t).*sin(3*pi*(lb-h));
%     ghostNode_right = @(t) exp(-pi^2*t).*sin(pi*(ub+h)) - exp(-9*pi^2*t).*sin(3*pi*(ub+h));
%     lambda = 1/6;
%     nt = 1;
%     lambda = delta_t/h^2
%     delta_t = h^2/10;
%     pexact=@(t)exp(-pi^2*t).*sin(pi*x_cent) -
%     exp(-9*pi^2*t).*sin(3*pi*x_cent);

