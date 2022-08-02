%close all
%clear all
for iter = 1:1
    N = [25 50 100];
%Create parabolic plots
addpath('HelpfulFunctions')
addpath('DerivativeDiscretizations')
addpath('Jameson')
tic;
alpha = 1;
BC_type = 'Dir';
timeScheme = 'ForwardEuler';
%timeScheme = 'BackwardEuler';
%timeScheme = 'RK2_TVD';
%k_deriv = 'Central_2nd';
k_deriv = 'Central_4th';
%k_deriv = 'D_plus';
%k_deriv = 'D_minus';
k_deriv = 'explicit';

conserv = true;
eigen = false;
N =  50; %dx = 0.08
theta = 0.5;
k_pow = -1; %for constan coeff heat use k_pow 0
nt_inner = 1;
timeMethod = 'regular';
ramp = 'arctan';
%ramp = 'sine';
%ramp = 'linear';
globalAtan = false;
if (strcmp(ramp, 'arctan'))
    globalAtan = true;
end
avg = 'fv harm';
%avg = 'upwind';
avg = 'arith';
MHM = false;
avg = 'int';
%avg = 'k_eff';
errorType = 'l2norm';
use_exact = false; %exatc mimetic k!
upwind = true;
weight = 1.0;
lb = 0.0; %for bdry heat eqtn test case
ub = 1.0; %8 for long time %for bdry heat eqtn test case
%plots = zeros(length(N),1);   
colorstring = 'mkrgrmxckbg';
colorstring2 = 'bmkrbcmbgk';
count = iter;
if (k_pow ~= -1)
    nt = 8e-2;
else
    nt = 0.05;
    %nt = 0.056; time when p_{i+1} vlaue in gap
end
%nt = time(4);
FDgrid = true;
modEqtn = false;
%flux equal assumption using simple k_eff with dx* equiv to integral average
%and so velocity
%use exact xi
slowDiff = false;
int = true;
exact = false;
adaptiveMesh = false;
levelset = true;
linear = false;
fluxEqual = true;
IC_evolve = true; %discont case
useOldv = false; %locks useexactV works
artDiffTime = false;
artDiffSpace = false;
newUpdateStep = true;
weightedTimeAvg = false;
weightedTimeAvgOld = false;
mimeticLock = false;
if (IC_evolve)
    t0 = 0.0479;
else
    t0 = 2.5e-5;
end
if mimeticLock 
    nt = 0.7; %for mimetic locking problem
end
%Simple = false for truncation error cancellation for p^n MHM
simple = true; %True for SAM
%MHM param
eps_px = 0;
eps_pxx = 0;
if (strcmp(avg, 'fv harm') && MHM)
   eps_px = 1;
   eps_pxx = 1;
end
noAvg = false;
Jameson = true;
upwindFlux = true;
JSolver = 'JST';
JSolver = 'SLIP';
JSolver = 'USLIP';
limiter = 'minmod';
limiter = 'superbee';
limiter = 'vanleer';
%problems when not enough pressure points in the gap
%oscillations caused by size of gap: ep sis width of gap
eps = 0.0; %works for arctan case 0.09 arctan to cancel but too smooth
k_lower = 0.0; %k_eff won't work for 0, harm locks when 0
%k_lower = 0;
if (N(iter) >= 100)
    x_coord = 0.13;
else %not present on N = 50 grid
    x_coord = 0.12;
end
if (k_pow ~= -1)
    x_coord = 0.1;
end
%for k = p^n x = 0.12 or 0.13 0.33 for discont nonlinear, 0.1 for 40 p^3
%harm
if (IC_evolve)
    x_coord = 0.3; %.33 old
    x_coord = 0.32; %25 gird .32
else
    x_coord = 0.1; %.13 old
    x_coord = 0.12; %Jak
end
if (k_pow ~= -1)
    x_coord = 0.12;
end
dx = 1 / N(iter); %pressure is size N+1
x = 0:dx:1;
eps_pxxx = 0;
eps_pxx2 = 0;
%p2 = zeros(2,1);
%p3 = zeros(2,1);
%count = 1;
num_peaks = [];
%For truncation error cancellation upper bound is 1
for i = 0:0
    for j = 0:0
        for k = 0:0
        if (i == 1)
            if (strcmp(k_deriv, 'D_minus'))
                if (j == 0 || k == 0)
                    continue;
                end
            elseif(k_pow ~= -1)
                continue;
            end
        end
        if (~conserv && k == 1 && ~strcmp(k_deriv, 'D_minus') && ~strcmp(k_deriv, 'D_plus'))
            continue;
        end
%         eps_pxx2 = i; 
%         eps_pxx = j;
%         eps_px = k;
%         eps_pxxx = i;
        if (conserv && k_pow == 2 && k == 1)
            continue;
        end
        if (k_pow == -1 && (i == 1))
            continue;
        end
        %[p p_t x_cent t p_t_exact p_exact pos5] 
        if (simple)
            [p p_t x_cent t p_t_exact p_exact pos5 alpha_star, xi_plot] = testProblem(N(iter),alpha, BC_type, timeScheme, ...
                      k_pow, nt_inner, timeMethod, avg, x_coord, ...
                      errorType, use_exact, upwind, weight, theta, nt, ...
                      k_deriv, eps_pxx, eps_px, eps_pxxx, eps_pxx2, conserv, eigen,...
                      eps, ramp, k_lower, FDgrid, modEqtn, globalAtan, simple, noAvg, ...
                      Jameson, JSolver, limiter, upwindFlux, fluxEqual, IC_evolve, t0, useOldv, ...
                      artDiffTime, artDiffSpace, newUpdateStep, weightedTimeAvg, ...
                      weightedTimeAvgOld, adaptiveMesh, exact, int, N, slowDiff, ...
                      mimeticLock, lb, ub, levelset, linear);
        else
            [p p_t x_cent t p_t_exact p_exact pos5] = testProblem(N(iter),alpha, BC_type, timeScheme, ...
                      k_pow, nt_inner, timeMethod, avg, x_coord, ...
                      errorType, use_exact, upwind, weight, theta, nt, ...
                      k_deriv, eps_pxx, eps_px, eps_pxxx, eps_pxx2, conserv, eigen,...
                      eps, ramp, k_lower, FDgrid, modEqtn, globalAtan, simple, noAvg, ...
                      Jameson, JSolver, limiter, upwindFlux, fluxEqual, IC_evolve, t0, useOldv, ...
                      artDiffTime, artDiffSpace, newUpdateStep, weightedTimeAvg, ...
                      weightedTimeAvgOld, adaptiveMesh, exact, int, N, slowDiff, ...
                      mimeticLock, lb, ub, levelset, linear);
        end
        %save('ref1600_time','p_t', 't')
        %postprocess oscillations: Find locations of minimum by passing in
        %-p_t
        [oscill loc] = findpeaks(-p_t);
        %oscill;
        %output for p^n case
        period = abs(pos5(loc(1:end-1)) - pos5(loc(2:end)))
        t(loc)
        if (isempty(oscill))
            fprintf('No oscillations found!\n')
        else
            fprintf('Number of Oscillation Peaks found: %i\n', length(oscill))
        end
        %max(abs(p-p_exact))
        %psi = alpha_star * sqrt(nt)
        %check convergenece to left of shock in smooth region
%         figure(10)
%         hold on
%         plot(x(x > psi), abs(p(x > psi) - p_exact(x > psi)), '-r')
        %max(abs(p(x < psi) - p_exact(x < psi)))
        sqrt(dx) * norm(p - p_exact);
        f2 = figure(2); 
        hold on
        %p2(count) = plot(x_cent, p, strcat('-o',colorstring(count)));
        plot(x_cent, p_exact, '-o', x_cent, p, '-o')
        %plot(x_cent, p, '-o')
        %save('var_int.mat', 'p')
        %legend('Exact',sprintf('%s',avg))
        %legend('Exact','Numerical')
        %legend('Exact','Int. Avg.', 'Shock Avg.', 'Location', 'NorthWest')
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca, 'FontSize',16,'FontWeight','bold');
        xlabel('position x')
        ylabel('pressure p')
        if(~conserv)
            title(sprintf('Plot of numerical pressure with %s averaging for N = %i for n = %i', k_deriv, N, k_pow))
        else
            title(sprintf('Plot of numerical pressure with %s averaging for N = %i for n = %i', avg, N, k_pow))
        end
        if (Jameson)
            title(sprintf('Pressure vs. position for %s with %s limiter', JSolver, limiter))
            if(strcmp(JSolver, 'JST'))
                 title(sprintf('Pressure vs. pos for %s', JSolver))
            end
        end
        if (noAvg)
            title('Pressure vs pos for Central 2nd Deriv (Explicit Integral Average)')
        end
        title('Pressure vs. position ')
        xlabel('x')
        ylabel('p')
        hold on;
        %saveas(gcf, '~/Desktop/pos.jpg')
        %Save file
        %save(sprintf('data/foam/intamr_nx%i.mat', N), ...
        %    'x_cent', 'p', 'p_t', 't')
        f3 = figure(3);
        hold on;
        %p3(count) = plot(t,p_t, colorstring(count));
        plot(t-0.0479, p_t_exact, t-0.0479, p_t, '-r');
        %plot(t, p_t);
        %legend('Exact', sprintf('%s',avg), 'Location', 'NorthWest')
        %legend('Exact','Numerical', 'Location', 'NorthWest')
        %legend('Exact','Int. Avg.', 'Shock Avg.', 'Location', 'NorthWest')
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca, 'FontSize',16,'FontWeight','bold');
        ylabel(sprintf('pressure at position x = %.2f', x_coord))
        if (~conserv)
            title(sprintf('Plot of pressure versus time at x = %.2f for N = %i with %s averaging for n = %i', x_coord, N, k_deriv, k_pow))
        else
            title(sprintf('Plot of pressure versus time at x = %.2f for N = %i with %s averaging for n = %i', x_coord, N, avg, k_pow))
        end
        xlabel('time t')
        if (Jameson)
            title(sprintf('Pressure vs. time for %s with %s limiter', JSolver, limiter))
            if(strcmp(JSolver, 'JST'))
                 title(sprintf('Pressure vs. time for %s', JSolver))
            end
        end
        if (noAvg)
            title('Pressure vs time for Central 2nd Deriv (Explicit Integral Average)')
        end
        xlabel('t')
        title('Pressure vs. time ')
        %ylabel(sprintf('p(t, x = %.2f) '), x_coord)
        hold on;
        %saveas(gcf, '~/Desktop/time.jpg')
        count = count + 1;
        end
   end
end
if (k_pow == -1)
    %legend(p2, sprintf('%s',avg), 'Exact')
    %saveas(f2, sprintf('~/Desktop/discontPerm/%s/nx%i_dte-7.jpg', avg, N));
    %saveas(f3, sprintf('~/Desktop/discontPerm/%s/nx%i_time_dte-7.jpg',avg, N));
    %legend(p2, 'Arith', 'Cancel p_x', 'Cancel p_{xx}', 'Cancel both', 'Cancel p_{xxx}', 'Cancel p_{xxx}, p_{x}','Cancel p_{xxx}, p_{xx}', 'Cancel all', 'Location','Northeast')
    %legend(p3, 'Arith', 'Cancel p_x', 'Cancel p_{xx}', 'Cancel both', 'Cancel p_{xxx}', 'Cancel p_{xxx}, p_{x}','Cancel p_{xxx}, p_{xx}', 'Cancel all', 'Location','Southeast')
end
if (~conserv)
    if (strcmp(k_deriv, 'Central_2nd')) %No oscillations add antidiffusion to create them
        legend(p2, '2nd Order Central Difference k_x', 'Canceled Positive Diffusion truncation', 'Location','Southwest')
        legend(p3, '2nd Order Central Difference k_x', 'Canceled Positive Diffusion truncation', 'Location','SouthEast')
    elseif (strcmp(k_deriv, 'Central_4th')) %No oscillations add antidiffusion to create them
        legend(p2, '4th Order Central Difference k_x', 'Canceled Positive Diffusion truncation', 'Location','Southwest')
        legend(p3, '4th Order Central Difference k_x', 'Canceled Positive Diffusion truncation', 'Location','SouthEast')
    elseif(strcmp(k_deriv, 'D_plus'))
        legend(p2, 'D_+ k_x', 'Cancel p_x', 'Cancel p_{xx}','Cancel both', 'Location', 'NorthEast')
        legend(p3, 'D_+k_x', 'Cancel p_x', 'Cancel p_{xx}','Cancel both', 'Location', 'SouthEast')
    elseif(strcmp(k_deriv, 'D_minus'))
        legend(p2, 'D_- k_x', 'Cancel 1^{st} Order p_x', 'Cancel 1^{st} Order p_{xx}', 'Cancel both', 'Cancel both & 2^{nd} Order p_{xx}', 'Location','Southwest')
        legend(p3, 'D_- k_x', 'Cancel 1^{st} Order p_x', 'Cancel 1^{st} Order p_{xx}', 'Cancel both', 'Cancel both & 2^{nd} Order p_{xx}', 'Location','Southeast')
    elseif (strcmp(k_deriv, 'explicit'))
        legend(p2, 'explicit k_x = np^{n-1}p_x', 'Add positive diffusion on order dx^2', 'Location','Northeast')
        legend(p3, 'explicit k_x = np^{n-1}p_x', 'Add positive diffusion on order dx^2', 'Location','Southeast')
    end
    %saveas(f2, sprintf('~/Desktop/nonconservForm/%s/nx%i.jpg', k_deriv, N));
    %saveas(f3, sprintf('~/Desktop/nonconservForm/%s/nx%i_time.jpg', k_deriv, N));
elseif (strcmp(avg, 'fv harm'))
    if (k_pow == 3 || k_pow == 1 || k_pow == -1)
        %legend(p2, 'Harmonic', 'Cancel p_x', 'Cancel p_{xx}','Cancel both', 'Location','Northeast')
        %legend(p3, 'Harmonic', 'Cancel p_x', 'Cancel p_{xx}','Cancel both', 'Location','Southeast')
    elseif (k_pow == 2)
        %legend(p2, 'Harmonic', 'Cancel p_{xx}', 'Location','Northeast')
        %legend(p3, 'Harmonic', 'Cancel p_{xx}', 'Location','Southeast')
    end
    %saveas(f2, sprintf('~/Desktop/conservForm/%s/n=%i/nx%i.jpg', avg, k_pow,N));
    %saveas(f3, sprintf('~/Desktop/conservForm/%s/n=%i/nx%i_time.jpg',avg, k_pow,N));
elseif (strcmp(avg, 'arith')) %create oscillations
    if (k_pow == 3 || k_pow == -1)
        %legend(p2, 'Arithmetic', 'Cancel p_x', 'Cancel p_{xx}','Cancel both', 'Location','Northeast')
        %legend(p3, 'Arithmetic', 'Cancel p_x', 'Cancel p_{xx}','Cancel both', 'Location','Southeast')
%     elseif(k_pow == 2)
%         legend(p2, 'Arithmetic', 'Add anti-diff p_{xx}', 'Location','Northeast')
%         legend(p3, 'Arithmetic', 'Add anti-diff p_{xx}', 'Location','Southeast')
     end
    %saveas(f2, sprintf('~/Desktop/conservForm/%s/n=%i/nx%i.jpg', avg, k_pow,N));
    %saveas(f3, sprintf('~/Desktop/conservForm/%s/n=%i/nx%i_time.jpg',avg, k_pow,N));
end
end
toc;
gcf = 3;
%p_exact = load('discont');
x_exact = 0:1/3200:1;
%p2(count+1) = plot(t, p_exact.p_t_exact, 'Linewidth', 2);
%legend('show', 'N = 25', 'N = 50', 'N = 100', 'Exact Sol.')
%legend('location', 'southeast')
%saveas(gcf,sprintf('~/Desktop/Results/p_%d/upwindFlux/%s_%s_time.jpg', k_pow, JSolver, limiter))
%legend(p2, 'Ref. Sol.', 'Arithmetic', 'Harmonic','Mod. Harmonic', 'Location','Northeast')
%legend(p3, 'Ref. Sol.', 'Arithmetic', 'Harmonic','Mod. Harmonic', 'Location','Southeast')
%gcf = 2;
%x_exact = 0:1/3200:1;
%title('Pressure vs. position for harmonic averaging')
%title('Pressure vs. time for harmonic averaging')
%p2(count+1) = plot(x_exact, p_exact.p_exact, 'Linewidth', 2);
%legend(p2, 'Harmonic', 'Cancel p_x', 'Cancel p_{xx}','Cancel both', 'Ref. Sol', 'Location','Northeast')
 %saveas(gcf,sprintf('~/Desktop/Results/p/%s_%s_pos.jpg', JSolver, limiter))
%saveas(f2, sprintf('~/Desktop/discontPerm/conv_pos.jpg', avg, N));
%saveas(f3, sprintf('~/Desktop/discontPerm/conv_time.jpg', avg, N));
%saveas(f2, sprintf('~/Desktop/NewResults/new_global_arctan_gap0.001/lambda1_2/%s_pos_nx%i.jpg', avg, N));
%saveas(f3,
%sprintf('~/Desktop/NewResults/new_global_arctan_gap0.001/lambda1_2/%s_time
%_nx%i.jpg', avg, N));
%image on top of iamge!!
%axes('pos',[.1 .6 .5 .3])
%imshow('test.png')
if (upwindFlux)
    save(sprintf('data/foam/slopeLimiters/upwindFlux/%s%sN%i.mat',JSolver,limiter,N),'x_cent','p','p_t','t')
else
    save(sprintf('data/foam/slopeLimiters/downwindFlux/%s%s%i.mat',JSolver,limiter,N),'x_cent','p','p_t','t')
end
