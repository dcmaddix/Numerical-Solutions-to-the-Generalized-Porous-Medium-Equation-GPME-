%convergence study for parabolic problem for k = p^3 and arithmetic
%averaging
N = [50 100 200 400 800];
n_grids = length(N);
h = 1 ./ N;
l2_error = zeros(3, length(N));
% dir_harm = '../data/fv harm/p_3/harmonic_pressure_';
% dir_arith = '../data/arith/p_3/arithmetic_pressure_';
% dir_modharm = '../data/mod_harm/p_3/central_2nd/harmonic_pressure_';
% dir_harm = '../data/fv harm/p_2/harmonic_pressure_';
% dir_arith = '../data/arith/p_2/arithmetic_pressure_';
% dir_modharm = '../data/mod_harm/p_2/harmonic_pressure_';
% dir_harm = '../data/fv harm/p_1/harmonic_pressure_';
% dir_arith = '../data/arith/p_1/arithmetic_pressure_';
% dir_modharm = '../data/mod_harm/p_1/harmonic_pressure_';
dir_harm = '../data/exp_p/fv harm/harmonic_pressure_';
dir_arith = '../data/exp_p/arith/arithmetic_pressure_';
dir_modharm = '../data/exp_p/mod_harm/harmonic_pressure_';
%exact_harm = load(strcat(dir_harm, 'refsol_p1'));
%exact_arith = load(strcat(dir_arith, 'refsol_p1'));
exact_mod = load('../data/exp_p/arith/arithmetic_nx1600.mat');
exact_arith = exact_mod;
lb = 0.0; %for bdry heat eqtn test case
ub = 1.0; %for bdry heat eqtn test case
N_ex = 3200;
h_ex = 1/3200;
x = linspace(lb + h_ex/2, ub-h_ex/2, N_ex);
%define x at cell-centers at cell centers to plot against pressure
    x_cent = [lb 0.5 * (x(1:end-1) + x(2:end)) ub];
%exact_arith = load(strcat(dir_arith, 'refsol_p3'));
for i = 1:n_grids-1
    %store harmonic error
    num = load(strcat(dir_harm, sprintf('nx%i', N(i))));
    l2_error(1,i) = h(i)^(0.5)*norm(num.p - exact_arith.p(1:2^(n_grids+1-i):end));
    %l2_error(1,i) = norm(num.p - exact_arith.p(1:2^(n_grids+1-i):end), 'inf');
    l2_error(1,i) = h(i)^(0.5) * norm(num.p - exact_arith.p(1:2^(n_grids+1-i):end), 1);
    %store arith error
    num = load(strcat(dir_arith, sprintf('nx%i', N(i))));
    l2_error(2,i) = h(i)^(0.5)*norm(num.p - exact_arith.p(1:2^(n_grids+1-i):end));
    %l2_error(2,i) = norm(num.p - exact_arith.p(1:2^(n_grids+1-i):end), 'inf');
    l2_error(2,i) = h(i)^(0.5) * norm(num.p - exact_arith.p(1:2^(n_grids+1-i):end), 1);
    num = load(strcat(dir_modharm, sprintf('nx%i', N(i))));
    %need to compare against new arithmetic reference solution since both
    %had different initial conditions!!!!
    %depends +-1
    l2_error(3,i) = h(i)^(0.5)*norm(num.p - exact_mod.p(1:2^(n_grids+1-i):end));
    %l2_error(3,i) = norm(num.p - exact_mod.p(1:2^(n_grids+1-i):end), 'inf')
    l2_error(3,i) = h(i)^(0.5) * norm(num.p - exact_mod.p(1:2^(n_grids+1-i):end), 1)
    
    
    %define spatial vector at cell faces to plot against velocity
    
    %define x at cell-centers at cell centers to plot against pressure
    %p_e(ind)
    %num.p(ind)
    %l2_error(2,i) = h(i)^(0.5)*norm(num.p - exact_arith.p_exact(1:2^(n_grids+1-i):end));
    %linear loglog plot in asymtotic convergence region
end
%legend(plots, 'N = 25', 'N = 50', 'N = 100', 'N = 200', 'N = 400', 'N = 800')
conv_harm_p3 = polyfit(log10(h(1:end-1)), log10(l2_error(1,1:end-1)),1);
conv_harm_p3(1)
conv_arith_p3 = polyfit(log10(h(1:end-1)), log10(l2_error(2,1:end-1)),1);
conv_arith_p3(1)
conv_modharm_p3 = polyfit(log10(h(1:end-1)), log10(l2_error(3,1:end-1)),1);
conv_modharm_p3(1)
second_order = 0.7 ./ 4.^(0:length(N)-1); %0.5*(l2_error(1,1) + l2_error(2,1)) ./ 4.^(0:length(N)-1);
sec_order = polyfit(log10(h), log10(second_order),1);
first_order = 0.1 ./ 2.^(0:length(N)-1);
sec_order(1)
figure(2)
%  plot(log10(h(3:end)), log10(second_order(3:end)),'-ko', log10(h(3:end)),log10(l2_error(1,3:end)), '-cd',...,
%      log10(h(3:end)), log10(l2_error(2,3:end)),'-rs', log10(h(3:end)), log10(l2_error(3,3:end)), '-o')%, ...
     %log10(h(3:end)), log10(first_order(3:end)),'-o')
plot(...     
    log10(h(1:end-1)), log10(l2_error(2,1:end-1)),'-s', ...
    log10(h(1:end-1)), log10(l2_error(3,1:end-1)), '-ro')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca, 'FontSize',16,'FontWeight','bold');
%legend('Second Order', 'Harmonic', 'Arithmetic' ,'Modified Harmonic', 'Location', 'southeast')
legend('Arithmetic' ,'Modified Harmonic', 'Location', 'southeast')
xlabel('log_{10}(\Deltax)')
ylabel('log_{10}(Error)')
%title('l_2 norm of error for k = p^2 ')
%title('Max. norm of error for k = p ')
title('l_1 norm of error for k = exp(-1/p) ')
%axis([-3.1 -1.8 -7 -2])
%axis([-3.5 -1.8 -6.5 -1])
%axis([-3.5 -1.8 -5 -0.5])
%axis([-3.5 -1.8 -5 0.1])
%axis([-3.5 -1.8 -4.5 1])
axis([-3 -1.5 -3.5 -0.25])
print(2,'-depsc2','test');
%fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/exp_p/loglog_l1.eps');
