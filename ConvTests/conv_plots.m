%create error plot using max norm err and 2 norm err_norm
close all
d = load('arith_conv_data');
d_h = load('harm_conv_data');

n_grids = 5;
dx = 1 ./(2.^(0:n_grids)*100);
size(dx)
%% lambda = 1/2
f1 = figure(1);
ax1 = axes();
set(ax1,'fontsize',12, 'fontname' ,'times')
p1 = loglog(dx, d.err_lambda(1,:), dx, d_h.err_lambda(1,:), 'Linewidth', 2);
set(p1,'marker','o','markerfacecolor','k','markeredgecolor','k');
ylabel('l_2 Norm Error', 'FontSize' ,14)
xlabel('\Delta x', 'FontSize' ,14)
%set('interpreter', latex')
%set(f1, 'PaperPositionMode', 'auto')
title('Convergence plot of l_2 Norm Error vs \Delta x for \lambda = 1/2')
saveas(f1, '~/Desktop/NewResults/Conv/l2_norm_conv_lambda_1_2.jpg');
legend('Arithmetic Averaging', 'Harmonic Averaging', 'Location', 'Northwest')
p = polyfit(log10(dx), log10(d.err_lambda(1,:)),1);
convRate_arith = p(1)
p = polyfit(log10(dx), log10(d_h.err_lambda(1,:)),1);
convRate_harm = p(1)

%% lambda = 1/4
f2 = figure(2);
ax2 = axes();
dx = dx(1:end-1);
set(ax2,'fontsize',12, 'fontname' ,'times')
p2 = loglog(dx, d.err_lambda(2,1:end-1), dx, d_h.err_lambda(2,1:end-1), 'Linewidth', 2);
set(p2,'marker','o','markerfacecolor','k','markeredgecolor','k');
ylabel('l_2 Norm Error', 'FontSize' ,14)
xlabel('\Delta x', 'FontSize' ,14)
%set('interpreter', latex')
%set(f1, 'PaperPositionMode', 'auto')
title('Convergence plot of l_2 Norm Error vs \Delta x for \lambda = 1/4')
saveas(f2, '~/Desktop/NewResults/Conv/l2_norm_conv_lambda_1_4.jpg');
legend('Arithmetic Averaging', 'Harmonic Averaging', 'Location', 'Northwest')
p = polyfit(log10(dx), log10(d.err_lambda(2,1:end-1)),1);
convRate_arith = p(1)
p = polyfit(log10(dx), log10(d_h.err_lambda(2,1:end-1)),1);
convRate_harm = p(1)

%% lambda = 1/8
f3 = figure(3);
ax3 = axes();
set(ax3,'fontsize',12, 'fontname' ,'times')
p3 = loglog(dx, d.err_lambda(3,1:end-1), dx, d_h.err_lambda(3,1:end-1), 'Linewidth', 2);
set(p3,'marker','o','markerfacecolor','k','markeredgecolor','k');
ylabel('l_2 Norm Error', 'FontSize' ,14)
xlabel('\Delta x', 'FontSize' ,14)
%set('interpreter', latex')
%set(f1, 'PaperPositionMode', 'auto')
title('Convergence plot of l_2 Norm Error vs \Delta x for \lambda = 1/8')
saveas(f3, '~/Desktop/NewResults/Conv/l2_norm_conv_lambda_1_8.jpg');
legend('Arithmetic Averaging', 'Harmonic Averaging', 'Location', 'Northwest')
p = polyfit(log10(dx), log10(d.err_lambda(3,1:end-1)),1);
convRate_arith = p(1)
p = polyfit(log10(dx), log10(d_h.err_lambda(3,1:end-1)),1);
convRate_harm = p(1)