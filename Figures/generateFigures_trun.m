close all
addpath('fixPSlinestyle')
%N = 100
m = 3;
trun_arith = load(sprintf('data/trun/trun_arith_m%iNx100.mat', m));
trun_harm = load(sprintf('data/trun/trun_harm_m%iNx100.mat', m));
createfigure_trun(trun_arith.x_cent(2:end-1), [trun_arith.trun trun_harm.trun])
figure(1)
ylabel('k + (A+B^*)p_x^2')
title('Coefficient of p_{xx} as a function of x')
legend('Arithmetic (B^* = B)', 'Harmonic (B^* = B^H)')
print(1,'-depsc2','test');
fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/trunm%i.eps',m));