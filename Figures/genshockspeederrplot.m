%generate shock speed eror plot
%load('xierr_longtimedx04.mat')
N = 25;
%N = 12.5;
dx = 1/ N;
if (N == 25)
    load('foam/xierr_longtimedx04.mat')
elseif(N == 12.5)
    load('foam/xierr_longtimedx08.mat')
end
addpath('fixPSlinestyle/')
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/Discontinous/Results/';
plot(xdata, ydata / dx)
set(findall(gca, 'Type', 'Line'),'LineWidth',4);
set(gca, 'FontSize',16,'FontWeight','bold');
title('Shock Position Error Relative to \Deltax as a function of time')
ylabel('(\xi - x^*) / \Delta x')
xlabel('t')
print(1,'-depsc2','test');
fixPSlinestyle('test.eps',[fileLoc 'shockpos_err.eps']);