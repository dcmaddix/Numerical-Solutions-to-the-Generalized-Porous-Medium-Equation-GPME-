%generateFigures
close all
addpath('../fixPSlinestyle')
N = 100;
pi = load(sprintf('data/foam/sam_nx%i.mat',N));
ref = load(sprintf('data/foam/exact_nx%i.mat',N));
createfigure_sam(pi.x_cent, [ref.p_exact pi.p])
figure(1)
legend('Exact Sol.', 'SAM')
print(1,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/Discontinous/Results/';
fixPSlinestyle('test.eps',[fileLoc sprintf('pos_sam%i.eps',N)]);
createfigure_foamsamt((pi.t - 0.0479), [ref.p_t_exact' pi.p_t'])
legend('Exact Sol.', 'SAM', 'Location', 'southeast')
 y_l = 0.45;
 y_u = 0.6;
 x_l = 0.015;
 x_u = 0.05;
 h = y_u - y_l;
 w = x_u - x_l;
r = rectangle('Position',[x_l y_l w h],'LineWidth',2,...
    'EdgeColor',[0.235294118523598 0.235294118523598 0.235294118523598]);
print(2,'-depsc2','test');
fixPSlinestyle('test.eps',[fileLoc sprintf('time_sam%i.eps',N)]);
ylim(gca,[y_l y_u]);
xlim(gca,[x_l x_u]);
legend('Exact Sol.', 'SAM', 'Location', 'southeast')
set(r,'Visible','off')
print(2,'-depsc2','test');
fixPSlinestyle('test.eps',[fileLoc sprintf('time_sam%i_zoom.eps',N)]);