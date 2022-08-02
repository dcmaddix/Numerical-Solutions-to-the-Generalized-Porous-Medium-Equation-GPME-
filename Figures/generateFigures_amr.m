%generateFigures
close all
%addpath('../fixPSlinestyle')
N = 50;
pa = load(sprintf('data/foam/int_nx%i.mat',N));
ph = load(sprintf('data/foam/intamr_nx%i.mat',N));
ref = load(sprintf('data/foam/exact_nx%i.mat',N));
createfigure_foamamr(pa.x_cent, [ref.p_exact pa.p ph.p])
legend('Exact Sol.', 'Integral', 'AMR Integral')
figure(1)
print(1,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/Discontinous/Results/';
fixPSlinestyle('test.eps',[fileLoc sprintf('pos_intamr%i.eps',N)]);
figure(2)
createfigure_foamamrt((pa.t - 0.0479), [ref.p_t_exact', pa.p_t'], (ph.t - 0.0479), ph.p_t')
legend('Exact Sol.', 'Integral', 'AMR Integral', 'Location', 'Southeast')
 y_l = 0.5;
 if (N == 50)
    y_u = 0.54;
 else
     y_u = 0.52;
 end
 x_l = 0.018;
 if (N == 50)
    x_u = 0.035;
 else
     x_u = 0.026;
 end
 h = y_u - y_l;
 w = x_u - x_l;
r = rectangle('Position',[x_l y_l w h],'LineWidth',2,...
    'EdgeColor',[0.235294118523598 0.235294118523598 0.235294118523598]);
print(3,'-depsc2','test');
fixPSlinestyle('test.eps',[fileLoc sprintf('time_intamr%i.eps',N)]);
ylim(gca,[y_l y_u]);
xlim(gca,[x_l x_u]);
set(r,'Visible','off')
legend('Exact Sol.', 'Integral', 'AMR Integral', 'Location', 'Southeast')
print(3,'-depsc2','test');
fixPSlinestyle('test.eps',[fileLoc sprintf('time_intamr%i_zoom.eps',N)]);