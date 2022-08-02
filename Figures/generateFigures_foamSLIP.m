%generateFigures
close all
addpath('fixPSlinestyle')
N = 50;
SLIP = false;
lim = true;
if (SLIP)
    pimm = load(sprintf('data/foam/slopeLimiters/upwindFlux/SLIPminmodN%i.mat',N));
    pisb = load(sprintf('data/foam/slopeLimiters/upwindFlux/SLIPsuperbeeN%i.mat',N));
    pivl = load(sprintf('data/foam/slopeLimiters/upwindFlux/SLIPvanleerN%i.mat',N));
elseif(lim)
    pimm = load(sprintf('data/foam/slopeLimiters/upwindFlux/JSTN%i.mat',N));
    pisb = load(sprintf('data/foam/slopeLimiters/upwindFlux/SLIPvanleerN%i.mat',N));
    pivl = load(sprintf('data/foam/slopeLimiters/upwindFlux/USLIPvanleerN%i.mat',N));
else
    pimm = load(sprintf('data/foam/slopeLimiters/upwindFlux/USLIPminmodN%i.mat',N));
    pisb = load(sprintf('data/foam/slopeLimiters/upwindFlux/USLIPsuperbeeN%i.mat',N));
    pivl = load(sprintf('data/foam/slopeLimiters/upwindFlux/USLIPvanleerN%i.mat',N));
end
ref = load(sprintf('data/foam/exact_nx%i.mat',N));
createfigure_SLIP(pimm.x_cent, [ref.p_exact pimm.p pisb.p pivl.p])
if (lim)
    legend('Exact Sol.', 'JST', 'SLIP van leer', 'USLIP van leer')
end
figure(1)
print(1,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/Discontinous/Results/';
if (SLIP)
    fixPSlinestyle('test.eps',[fileLoc sprintf('pos_slip%i.eps',N)]);
elseif (lim)
    fixPSlinestyle('test.eps',[fileLoc sprintf('pos_lim%i.eps',N)]);
else
    fixPSlinestyle('test.eps',[fileLoc sprintf('pos_uslip%i.eps',N)]);
end
createfigure_foamjstt((pimm.t - 0.0479), [ref.p_t_exact' pimm.p_t' pisb.p_t' pivl.p_t'])
if (SLIP)
    legend('Exact Sol.', 'minmod', 'superbee', 'van leer', 'Location', 'Southeast')
else
    legend('Exact Sol.', 'JST', 'SLIP van leer', 'USLIP van leer', 'Location', 'Southeast')
end
 y_l = 0.45;
 y_u = 0.6;
 x_l = 0.015;
 x_u = 0.05;
 h = y_u - y_l;
 w = x_u - x_l;
r = rectangle('Position',[x_l y_l w h],'LineWidth',2,...
    'EdgeColor',[0.235294118523598 0.235294118523598 0.235294118523598]);
print(2,'-depsc2','test');
if (SLIP)
    fixPSlinestyle('test.eps',[fileLoc sprintf('time_slip%i.eps',N)]);
elseif (lim)
    fixPSlinestyle('test.eps',[fileLoc sprintf('time_lim%i.eps',N)]);
else
    fixPSlinestyle('test.eps',[fileLoc sprintf('time_uslip%i.eps',N)]);
end
ylim(gca,[y_l y_u]);
xlim(gca,[x_l x_u]);
set(r,'Visible','off')
if (SLIP)
    legend('Exact Sol.', 'minmod', 'superbee', 'van leer', 'Location', 'Southeast')
else
    legend('Exact Sol.', 'JST', 'SLIP van leer', 'USLIP van leer', 'Location', 'Southeast')
end
print(2,'-depsc2','test');
if (SLIP)
    fixPSlinestyle('test.eps',[fileLoc sprintf('time_slip%i_zoom.eps',N)]);
elseif (lim)
    fixPSlinestyle('test.eps',[fileLoc sprintf('time_lim%i_zoom.eps',N)]);
else
    fixPSlinestyle('test.eps',[fileLoc sprintf('time_uslip%i_zoom.eps',N)]);
end