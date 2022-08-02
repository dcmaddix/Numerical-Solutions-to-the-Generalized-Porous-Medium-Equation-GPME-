%generateFiguresexp_p
close all
addpath('../fixPSlinestyle')
N = 400;
pa = load(sprintf('../data/exp_p/arith/arithmetic_pressure_nx%i.mat',N));
ph = load(sprintf('../data/exp_p/fv harm/harmonic_pressure_nx%i.mat',N));
pmod = load(sprintf('../data/exp_p/mod_harm/harmonic_pressure_nx%i.mat',N));
ref = load('../data/exp_p/arith/arithmetic_nx1600.mat');
createfigurepos2(0:1/1600:1, ref.p, pmod.x_cent, [pa.p ph.p pmod.p])
figure(1)
print(1,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/exp_p/';
fixPSlinestyle('test.eps',[fileLoc sprintf('pos%i.eps',N)]);
createfiguretime(ref.t, ref.p_t, pmod.t, [pa.p_t' ph.p_t' pmod.p_t'])
legend('Ref. Sol.', 'Arithmetic', 'Harmonic', 'Mod. Harmonic', 'Location', 'southeast')
print(2,'-depsc2','test');
fixPSlinestyle('test.eps',[fileLoc sprintf('time%i.eps',N)]);