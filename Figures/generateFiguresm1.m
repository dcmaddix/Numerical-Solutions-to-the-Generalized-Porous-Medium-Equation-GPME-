%generateFigures
close all
addpath('../fixPSlinestyle')
%N = 50
pa = load('../data/p^m/arith/p_1/new/arithmetic_pressure_nx50.mat');
ph = load('../data/p^m/fv harm/p_1/new/harmonic_pressure_nx50.mat');
pmod = load('../data/p^m/mod_harm/p_1/harmonic_pressure_nx50.mat');
ref = load('../data/p^m/ref1600_spacem1.mat');
createfigurepos2(ref.x_cent, ref.p, pmod.x_cent, [pa.p ph.p pmod.p])
figure(1)
print(1,'-depsc2','test');
fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=1/pos50.eps');
ref_t = load('../data/p^m/ref1600_timem1.mat');
createfiguretime(ref_t.t, ref_t.p_t, pmod.t, [pa.p_t' ph.p_t' pmod.p_t'])
legend('Ref. Sol.', 'Arithmetic', 'Harmonic', 'Mod. Harmonic', 'Location', 'southeast')
print(2,'-depsc2','test');
fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=1/time50.eps');