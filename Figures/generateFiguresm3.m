%generateFigures
close all
addpath('../fixPSlinestyle')
%N = 50
pa = load('../data/p^m/arith/p_3/new/arithmetic_pressure_nx50.mat');
ph = load('../data/p^m/fv harm/p_3/new/harmonic_pressure_nx50.mat');
pmod = load('../data/p^m/mod_harm/p_3/harmonic_pressure_nx50.mat');
ref = load('../data/p^m/ref1600_spacem3.mat');
createfigurepos2(ref.x_cent, ref.p, pmod.x_cent, [pa.p ph.p pmod.p])
figure(1)
print(1,'-depsc2','test');
%fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/pos50.eps');
ref_t = load('../data/p^m/ref1600_timem3.mat');
createfiguretime(ref_t.t, ref_t.p_t, pmod.t, [pa.p_t' ph.p_t' pmod.p_t'])
legend('Ref. Sol.', 'Arithmetic', 'Harmonic', 'Mod. Harmonic', 'Location', 'southeast')
print(2,'-depsc2','test');
%fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/time50.eps');
%N = 100
pa = load('../data/p^m/arith/p_3/new/arithmetic_pressure_nx100.mat');
ph = load('../data/p^m/fv harm/p_3/new/harmonic_pressure_nx100.mat');
pmod = load('../data/p^m/mod_harm/p_3/harmonic_pressure_nx100.mat');
createfigurepos2(ref.x_cent, ref.p, pmod.x_cent, [pa.p ph.p pmod.p])
figure(1)
print(3,'-depsc2','test');
%fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/pos100.eps');
createfiguretime(ref_t.t, ref_t.p_t, pmod.t, [pa.p_t' ph.p_t' pmod.p_t'])
legend('Ref. Sol.', 'Arithmetic', 'Harmonic', 'Mod. Harmonic', 'Location', 'southeast')
print(4,'-depsc2','test');
%fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/time100.eps');
%N = 200
pa = load('../data/p^m/arith/p_3/new/arithmetic_pressure_nx200.mat');
ph = load('../data/p^m/fv harm/p_3/new/harmonic_pressure_nx200.mat');
pmod = load('../data/p^m/mod_harm/p_3/harmonic_pressure_nx200.mat');
createfigurepos2(ref.x_cent, ref.p, pmod.x_cent, [pa.p ph.p pmod.p])
figure(1)
print(5,'-depsc2','test');
%fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/pos200.eps');
createfiguretime(ref_t.t, ref_t.p_t, pmod.t, [pa.p_t' ph.p_t' pmod.p_t'])
legend('Ref. Sol.', 'Arithmetic', 'Harmonic', 'Mod. Harmonic', 'Location', 'southeast')
print(6,'-depsc2','test');
%fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/time200.eps');