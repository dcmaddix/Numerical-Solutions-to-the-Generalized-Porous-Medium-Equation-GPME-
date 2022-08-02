close all
addpath('../fixPSlinestyle')
%N = 100
px4 = load('../data/p^m/mh_px4100.mat');
pxx = load('../data/p^m/mh_pxx100.mat');
pmod = load('../data/p^m/mod_harm/p_3/harmonic_pressure_nx100.mat');
ph = load('../data/p^m/fv harm/p_3/new/harmonic_pressure_nx100.mat');
ref = load('../data/p^m/ref1600_spacem3.mat');
createfigure_pxpxx(ref.x_cent, ref.p, pmod.x_cent, [ph.p pxx.p px4.p pmod.p])
figure(1)
print(1,'-depsc2','test');
%fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/mhm100_px.eps');
ref_t = load('../data/p^m/ref1600_timem3.mat');
createfigure_timepxpxx(ref_t.t, ref_t.p_t, pmod.t, [ph.p_t' pxx.p_t' px4.p_t' pmod.p_t'])
print(2,'-depsc2','test');
%fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/mhm100_pxt.eps');