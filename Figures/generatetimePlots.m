%generate time plots
close all
N = 50;
figure(1)
fileLoc = 'data/p^m/fv harm/timeSchemes/';
BE_sqr = load([fileLoc sprintf('BEsqr_%i.mat', N)]);
BE = load([fileLoc sprintf('BE_%i.mat', N)]);
TVD = load([fileLoc sprintf('TVD_%i.mat', N)]);
fileLoc = 'data/p^m/fv harm/p_3/new/';
FE = load([fileLoc sprintf('harmonic_pressure_nx%i.mat', N)]);
p = load('data/p^m/ref1600_spacem3.mat');
p_t = load('data/p^m/ref1600_timem3.mat');
createfigure_timeCompt(p_t.t, p_t.p_t, FE.t, [FE.p_t' BE_sqr.p_t' TVD.p_t'], BE.t, BE.p_t)
print(2,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/timeSchemes/';
fixPSlinestyle('test.eps',[fileLoc sprintf('time%i.eps',N)]);
p_t = load('data/p^m/ref1600_spacem3.mat');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/timeSchemes';
%fixPSlinestyle('test.eps',[fileLoc sprintf('time%i.eps',N)]);
figure(2)
createfigure_timeComp(p.x_cent, p.p, FE.x_cent, [FE.p BE_sqr.p BE.p TVD.p])
print(3,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/timeSchemes/';
fixPSlinestyle('test.eps',[fileLoc sprintf('space%i.eps',N)]);