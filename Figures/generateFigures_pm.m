close all 
fileLoc = 'data/p^m/ref1600_spacem';
 p1 = load([fileLoc '1.mat']);
 p2 = load([fileLoc '2.mat']);
 p3 = load([fileLoc '3.mat']);
createfigurepm(p1.x_cent, [p1.p p2.p p3.p]);
print(1,'-depsc2','test');
fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m.eps');
 