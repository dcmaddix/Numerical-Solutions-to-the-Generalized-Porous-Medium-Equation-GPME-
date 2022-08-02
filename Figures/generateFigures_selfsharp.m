close all 
nolabel = true;
fileLoc = 'data/p^m/selfsharpeningm3/selfsharpening_nx100';
 pt0 = load([fileLoc 't0.mat']);
 pt02 = load([fileLoc 't02.mat']);
 pt04 = load([fileLoc 't04.mat']);
 pt08 = load([fileLoc 't08.mat']);
 createfigureselfsharp(pt0.x_cent, [pt0.p pt02.p pt04.p pt08.p], nolabel);
 print(1,'-depsc2','test');
 if (nolabel)
     fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/selfsharpeningno.eps');
 else
    fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/selfsharpening.eps');
 end
 