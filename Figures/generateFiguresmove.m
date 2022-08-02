close all
nolabel = true;
fileLoc = 'data/p^m/movinginterface/movinginterface_nx400';
 pt0 = load([fileLoc 't0.mat']);
 pt02 = load([fileLoc 't02.mat']);
 pt04 = load([fileLoc 't04.mat']);
 pt08 = load([fileLoc 't08.mat']);
 createfiguremove(pt0.x_cent, [pt0.p pt02.p pt04.p pt08.p], nolabel);
 print(1,'-depsc2','test');
 if (nolabel)
   fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/movinginterfaceno.eps');
 else
    fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/movinginterface.eps');
 end
 figure(2)
 plot1 = plot(pt0.x_cent,pt08.p,'r','LineWidth',4);
 % Create xlabel
xlabel('x','FontWeight','bold','FontSize',16);
% Create ylabel
ylabel('p','FontWeight','bold','FontSize',16);
set(gca, 'Xtick', [])
set(gca, 'Ytick', [])
print(2,'-depsc2','test');
 fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/intropme.eps')
 figure(2)
 plot1 = plot(pt0.x_cent,pt0.p,'LineWidth',4);
 % Create xlabel
xlabel('x','FontWeight','bold','FontSize',16);
% Create ylabel
ylabel('p','FontWeight','bold','FontSize',16);
set(gca, 'Xtick', [])
set(gca, 'Ytick', [])
print(2,'-depsc2','test');
fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/intro.eps');
pt08 = load('data/heat_eqtn.mat');
plot1 = plot(pt08.x,pt08.p,'k', 'LineWidth',4);
xlabel('x','FontWeight','bold','FontSize',16);
% Create ylabel
ylabel('p','FontWeight','bold','FontSize',16);
set(gca, 'Xtick', [])
set(gca, 'Ytick', [])
print(2,'-depsc2','test');
 fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/introheat.eps')
 
