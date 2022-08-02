close all
N = 100;
fileLoc = sprintf('data/p^m/temp_oscill/nx%i/', N);
pt1 = load([fileLoc 'ptime1.mat']);
pt2 = load([fileLoc 'ptime2.mat']);
pt3 = load([fileLoc 'ptime3.mat']);
pt4 = load([fileLoc 'ptime4.mat']);
time = load([fileLoc 'time.mat']);
time = time.time;
x = pt1.x;
% plot(x, pt1.p, x, pt2.p, x,pt3.p, x,pt4.p);
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         set(gca, 'FontSize',16,'FontWeight','bold');
% xlabel('x')
% ylabel('p')
% title('Solution profile as a function of x')
createfigure_temposcill(x, [pt1.p pt2.p pt3.p pt4.p]);
if (N == 50)
 legend(sprintf('t = %.3f', time(1)), sprintf('t = %.3f', time(2)), ...
     sprintf('t = %.3f', time(3)),sprintf('t = %.3f', time(4)))
else
    legend(sprintf('t = %.4f', time(1)), sprintf('t = %.4f', time(2)), ...
     sprintf('t = %.4f', time(3)),sprintf('t = %.4f', time(4)))
end
axis([0.1 0.3 0 2])
print(1,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/temposcill/';
fixPSlinestyle('test.eps',[fileLoc sprintf('nx%i.eps', N)]);
