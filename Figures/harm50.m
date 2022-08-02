close all
N = 100;
fileLoc = 'data/p^m/fv harm/p_3/new/';
p = load([fileLoc sprintf('harmonic_pressure_nx%i.mat',N)]);
createfigure_harmtemposcill(p.t, p.p_t)
title('Solution profile as a function of t')
xlabel('t')
ylabel('p')
if (N == 50)
    axis([0.014 0.038 0.8 2])
else
    axis([0.0067 0.013 0.8 2])
    %xlim(gca,[0.0067 0.013]);
    %set(gca, 'XGrid', 'on')
    %axes('XTick',[0.0067 0.0083 0.0098 0.0114 0.0130]);
    set(gca, 'XTick', [0.0067 0.0083 0.0098 0.0114 0.0130]);
end
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',16,'FontWeight','bold');
print(1,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/temposcill/';
fixPSlinestyle('test.eps',[fileLoc sprintf('harm%i.eps', N)]);
%axis([0.0067 0.013 0.8 2])
