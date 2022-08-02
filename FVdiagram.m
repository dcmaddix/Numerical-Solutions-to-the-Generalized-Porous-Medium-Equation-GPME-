%create finite volume block
close all
%line([0 0], [0 1], 'Linewidth', 2, 'Color', 'k')
%line([1 1], [0 1], 'Linewidth', 2, 'Color', 'k')
%line([0.5 0.5], [0 1], 'Linewidth', 2, 'Color', 'k')
%set(gca,'XTickLabel','')
%set(gca,'YTickLabel','')
%axis off
%months = ['i';'j'];
%set(gca, 'XtickLabel', 'i, i+1')
%set(gca,'XTick',1:2,'XTickLabel',months)
%hold on;
%plot(0.25, 0.5, '.', 'Markersize', 30)
%plot(0.75, 0.5, '.', 'Markersize', 30)
%only create 1 image
%create horizontal line in middl eand dotted lines around star
%createfigure_FVnew([0 0], [0 1], [1 1], [0.5 0.5], 0.25, 0.5, 0.75, 0.6)
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/Discontinous/Results/';
%createfigure_FVfinal([0 0], [0 1], [1 1], [0.5 0.5], 0.25, 0.5, 0.75, 0.6)
fvgrid([0 0], [0 1], [1 1], [0.5 0.5], 0.25, 0.5, 0.75, 0.6)
%createfigure_FVnew([0 0], [0 1], [1 1], [0.5 0.5], 0.25, 0.5, 0.75, 0.4)
print(1,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPresentations/Results/';
%fixPSlinestyle('test.eps',[fileLoc 'FVdiagram_new.eps']);
fixPSlinestyle('test.eps',[fileLoc 'fvgrid2.eps']);
%set(gca,'XTick',[0.25 0.75],'XTickLabel',['i'; 'i+1'])
%neworder = {
%    'i + 1'  [0.600]};
%set(gca,'XtickL',neworder(:,1))
%set(gca,'XTick',0.75,'XTickLabel','j')
%str = 'Straight Line Plot from 1 to 10';
%dim = [.1 .5 0 .3];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');