close all
addpath('fixPSlinestyle')
N = 100;
nolabel = false;
addpath('data/foam/piecewiselinear_ic/')
pharm02 = load(sprintf('harm_nx%int02.mat',N));
parith02 = load(sprintf('arith_nx%int02.mat',N));
pint02 = load(sprintf('int_nx%int02.mat',N));
psam02 = load(sprintf('sam_nx200nt02.mat',N));
parith1 = load(sprintf('arith_nx%int1.mat',N));
pint1 = load(sprintf('int_nx%int1.mat',N));
psam1 = load(sprintf('sam_nx200nt1.mat',N));
pharm1 = load(sprintf('harm_nx%int1.mat',N));
figure(1)
x02 = parith02.x_cent;
t1 = parith1.t;
x1 = parith1.x_cent;
if (nolabel)
    p = plot(psam02.x_cent, psam02.p);
else
    p = plot(psam02.x_cent, psam02.p, x02, parith02.p, x02, pharm02.p, x02, pint02.p);
end
set(gca, 'FontSize',16,'FontWeight','bold');
set(findall(gca, 'Type', 'Line'),'LineWidth',4);
if (nolabel)
     set(p(1),'LineStyle','--','Color','r','Linewidth', 6,...
    'DisplayName','Ref. Sol.');
end
if (~nolabel)
    set(p(1),'LineStyle','-','Color','k','Linewidth', 12,...
    'DisplayName','Ref. Sol.');
    set(p(2),'Marker','o','LineStyle',':','Color',[1 0 0],...
        'DisplayName','Arithmetic');
    set(p(3),'Marker','^','LineStyle','--','Color',[0 0 1],...
        'DisplayName','Harmonic');
    set(p(4),'Marker','x','LineStyle','-.',...
        'Color',[0 0.498039215803146 0],...
        'DisplayName','Integral', 'Markersize', 6);
end
% set(p(3),'MarkerSize',6,'Marker','pentagram','LineWidth',3,...
%     'LineStyle','-.',...
%     'Color',[0.47843137383461 0.062745101749897 0.894117653369904],...
%     'DisplayName','Ref. Sol.');
if (~nolabel)
    legend('Ref. Sol.','Arithmetic', 'Harmonic','Integral')
end
xlabel('x')
ylabel('p')
if (~nolabel)
    title('Solution profile as a function of x')
end
hold on
if (nolabel)
    p = plot(psam1.x_cent, psam1.p);
else
    p = plot(psam1.x_cent, psam1.p, x1, parith1.p, x1, pharm1.p, x1, pint1.p);
end
%set(findall(gca, 'Type', 'Line'),'LineWidth',4);
set(gca, 'FontSize',16,'FontWeight','bold');
if (nolabel)
    set(p(1),'LineStyle',':','Color',[0 0.498039215803146 0],'Linewidth', 4,...
    'DisplayName','Ref. Sol.');
end
if (~nolabel)
    set(p(1),'LineStyle','-','Color','k','Linewidth', 12,...
    'DisplayName','Ref. Sol.');
    set(p(2),'Marker','o','LineStyle',':','Color',[1 0 0],'Linewidth', 4,...
        'DisplayName','Arithmetic');
    set(p(3),'Marker','^','LineStyle','--','Color',[0 0 1],'Linewidth', 4,...
        'DisplayName','Harmonic');
    set(p(4),'Marker','x','LineStyle','-.','Linewidth', 4,...
        'Color',[0 0.498039215803146 0],...
        'DisplayName','Integral', 'Markersize', 6);
end
if (nolabel)
%     annotation(gcf,'arrow',[0.687548942834769 0.756460454189507],...
%     [0.816725752508361 0.816053511705686],'LineWidth',2);
% 
% % Create textbox
% annotation(gcf,'textbox',...
%     [0.704210649960846 0.834448160535117 0.0342388410336725 0.0317725752508361],...
%     'String',{'t'},...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',16,...
%     'FitBoxToText','off',...
%     'LineStyle','none');
end
% set(p(3),'MarkerSize',6,'Marker','pentagram','LineWidth',3,...
%     'LineStyle','-.',...
%     'Color',[0.47843137383461 0.062745101749897 0.894117653369904],...
%     'DisplayName','SAM');
p = zeros(size(x1));
p(x1 <= 0.5) = 1 - 2*x1(x1 <= 0.5)';
p = plot(x1,p);
if (~nolabel)
set(p,'LineStyle','-','Color','k','Linewidth',6,...
    'DisplayName','Ref. Sol.');
else
    set(p,'LineStyle','-','Color','b','Linewidth',4,...
    'DisplayName','Ref. Sol.');
end
   
% set(p,'MarkerSize',6,'Marker','pentagram','LineWidth',3,...
%     'LineStyle','-.',...
%     'Color',[0.47843137383461 0.062745101749897 0.894117653369904],...
%     'DisplayName','SAM');
% Create textbox
if (~nolabel)
    annotation('textbox',...
        [0.554539450634153 0.39625850340136 0.104975474738982 0.0799319727891148],...
        'String',{'t = 0.1'},...
        'FontWeight','bold',...
        'FontSize',16,...
        'FitBoxToText','off',...
        'LineWidth',2);

    % Create textbox
    annotation('textbox',...
        [0.425373134328358 0.397959183673469 0.108208955223881 0.0820408163265285],...
        'String',{'t = 0.02'},...
        'FontWeight','bold',...
        'FontSize',16,...
        'FitBoxToText','off',...
        'LineWidth',2);

    % Create textbox
    annotation('textbox',...
        [0.216828358208955 0.409863945578231 0.0994029850746275 0.0790476190476162],...
        'String',{'t = 0.0'},...
        'FontWeight','bold',...
        'FontSize',16,...
        'FitBoxToText','off',...
        'LineWidth',2);
end
if (nolabel)
    set(gca,'Xtick',[])
    set(gca,'Ytick',[])
end
% p = plot(x1, p);
% set(p,'MarkerSize',6,'Marker','pentagram','LineWidth',3,...
%     'LineStyle','-.',...
%     'Color',[0.47843137383461 0.062745101749897 0.894117653369904],...
%     'DisplayName','SAM');
print(1,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/Discontinous/Results/';
if (nolabel)
    fixPSlinestyle('test.eps',[fileLoc sprintf('plinearic_posintro%ino.eps',N)]);
else
    %fixPSlinestyle('test.eps',[fileLoc sprintf('plinearic_posintro%i.eps',N)]);
end
figure(2)
t1 = t1 -0.0479;
if (~nolabel)
    p = plot(psam1.t-0.0479, psam1.p_t, t1, parith1.p_t,t1, pharm1.p_t,t1,pint1.p_t);
else
    p = plot(psam1.t-0.0479, psam1.p_t);
end
%add in exact solution for reference
set(findall(gca, 'Type', 'Line'),'LineWidth',4);
set(gca, 'FontSize',16,'FontWeight','bold');
% set(p(1),'LineStyle','--',...
%     'Color',[0.47843137383461 0.062745101749897 0.894117653369904],...
%     'DisplayName','Integral');
set(p(1),'LineStyle','-',...
     'Color','k',...
     'DisplayName','Ref. Sol.');
if (~nolabel)
    set(p(2),'LineStyle',':','Color',[1 0 0],...
        'DisplayName','Arithmetic');
    set(p(3),'LineStyle','--','Linewidth', 8, 'Color',[0 0 1],...
        'DisplayName','Harmonic');
    set(p(4),'LineStyle','-.','Color',[0 0.498039215803146 0],...
        'DisplayName','Integral');
end
xlabel('t')
ylabel('p')
x_l = 0.005;
x_u = 0.025;
y_l = 0.48;
y_u = 0.56;
h = y_u - y_l;
w = x_u - x_l;
r = rectangle('Position',[x_l y_l w h],'LineWidth',2);
title('Solution profile as a function of t')
legend('Ref. Sol', 'Arithmetic', 'Harmonic','Integral', 'Location', 'northwest')
print(2,'-depsc2','test');
fixPSlinestyle('test.eps',[fileLoc 'pclinic_intro.eps']);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim([0.48 0.56]);
% %box('on');
% set(findall(gca, 'Type', 'Line'),'LineWidth',6);
% xlabel('t')
% ylabel('p')
% %axes1 = axes('FontWeight','bold','FontSize',16);
% title('Solution profile as a function of t')
% print(2,'-depsc2','test');
% if (nolabel)
%     fixPSlinestyle('test.eps',[fileLoc 'pclinic_zoomintrono.eps']);
% else
%     fixPSlinestyle('test.eps',[fileLoc 'pclinic_zoomintro.eps']);
% end

