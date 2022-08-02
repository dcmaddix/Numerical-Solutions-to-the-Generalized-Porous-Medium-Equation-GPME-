%generate plots for the test case showing harmonic is better
N = 25;
fileLoc = 'data/p^m/';
parith = load([fileLoc sprintf('arith/p_3/arithmetic_pressure_nx%i.mat', N)]);
pharm = load([fileLoc sprintf('fv harm/p_3/harmonic_pressure_nx%i.mat',N)]);
p = load([fileLoc 'ref1600_spacem3.mat']);
%stored at four times give it at a couple on same plot
x = 0:1/N:1;
xexact = 0:1/1600:1;
plot1 = plot(xexact, p.p, '-k', x, parith.p, '-r', x, pharm.p, '-b');
legend('Ref. Sol.', 'Arithmetic', 'Harmonic', 'Location', 'Northeast')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',16,'FontWeight','bold');
title('Solution profile as a function of x')
xlabel('x')
ylabel('p')
set(plot1(1), 'MarkerSize',8,'LineWidth',4)
set(plot1(2),'Marker','o','LineStyle','-.','Color',[1 0 0],...
'DisplayName','Arithmetic', 'Markersize', 8, 'Linewidth', 2);
set(plot1(3),'Marker','^','LineStyle','--','Color',[0 0 1],...
'DisplayName','Harmonic', 'Markersize', 8, 'Linewidth', 2);
set(gca, 'FontSize',16,'FontWeight','bold');
print(1,'-depsc2','test');
fileLoc2 = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/multi-material/';
fixPSlinestyle('test.eps',[fileLoc2 'pos_pme.eps']);
%load corresponding time plots as well
figure(2)
parith = load([fileLoc sprintf('arith/p_3/p_tarith%i.mat', N)]);
pharm = load([fileLoc sprintf('fv harm/p_3/p_tharm%i.mat',N)]);
p = load([fileLoc 'ref1600_timem3.mat']);
%stored at four times give it at a couple on same plot
x = 0:1/N:1;
xexact = 0:1/1600:1;
plot1 = plot(p.t, p.p_t, '-k', parith.t, parith.p_t, '-r', parith.t, pharm.p_t, '-b');
legend('Ref. Sol.', 'Arithmetic', 'Harmonic', 'Location', 'Southeast')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(plot1(1), 'LineWidth',4);
set(plot1(2),'LineStyle','-.','Color',[1 0 0],'DisplayName','Arithmetic', 'LineWidth',4);
set(plot1(3),'LineStyle','--','Color',[0 0 1],...
    'DisplayName','Harmonic', 'LineWidth',4);
%set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',16,'FontWeight','bold');
title('Solution profile as a function of t')
xlabel('t')
ylabel('p')
print(2,'-depsc2','test');
fileLoc2 = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/multi-material/';
fixPSlinestyle('test.eps',[fileLoc2 'pos_pme_time.eps']);