%generate plots
figure(1);
N = 50;
smallfont = true;
fileLoc = 'data/p^m/mimeticp3/';
pa3 = load([fileLoc 'p50_arith3.mat']);
pa5 = load([fileLoc 'p50_arith5.mat']);
pa7 = load([fileLoc 'p50_arith7.mat']);%* in mimetic
if (smallfont)
    p = plot(pa3.x_cent,pa3.p, 'o', pa5.x_cent,pa5.p,'rsquare', pa7.x_cent, pa7.p, '^')
    set(p(3),'Color',[0 0.498039215803146 0]);
else
    plot(pa3.x_cent,pa3.p, 'o', pa5.x_cent,pa5.p,'square', pa7.x_cent, pa7.p, '^')
end
set(findall(gca, 'Type', 'Line'),'LineWidth',4);
if (smallfont)
    set(gca, 'FontSize',16,'FontWeight','bold');
else
    set(gca, 'FontSize',20,'FontWeight','bold');
end
%title('Plot of numerical solution at various times')
title('Solution profile as a function of x')
legend('t = 0.3', 't = 0.5', 't = 0.7')
xlabel('x')
ylabel('p') 
figure(2);
pa3 = load([fileLoc 'p50_harm3.mat']);
pa5 = load([fileLoc 'p50_harm5.mat']);
pa7 = load([fileLoc 'p50_harm7.mat']);
if (smallfont)
    p = plot(pa3.x_cent,pa3.p, 'o', pa5.x_cent,pa5.p,'rsquare', pa7.x_cent, pa7.p, '^');
    set(p(3),'Color',[0 0.498039215803146 0]);
else
    plot(pa3.x_cent,pa3.p, 'o', pa5.x_cent,pa5.p,'square', pa7.x_cent, pa7.p, '^')
end
set(findall(gca, 'Type', 'Line'),'LineWidth',4);
if (smallfont)
    set(gca, 'FontSize',16,'FontWeight','bold');
else
    set(gca, 'FontSize',20,'FontWeight','bold');
end
%title('Plot of numerical solution at various times')
title('Solution profile as a function of x')
legend('t = 0.3', 't = 0.5', 't = 0.7')
xlabel('x')
ylabel('p') 
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/mimetic/';
print(1,'-depsc2','test');
if (smallfont)
    fixPSlinestyle('test.eps',[fileLoc 'arith50small.eps']);
else
    fixPSlinestyle('test.eps',[fileLoc 'arith50.eps']);
end
print(2,'-depsc2','test');
if (smallfont)
    fixPSlinestyle('test.eps',[fileLoc 'harm50small.eps']);
else
    fixPSlinestyle('test.eps',[fileLoc 'harm50.eps']);
end
figure(3);
[x_cent p_exact] = gentrueSolMim(N);
if (smallfont)
    p = plot(x_cent,p_exact(1,:), 'o', x_cent,p_exact(2,:),'rsquare', x_cent, p_exact(3,:), '^');
    set(p(3),'Color',[0 0.498039215803146 0]);
else
    plot(x_cent,p_exact(1,:), 'o', x_cent,p_exact(2,:),'rsquare', x_cent, p_exact(3,:), '^')
end
set(findall(gca, 'Type', 'Line'),'LineWidth',4);
if (smallfont)
    set(gca, 'FontSize',16,'FontWeight','bold');
else
    set(gca, 'FontSize',20,'FontWeight','bold');
end
title('Solution profile as a function of x')
legend('t = 0.3', 't = 0.5', 't = 0.7')
xlabel('x')
ylabel('p') 
%fileLoc = 'data/p^m/mimeticp3/';
print(3,'-depsc2','test');
if (smallfont)
fixPSlinestyle('test.eps',[fileLoc 'exact50small.eps']);
else
    fixPSlinestyle('test.eps',[fileLoc 'exact50.eps']);
end