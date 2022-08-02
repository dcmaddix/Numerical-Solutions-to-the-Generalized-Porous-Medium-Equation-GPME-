%%Plot k_{eff}
close all
eps = 1e-4;
dy = 1e-6;
y = 0:dy:1;
harm = (2*eps)/(1+eps);
f = (eps) ./ ((eps -1)*y+1);
ind_H = find(f == harm);
arith = (1+eps)/2;
ind_A = find(abs(arith - f) < 1e-4);
plot(y,f, '-o');
hold on;
plot(y(1), f(1), '-or',y(end), f(end), '-og', y(ind_H), f(ind_H), '-om', y(ind_A), f(ind_A), '-oc','Linewidth',5,'Markersize',15)
legend('k_{eff}', 'k_{min}','k_{max}', 'k_{harm}', 'k_{arith}', 'Location', 'Northwest')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',16,'FontWeight','bold');
eps_string = sprintf('for k_{min} = 0.01',eps);
title(['Plot of k_{eff} vs. \Delta x* / \Delta x '])
xlabel('\Delta x* / \Delta x')
ylabel('k_{eff}')
axis([0 1 0 1.01])
%createfigure_keff(X1, Y1, X2, Y2, Y3, X3, Y4, X4, Y5)
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/Discontinous/Results/';
createfigure_keff(y, f, y(1), f(1), f(end), y(ind_H), f(ind_H), y(ind_A), f(ind_A), y(end))
%print(2,'-depsc2','test');
%fixPSlinestyle('test.eps',[fileLoc 'keff.eps']);
