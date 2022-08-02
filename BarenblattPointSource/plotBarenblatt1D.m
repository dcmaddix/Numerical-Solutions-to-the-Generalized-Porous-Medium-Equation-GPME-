%script to plot 1D barenblatt solution
close all
m = [1 2 3 10 20 50 100];
m = 3;
radial = false;
bound = 5;
x = -bound:0.01:bound; %interface with positive region and compact support
t = 1.0;
t = 0.08;
d = 1; %1D solution
U = zeros(length(x), length(m));
for i = 1:length(m)
    U(:,i) = barenblatt_pme(m(i),t,d,x, x,radial); %defined in 
end
U_heat = max(U)*(4*pi*t)^(-1/2)*exp(-abs(x).^2/(4*t));
figure(1)
plot(x,U, '-o')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',16,'FontWeight','bold');
xlabel('x')
ylabel('p')
title('Solution profile as a function of x')
figure(2)
plot(x, U_heat, '-o')
%legend('m = 1', 'm = 2', 'm = 3', 'm = 10', 'm = 20', 'm = 50', 'm = 100')
%title('Self-similar Barenblatt-Prattle Solution to PME for various m')
xlabel('x')
ylabel('p')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',16,'FontWeight','bold');
xlabel('x')
ylabel('p')
title('Solution profile as a function of x')
figure(3)
plot(x, U, '-b',x,U_heat,'--r')
legend('PME', 'HE')
set(findall(gca, 'Type', 'Line'),'LineWidth',5);
set(gca, 'FontSize',16,'FontWeight','bold');
xlabel('x')
ylabel('p')
title('Solution profile as a function of x')
print(3,'-depsc2','test');
fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/bar_he.eps'));
%legend('m = 1', 'm = 2', 'm = 3', 'm = 4', 'm = 5', 'm = 10', 'm = 20', 'm = 50')