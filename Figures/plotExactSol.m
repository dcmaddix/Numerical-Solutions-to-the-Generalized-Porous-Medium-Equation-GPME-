%compute Exact solution at various times
addpath('HelpfulFunctions/')
t = [0.01, .08, .2, .5];
nolabel = true;
dx = 0.005;
x_cent = 0:dx:1;
lim = true; %false
k_lower = 0; %0.01
close all
str = {'-b', '-r', '-g', '-k'};
p = zeros(length(x_cent), length(t));
for nt = 1:length(t)
    p(:,nt) = computeExactIC(x_cent, lim, k_lower, t(nt));
%     hold on;
%     plot(x_cent, p, str{nt});
%     set(findall(gca, 'Type', 'Line'),'LineWidth',4);
%     set(gca, 'FontSize',16,'FontWeight','bold');
%     title('Solution profile as a function of x')
%     xlabel('x')
%     ylabel('p')
%     leg{nt} = sprintf('t = %.2f', t(nt));
end
createfigure_exactfoam(x_cent,[p 0.5*ones(size(x_cent'))], nolabel);
%plot(x_cent, 0.5*ones(size(x_cent)));
%leg{nt+1} = 'p^* = 0.5';
%legend(leg);
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/Discontinous/Results/';
print(1,'-depsc2','test');
if (nolabel)
    fixPSlinestyle('test.eps',[fileLoc 'exactsol_timeno.eps']);
else
    fixPSlinestyle('test.eps',[fileLoc 'exactsol_time1.eps']);
end