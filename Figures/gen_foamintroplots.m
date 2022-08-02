close all
addpath('data/foam/')
addpath('IC_0.0479/')
addpath('fixPSlinestyle')
N = 25;
pa = load(sprintf('arith_nx%i.mat', N));
ph = load(sprintf('harm_nx%i.mat', N));
pi = load(sprintf('int_nx%i.mat', N));
addpath('HelpfulFunctions/')
t = 0.0979;
dx = 1/N;
x_cent = 0:dx:1;
lim = true; %false
k_lower = 0; %0.01
p_exact = computeExactIC(x_cent, lim, k_lower, t);
figure(1)
p = plot(pa.x_cent, p_exact, '-k',pa.x_cent,pa.p, '-or', ph.x_cent, ph.p, '-ob', pi.x_cent, pi.p, '-.');
%add in exact solution for reference
set(findall(gca, 'Type', 'Line'),'LineWidth',5);
set(gca, 'FontSize',16,'FontWeight','bold');
set(p(1), 'Linewidth', 9)
set(p(2),'Marker','o','LineStyle',':','Color',[1 0 0],...
    'DisplayName','Arithmetic');
set(p(3),'Marker','^','LineStyle','--','Color',[0 0 1],...
    'DisplayName','Harmonic');
set(p(4),'Marker','x','LineStyle','-.',...
    'Color',[0 0.498039215803146 0],...
    'DisplayName','Integral', 'Markersize', 10);
legend('Exact Sol.', 'Arithmetic', 'Harmonic', 'Integral')
xlabel('x')
ylabel('p')
title('Solution profile as a function of x')
print(1,'-depsc2','test');
fixPSlinestyle('test.eps',[fileLoc 'art.eps']);
p_star = 0.5;
[alpha_star A B] = exact_sol_speed(lim, p_star, k_lower);
figure(2)
x_coord = 0.32;
cell = find(x_cent < x_coord + 1e-6 & x_cent > x_coord - 1e-6, 1, 'first');
p_t = zeros(1, length(pa.t));
for n = 1:length(pa.t)
   p_t(n) = exact_solution(pa.t(n), x_cent, cell, alpha_star,...
                                              A, B, k_lower, lim);
end
t = pa.t - 0.0479;
p = plot(t, p_t, '-k',t,pa.p_t, '-r', t, ph.p_t, '-b', t, pi.p_t, '-.');
%add in exact solution for reference
set(findall(gca, 'Type', 'Line'),'LineWidth',6);
set(gca, 'FontSize',16,'FontWeight','bold');
set(p(2),'LineStyle',':','Color',[1 0 0],...
    'DisplayName','Arithmetic');
set(p(3),'LineStyle','--','Linewidth', 10, 'Color',[0 0 1],...
    'DisplayName','Harmonic');
set(p(4),'LineStyle','-.',...
    'Color',[0 0.498039215803146 0],...
    'DisplayName','Integral');
legend('Exact Sol.', 'Arithmetic', 'Harmonic', 'Integral', 'Location', 'northwest')
xlabel('t')
ylabel('p')
%axes1 = axes('FontWeight','bold','FontSize',16);
annotation('rectangle',...
    [0.391858208955223 0.643945578231293 0.514597014925373 0.14455782312925],...
    'LineWidth',2,...
    'FaceColor','flat');
title('Solution profile as a function of t')
print(2,'-depsc2','test');
fixPSlinestyle('test.eps',[fileLoc 'art_t.eps']);
figure(3)
p = plot(t, p_t, '-k',t,pa.p_t, '-r', t, pi.p_t, '-.');
%add in exact solution for reference
set(findall(gca, 'Type', 'Line'),'LineWidth',6);
set(gca, 'FontSize',16,'FontWeight','bold');
set(p(2),'LineStyle',':','Color',[1 0 0],...
    'DisplayName','Arithmetic');
set(p(3),'LineStyle','-.',...
    'Color',[0 0.498039215803146 0],...
    'DisplayName','Integral');
% Uncomment the following line to preserve the X-limits of the axes
xlim([0.015 0.05]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim([0.45 0.6]);
%box('on');
legend('Exact Sol.', 'Arithmetic', 'Integral', 'Location', 'northwest')
xlabel('t')
ylabel('p')
%axes1 = axes('FontWeight','bold','FontSize',16);
title('Solution profile as a function of t')
print(3,'-depsc2','test');
fixPSlinestyle('test.eps',[fileLoc 'art_tzoom.eps']);