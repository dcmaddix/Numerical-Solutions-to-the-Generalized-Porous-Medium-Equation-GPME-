%convergence study for parabolic problem for k = p^3 and arithmetic
%averaging
N = 100;
n_grids = length(N);
h = 1 ./ N;
l2_error = zeros(3, length(N));
kpow = 3;
dir = sprintf('data/2D/m%i/planar/',kpow);
saveon = true;
for i = 1:n_grids
    exact = load(strcat(dir, sprintf('exact%i.mat',N(i))));
    ph = load(strcat(dir, sprintf('harm%i.mat',N(i))));
    l2_error(1,i) = h(i) * norm(ph.p(:) - exact.p(:));
    pa = load(strcat(dir, sprintf('arith%i.mat',N(i))));
    l2_error(2,i) = h(i) * norm(pa.p(:) - exact.p(:));
    pmhm = load(strcat(dir, sprintf('mhm%i.mat',N(i))));
   l2_error(3,i) = h(i) * norm(pmhm.p(:) - exact.p(:));
end
conv_harm_p3 = polyfit(log10(h), log10(l2_error(1,:)),1);
conv_harm_p3(1)
conv_arith_p3 = polyfit(log10(h), log10(l2_error(2,:)),1);
conv_arith_p3(1)
conv_modharm_p3 = polyfit(log10(h), log10(l2_error(3,:)),1);
conv_modharm_p3(1)
f = figure(1);
plot(log10(h),log10(l2_error(1,:)), '-kd',...,
     log10(h), log10(l2_error(2,:)),'-bs', log10(h), log10(l2_error(3,:)), '-or')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca, 'FontSize',20,'FontWeight','bold');
legend('Harmonic', 'Arithmetic' ,'Modified Harmonic', 'Location', 'southeast')
xlabel('log_{10}(\Deltax)')
ylabel('log_{10}(Error)')
title(sprintf('l_2 norm of error for k = p^%i in 2D',kpow))
if (kpow == 1)
    axis([-2.25 -1.25 -4.5 -1.5])
elseif (kpow == 2)
    axis([-2.25 -1.25 -4 -0.5])
else
    axis([-3 -1.25 -4 0.5])
end
print(1,'-depsc2','test');
if (saveon)
    fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/planar/loglog.eps'));
end
x = 0:h(i):1;
y = x;
[x y] = meshgrid(x,y);
figure(5);
h1 = surf(x,y, abs(pmhm.p - exact.p));
view([80.5 28]);
set(h1,'edgecolor', 'none')
%grid off
colorbar
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',20,'FontWeight','bold');
xlabel('x')
ylabel('y')
title('Absolute error as a function of (x,y)')
zlabel('|p-p_{exact}|')
colormap(paruly)
set(gca, 'zticklabel',[]);
print(5,'-depsc2','test');
if (saveon)
    fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/planar/mhm/error%i.eps', N));
end
figure(8);
h1 = surf(x,y, abs(ph.p - exact.p));
view([80.5 28]);
set(h1,'edgecolor', 'none')
%grid off
colorbar
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',20,'FontWeight','bold');
xlabel('x')
ylabel('y')
zlabel('|p-p_{exact}|')
set(gca, 'zticklabel',[]);
colormap(paruly)
title('Absolute error as a function of (x,y)')
print(8,'-depsc2','test');
if (saveon)
    fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/planar/harm/error%i.eps', N));
end
figure(6);
h1 = surf(x,y, abs(pa.p - exact.p));
view([80.5 28]);
set(h1,'edgecolor', 'none')
%grid off
colorbar
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',20,'FontWeight','bold');
xlabel('x')
ylabel('y')
zlabel('|p-p_{exact}|')
title('Absolute error as a function of (x,y)')
set(gca, 'zticklabel',[]);
colormap(paruly)
print(6,'-depsc2','test');
if (saveon)
    fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/planar/arith/error%i.eps', N));
end
figure(7);
h1 = surf(x,y,ph.p);
view([80.5 28]);
set(h1,'edgecolor', 'none')
%grid off
colorbar
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',20,'FontWeight','bold');
xlabel('x')
ylabel('y')
zlabel('p')
set(gca, 'zticklabel',[]);
colormap(paruly)
title('Solution profile as a function of (x,y)')
print(7,'-depsc2','test');
if (saveon)
    fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/planar/harm/p_final%i.eps', N));
end
figure(1);
h1 = surf(x,y,pmhm.p);
view([80.5 28]);
set(h1,'edgecolor', 'none')
%grid off
colorbar
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',20,'FontWeight','bold');
xlabel('x')
ylabel('y')
zlabel('p')
set(gca, 'zticklabel',[]);
title('Solution profile as a function of (x,y)')
colormap(paruly)
print(1,'-depsc2','test');
if (saveon)
    fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/planar/mhm/p_final%i.eps', N));
end
figure(2);
h1 = surf(x,y,pa.p);
view([80.5 28]);
set(h1,'edgecolor', 'none')
%grid off
set(gca, 'zticklabel',[]);
colorbar
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',20,'FontWeight','bold');
xlabel('x')
ylabel('y')
zlabel('p')
title('Solution profile as a function of (x,y)')
colormap(paruly)
print(2,'-depsc2','test');
if (saveon)
    fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/planar/arith/p_final%i.eps', N));
end
figure(3);
h1 = surf(x,y,exact.p);
view([80.5 28]);
set(h1,'edgecolor', 'none')
%grid off
colorbar
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',20,'FontWeight','bold');
xlabel('x')
ylabel('y')
zlabel('p')
set(gca, 'zticklabel',[]);
title('Solution profile as a function of (x,y)')
colormap(paruly)
print(3,'-depsc2','test');
if (saveon)
    fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/planar/exact/p_final%i.eps', N));
end
f = figure(4);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
plot1 = plot(exact.t, exact.p_t, '-k', pa.t(1:end-1), pa.p_t, pmhm.t(1:end-1), pmhm.p_t, 'LineWidth',4);
set(plot1(2),'LineStyle','-.','Color',[1 0 0],'DisplayName','Arithmetic');
set(plot1(3),'LineStyle',':','Color',[0 0.498039215803146 0],'DisplayName','Mod. Harmonic');
legend('Exact', 'Arithmetic', 'Mod. Harm', 'Location', 'Southeast')
set(gca, 'FontSize',20,'FontWeight','bold');
xlabel('t')
ylabel('p') 
title('Solution profile as a function of t')
print(4,'-depsc2','test');
legend('Arithmetic', 'Exact','MHM')
if (saveon)
    fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/planar/p_t%i.eps', N));
end