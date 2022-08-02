%convergence study for parabolic problem for k = p^3 and arithmetic
%averaging
N = [25 50 100 200];
n_grids = length(N);
h = 1 ./ N;
l2_error = zeros(3, length(N));
kpow = 1;
if (kpow == 3)
    Nexact = 400;
else
    Nexact = 200;
end
dir_arith = sprintf('data/2D/m%i/arith',kpow);
dir_harm = sprintf('data/2D/m%i/harm',kpow);
dir_modharm = sprintf('data/2D/m%i/mhm',kpow);
exact_arith = load(sprintf('data/2D/m%i/arith%i.mat',kpow,Nexact));
exact_mod = exact_arith;
stch = 1;
if (n_grids >= 4) %4
    stch = 0;
end
saveon = false;
%define x at cell-centers at cell centers to plot against pressure
%exact_arith = load(strcat(dir_arith, 'refsol_p3'));
for i = 1:length(N)
    %store harmonic error
    num = load(strcat(dir_harm, sprintf('%i', N(i))));
    err_harm = num.p - exact_arith.p(1:2^(n_grids+stch-i):end, 1:2^(n_grids+stch-i):end);
    l2_error(1,i) = h(i)*norm(err_harm(:));
    %store arith error
    num = load(strcat(dir_arith, sprintf('%i', N(i))));
    err_arith = num.p - exact_arith.p(1:2^(n_grids+stch-i):end,1:2^(n_grids+stch-i):end);
    l2_error(2,i) = h(i)*norm(err_arith(:));
    %store mhm error
    num = load(strcat(dir_modharm, sprintf('%i', N(i))));
     err_mhm = num.p - exact_mod.p(1:2^(n_grids+stch-i):end,1:2^(n_grids+stch-i):end);
    l2_error(3,i) = h(i)*norm(err_mhm(:));
end
M = length(h);
ind = 1:M;
ind_a = ind;
if (kpow == 2)
    ind = 3:M;
    ind_a = 1:M-1;
elseif (kpow == 3)
    ind = 4:M;
    ind_a = M-2:M-1;
end
conv_harm_p3 = polyfit(log10(h(ind)), log10(l2_error(1,ind)),1);
conv_harm_p3(1)
conv_arith_p3 = polyfit(log10(h(ind_a)), log10(l2_error(2,ind_a)),1);
conv_arith_p3(1)
conv_modharm_p3 = polyfit(log10(h(ind)), log10(l2_error(3,ind)),1);
conv_modharm_p3(1)
f = figure(1);
plot(log10(h),log10(l2_error(1,:)), '-kd',...,
     log10(h), log10(l2_error(2,:)),'-bs', log10(h), log10(l2_error(3,:)), '-or')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca, 'FontSize',16,'FontWeight','bold');
legend('Harmonic', 'Arithmetic' ,'Modified Harmonic', 'Location', 'southeast')
xlabel('log_{10}(\Deltax)')
ylabel('log_{10}(Error)')
title(sprintf('l_2 norm of error for k = p^%i in 2D',kpow))
if (kpow == 1)
    axis([-2.25 -1.25 -4.5 -1.5])
elseif (kpow == 2)
    %axis([-2.25 -1.25 -4 -0.5])
else
    axis([-3 -1.25 -4 0.5])
end
print(1,'-depsc2','test');
if (saveon)
    fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/loglog.eps', kpow));
end
%Harmonic error
i = 3;
dx = 1/N(i);
x = -1:dx:1;
[x, y] = meshgrid(x,x);
f = figure(2);
num = load(strcat(dir_arith, sprintf('%i', N(i))));
err_arith = num.p - exact_arith.p(1:2^(n_grids+stch-i):end, 1:2^(n_grids+stch-i):end);
h =surf(x,y,abs(err_arith));
view([0.5 -90]);
shading interp
xlabel('x')
ylabel('y')
title('Absolute Arithmetic Average Error')
colorbar
if (kpow == 1)
    caxis([0 .012])
elseif (kpow == 2)
    if (i ~= 2)
        caxis([0 .3])
    else
        caxis([0 .12])
    end
else
    if (i == 2)
        caxis([0 .8])
    elseif (i == 1)
        caxis([0 1])
    else
        caxis([0 0.53])
    end
end
print(2,'-depsc2','test');
set(h,'edgecolor','none')
set(h,'LineStyle','none')
if (saveon)
    saveas(f,sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/arith/error%i.jpg', kpow,N(i)));
end
f = figure(4);
num = load(strcat(dir_modharm, sprintf('%i', N(i))));
err_mhm = num.p - exact_arith.p(1:2^(n_grids+stch-i):end, 1:2^(n_grids+stch-i):end);
surf(x,y,abs(err_mhm))
view([0.5 -90]);
shading interp
xlabel('x')
ylabel('y')
title('Absolute MHM Error')
colorbar
if (kpow == 1)
    caxis([0 .012])
elseif (kpow == 2)
    if (i ~= 2)
        caxis([0 .3])
    else
        caxis([0 .12])
    end
else
    if (i == 2)
        caxis([0 .8])
    elseif (i == 1)
        caxis([0 1])
    else
        caxis([0 0.53])
    end
end
grid off
print(4,'-depsc2','test');
if (saveon)
   saveas(f,sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/mhm/error%i.jpg', kpow,N(i)));
end
num = load(strcat(dir_harm, sprintf('%i', N(i))));
err_harm = num.p - exact_arith.p(1:2^(n_grids+stch-i):end, 1:2^(n_grids+stch-i):end);
f = figure(3);
surf(x,y,abs(err_harm))
view([0.5 -90]);
%shading interp
colorbar
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',16,'FontWeight','bold');
xlabel('x')
ylabel('y')
title('Absolute Harmonic Average Error')
if (kpow == 1)
    caxis([0 .012])
elseif (kpow == 2)
    if (i ~= 2)
        caxis([0 .3])
    else
        caxis([0 .12])
    end
else
    if (i == 2)
        caxis([0 .8])
    elseif (i == 1)
        caxis([0 1])
    else
        caxis([0 0.53])
    end
end
axis tight
print(3,'-depsc2','test');
if (saveon)
    saveas(f,sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/harm/error%i.jpg', kpow,N(i)));
end