kpow = 3;
N = 25;
crossterms2d = false;
timeScheme = 'ForwardEuler';
nolabel = false;
close all
dir = sprintf('data/2D/m%i/',kpow);
if (~strcmp(timeScheme, 'ForwardEuler')) 
    dir = [dir timeScheme '/'];
end
if (kpow == 3)
    Nexact = 400;
else
    Nexact = 200;
end
dtdx = false;
exact = load(sprintf('data/2D/m%i/arith%i.mat',kpow,Nexact));
avg = {'arith', 'harm', 'mhm'};
f = figure(1);
backEuler = false;
for i = 1:3
    if (dtdx)
         num = load(strcat(dir, sprintf('%s%idx.mat', avg{i},N)));
    else
        if (crossterms2d)
            num = load(strcat(dir, sprintf('%s%ict.mat', avg{i},N)));
        else
            num = load(strcat(dir, sprintf('%s%i.mat', avg{i},N)));
        end
    end
    h1 = surf(num.x,num.y,num.p);
    %h1 = surf(num.x,num.y,peaks(num.x,num.y));
    view([0.5 90]);
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca, 'FontSize',20,'FontWeight','bold');
    if (~nolabel)
        xlabel('x')
        ylabel('y')
        zlabel('p')
        colorbar
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca, 'FontSize',20,'FontWeight','bold');
        title('Solution profile as a function of (x,y)')
    end
    if (nolabel)
        set(gca, 'Xtick', [])
        set(gca, 'Ytick', [])
    end
    set(h1,'edgecolor', 'none')
    colormap(paruly)
    axis image
    %grid off
    %shading interp
    print(1,'-depsc2','test');
    if (backEuler)
        fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/p_final%ibackeul.eps',kpow,avg{i},N));
    else
        if (nolabel)
            fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/p_final%ino.eps',kpow,avg{i},N));
        else
            fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/p_final%i.eps',kpow,avg{i},N));
        end
    end
    %saveas(f, sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/p_final%i.jpg',kpow,avg{i},N));
    pow = Nexact / N;
    h1 = surf(num.x,num.y,abs(num.p - exact.p(1:pow:end, 1:pow:end)));
    view([0.5 90]);
    colorbar
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca, 'FontSize',20,'FontWeight','bold');
    xlabel('x')
    ylabel('y')
    zlabel('p')
    axis image
    set(h1,'edgecolor', 'none')
    title('Absolute error as a function of (x,y)')
    %grid off
    %shading interp
%     if (kpow == 1)
%     caxis([0 .012])
% elseif (kpow == 2)
%     if (i ~= 2)
%         caxis([0 .3])
%     else
%         caxis([0 .12])
%     end
% else
%     if (i == 2)
%         caxis([0 .8])
%     elseif (i == 1)
%         caxis([0 1])
%     else
%         caxis([0 0.53])
%     end
% end
    print(1,'-depsc2','test');
    if(~strcmp(timeScheme, 'ForwardEuler'))
        fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/error%i%s.eps',kpow,avg{i},N,timeScheme));
    else
        if (~crossterms2d)   
            fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/error%i.eps',kpow,avg{i},N));
        else
            fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/error%ict.eps',kpow,avg{i},N));
        end
    end
    %saveas(f, sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/error%i.jpg',kpow,avg{i},N));
end
