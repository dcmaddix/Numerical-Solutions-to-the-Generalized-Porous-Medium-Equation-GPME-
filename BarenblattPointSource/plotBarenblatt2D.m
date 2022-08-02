%script to plot 2D barenblatt solution
close all
radial = true;
cross_sect = false;
m = [1 2 3 10 20 50 100];
m = 10;
bound = 7;
x = -bound:.1:bound; %interface with positive region and compact support
y = x;
[x y] = meshgrid(x,y);
t = 1.0;
d = 2; %2D solution
for i = 1:length(m)
    [U,r] = barenblatt_pme(m(i),t,d,x,y,radial); %defined in 
end
if (d == 1)
    plot(x,U, '-o')
elseif (d ==2)
    if (~radial)
        surf(x,y,U)
        shading interp
        colorbar
    elseif(cross_sect)
        plot(x(y == 0), U(y == 0)) %cross-section y = 0
    else
        plot(r,U)
    end
end
if (d == 1)
    legend('m = 1', 'm = 2', 'm = 3', 'm = 10', 'm = 20', 'm = 50', 'm = 100')
    title('Self-similar 1D Barenblatt-Prattle Solution to PME for various m')
    xlabel('x')
    ylabel('p')
elseif(~radial)
    title(sprintf('Self-similar 2D Barenblatt-Prattle Solution to PME for m = %i', m))
    xlabel('x')
    ylabel('y')
    zlabel('p')
else
    title(sprintf('Self-similar Radially Symmetric Barenblatt-Prattle Solution to PME for m = %i', m))
    xlabel('r')
    ylabel('p')
end
figure(2) %stationary solution
surf(x,y,x.^2-y.^2-1)
shading interp
%axis equal