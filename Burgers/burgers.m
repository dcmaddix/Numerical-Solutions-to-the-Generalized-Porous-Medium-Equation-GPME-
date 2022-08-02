function p = burgers(p0, dt, nt, N, x_cent, x_coord)
    p = p0;
    N_int = length(p(2:end-1));
    i = 2:N_int;
    h = 1/N;
    ind = find(x_cent < x_coord+1e-6 & x_cent > x_coord - 1e-6);
    t = 0:dt:nt;
    p_t = zeros(1,length(t)); %time vector where initial condition is given at t=0
    p_t(1) = p(ind);
    dt
    if (length(ind) > 1)
        ind = ind(1);
    elseif (isempty(ind)) %make sure xcoord exists
        ind = 1;
        fprintf('Ind warning\n')
    end
    lambda=dt/h;
    for n = 2:length(t) %have initial condition time t = 0 => start at time dt
        F=p.^2;%create flux function 1/2u^2
        %p(2:end-1)=0.5*(p(1:end-2)+p(3:end)-lambda*(F(3:end)-F(1:end-2)));
        p(2:end-1) = p(2:end-1) - 0.5*lambda*(F(3:end) - F(1:end-2));
        p_t(n) = p(ind);
%          plot(x_cent, p);
%          pause(1e-6)
    end
    figure(1);
    plot(x_cent, p, '-or');
    xlabel('position x')
    ylabel('pressure p')
    %title(sprintf('Plot of numerical pressure with %s averaging for n = %i versus position', avg, k_pow))
    %legend('pressure', 'Coefficient of p_{xx} in k - truncation error term')
    %saveas(gcf, sprintf('~/Desktop/nonconservform.jpg'));
    hold on;

    h3 = figure(3);
    plot(t,p_t, '-r')
    %ylabel(sprint('pressure at positions x'))
    ylabel(sprintf('pressure at position x = %.1f', x_coord))
    %title(sprintf('Plot of pressure versus time at x = 0.1, 0.2, 0.3, 0.4,0.5, 0.6, respectively'))
    %title(sprintf('Plot of pressure versus time at x = %.1f for N = %i with %s averaging for n = %i', x_coord, N, avg, k_pow))
    xlabel('time t')
    %legend(h3, 'x = 0.1', 'x= 0.2')
    %saveas(gcf, sprintf('~/Desktop/nonconservform_time.jpg'));
    hold on;
end